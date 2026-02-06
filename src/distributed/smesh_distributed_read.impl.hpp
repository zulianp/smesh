#include "smesh_distributed_read.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_path.hpp"
#include "smesh_read.hpp"
#include "smesh_sort.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <assert.h>
#include <limits>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "matrixio_array.h"

#include "smesh_base.hpp"
#include "smesh_distributed_base.hpp"

namespace smesh {

template <typename Index, typename T>
SMESH_INLINE void array_remap_scatter(const ptrdiff_t n,
                                      const Index *SMESH_RESTRICT const mapping,
                                      T *SMESH_RESTRICT const array,
                                      void *SMESH_RESTRICT const temp) {
  T *const temp_actual = static_cast<T *const>(temp);
  memcpy((void *)temp_actual, (const void *)array, (size_t)n * sizeof(T));
  for (ptrdiff_t i = 0; i < n; ++i) {
    array[mapping[i]] = temp_actual[i];
  }
}

template <typename Index, typename T>
SMESH_INLINE void array_remap_gather(const ptrdiff_t n,
                                     const Index *SMESH_RESTRICT const mapping,
                                     T *SMESH_RESTRICT const array,
                                     void *SMESH_RESTRICT const temp) {
  T *const temp_actual = static_cast<T *const>(temp);
  memcpy((void *)temp_actual, (const void *)array, (size_t)n * sizeof(T));
  for (ptrdiff_t i = 0; i < n; ++i) {
    array[i] = temp_actual[mapping[i]];
  }
}

template <typename idx_t>
static SMESH_INLINE int
find_owner_rank(const idx_t idx, const ptrdiff_t n_local_nodes, const int size,
                const idx_t *const SMESH_RESTRICT input_node_partitions) {
  int owner = std::min(size - 1, static_cast<int>(idx / n_local_nodes));

  assert(owner >= 0);
  assert(owner < size);

  if (idx == input_node_partitions[owner]) {
    // Do nothing
  } else if (input_node_partitions[owner + 1] <= idx) {
    while (input_node_partitions[owner + 1] <= idx) {
      owner++;
      assert(owner < size);
    }
  } else if (input_node_partitions[owner] > idx) {
    while (input_node_partitions[owner] > idx) {
      --owner;
      assert(owner >= 0);
    }
  }

  return owner;
}

template <typename FileType, typename TargetType>
int array_create_from_file_convert(MPI_Comm comm, const Path &path,
                                   TargetType **data,
                                   ptrdiff_t *n_local_elements,
                                   ptrdiff_t *n_global_elements) {

  FileType *temp = nullptr;
  if (array_create_from_file(comm, path.c_str(), smesh::mpi_type<FileType>(),
                             (void **)&temp, n_local_elements,
                             n_global_elements) != SMESH_SUCCESS) {
    *data = nullptr;
    *n_local_elements = 0;
    *n_global_elements = 0;
    return SMESH_FAILURE;
  }

  *data = (TargetType *)malloc(*n_local_elements * sizeof(TargetType));
  for (ptrdiff_t i = 0; i < *n_local_elements; i++) {
    (*data)[i] = (TargetType)temp[i];
  }

  free(temp);
  return SMESH_SUCCESS;
}

template <typename T>
int array_create_from_file_convert_from_extension(
    MPI_Comm comm, const Path &path, T **data, ptrdiff_t *n_local_elements,
    ptrdiff_t *n_global_elements) {
  auto ext = path.extension();
  if (ext == ".raw") {
    // We trust the user that the raw file is of the correct type.
    return array_create_from_file(comm, path.c_str(), smesh::mpi_type<T>(),
                                  (void **)data, n_local_elements,
                                  n_global_elements);
  } else if (ext == "float16") {
    return array_create_from_file_convert<f16, T>(
        comm, path, data, n_local_elements, n_global_elements);
  } else if (ext == "float32") {
    return array_create_from_file_convert<f32, T>(
        comm, path, data, n_local_elements, n_global_elements);
  } else if (ext == "float64") {
    return array_create_from_file_convert<f64, T>(
        comm, path, data, n_local_elements, n_global_elements);
  } else if (ext == "int16") {
    return array_create_from_file_convert<i16, T>(
        comm, path, data, n_local_elements, n_global_elements);
  } else if (ext == "int32") {
    return array_create_from_file_convert<i32, T>(
        comm, path, data, n_local_elements, n_global_elements);
  } else if (ext == "int64") {
    return array_create_from_file_convert<i64, T>(
        comm, path, data, n_local_elements, n_global_elements);
  } else {
    SMESH_ERROR("Unsupported file extension %s for file %s\n", ext.c_str(),
                path.c_str());
    return SMESH_FAILURE;
  }
}

template <typename idx_t>
int mesh_build_global_ids(MPI_Comm comm, const ptrdiff_t n_nodes,
                          const ptrdiff_t n_owned_nodes,
                          const ptrdiff_t n_owned_nodes_with_ghosts,
                          idx_t *node_mapping,
                          // int *node_owner,
                          idx_t **node_offsets_out, idx_t **ghosts_out
                          // ,
                          // ptrdiff_t *n_owned_nodes_with_ghosts_out
) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  long n_gnodes = n_owned_nodes;
  long global_node_offset = 0;
  MPI_Exscan(&n_gnodes, &global_node_offset, 1, MPI_LONG, MPI_SUM, comm);

  n_gnodes = global_node_offset + n_owned_nodes;
  MPI_Bcast(&n_gnodes, 1, MPI_LONG, size - 1, comm);

  idx_t *node_offsets = (idx_t *)malloc((size + 1) * sizeof(idx_t));
  SMESH_MPI_CATCH(MPI_Allgather(&global_node_offset, 1,
                                smesh::mpi_type<idx_t>(), node_offsets, 1,
                                smesh::mpi_type<idx_t>(), comm));

  node_offsets[size] = n_gnodes;

  const ptrdiff_t n_lnodes_no_reminder = n_gnodes / size;
  const ptrdiff_t begin = n_lnodes_no_reminder * rank;

  ptrdiff_t n_lnodes_temp = n_lnodes_no_reminder;
  if (rank == size - 1) {
    n_lnodes_temp = n_gnodes - begin;
  }

  const ptrdiff_t begin_owned_with_ghosts =
      n_owned_nodes - n_owned_nodes_with_ghosts;
  const ptrdiff_t extent_owned_with_ghosts = n_owned_nodes_with_ghosts;
  const ptrdiff_t n_ghost_nodes = n_nodes - n_owned_nodes;

  idx_t *ghost_keys = &node_mapping[begin_owned_with_ghosts];
  idx_t *ghost_ids = (idx_t *)malloc(
      std::max(extent_owned_with_ghosts, n_ghost_nodes) * sizeof(idx_t));

  int *send_displs = (int *)malloc((size + 1) * sizeof(int));
  int *recv_displs = (int *)malloc((size + 1) * sizeof(int));
  int *send_count = (int *)calloc(size, sizeof(int));
  int *recv_count = (int *)malloc(size * sizeof(int));

  for (ptrdiff_t i = 0; i < extent_owned_with_ghosts; i++) {
    ghost_ids[i] = global_node_offset + begin_owned_with_ghosts + i;
    assert(ghost_ids[i] < n_gnodes);
  }

  for (ptrdiff_t i = 0; i < extent_owned_with_ghosts; i++) {
    const idx_t idx = ghost_keys[i];
    int dest_rank =
        std::min(static_cast<int>(idx / n_lnodes_no_reminder), size - 1);
    assert(dest_rank < size);
    assert(dest_rank >= 0);

    send_count[dest_rank]++;
  }

  SMESH_MPI_CATCH(
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm));

  send_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    send_displs[r + 1] = send_displs[r] + send_count[r];
  }

  recv_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  idx_t *recv_key_buff = (idx_t *)malloc(recv_displs[size] * sizeof(idx_t));

  SMESH_MPI_CATCH(MPI_Alltoallv(
      ghost_keys, send_count, send_displs, smesh::mpi_type<idx_t>(),
      recv_key_buff, recv_count, recv_displs, smesh::mpi_type<idx_t>(), comm));

  idx_t *recv_ids_buff = (idx_t *)malloc(recv_displs[size] * sizeof(idx_t));

  SMESH_MPI_CATCH(MPI_Alltoallv(
      ghost_ids, send_count, send_displs, smesh::mpi_type<idx_t>(),
      recv_ids_buff, recv_count, recv_displs, smesh::mpi_type<idx_t>(), comm));

  idx_t *mapping = (idx_t *)malloc(n_lnodes_temp * sizeof(idx_t));

#ifndef NDEBUG
  for (ptrdiff_t i = 0; i < n_lnodes_temp; i++) {
    mapping[i] = smesh::invalid_idx<idx_t>();
  }
#endif

  // Fill mapping
  for (int r = 0; r < size; r++) {
    int proc_begin = recv_displs[r];
    int proc_extent = recv_count[r];

    idx_t *keys = &recv_key_buff[proc_begin];
    idx_t *ids = &recv_ids_buff[proc_begin];

    for (int k = 0; k < proc_extent; k++) {
      idx_t iii = keys[k] - begin;

      assert(iii >= 0);
      assert(iii < n_lnodes_temp);

      mapping[iii] = ids[k];
    }
  }

  /////////////////////////////////////////////////
  // Gather query ghost nodes
  memset(send_count, 0, size * sizeof(int));

  // Get the query nodes
  ghost_keys = &node_mapping[n_owned_nodes];

  idx_t *recv_idx = (idx_t *)malloc(n_ghost_nodes * sizeof(idx_t));
  idx_t *exchange_buff = (idx_t *)malloc(n_ghost_nodes * sizeof(idx_t));
  {
    for (ptrdiff_t i = 0; i < n_ghost_nodes; i++) {
      const idx_t idx = ghost_keys[i];
      int dest_rank =
          std::min(static_cast<int>(idx / n_lnodes_no_reminder), size - 1);

      assert(dest_rank < size);
      assert(dest_rank >= 0);

      send_count[dest_rank]++;
    }

    send_displs[0] = 0;
    for (int r = 0; r < size; r++) {
      send_displs[r + 1] = send_displs[r] + send_count[r];
    }

    memset(send_count, 0, sizeof(int) * size);

    for (ptrdiff_t i = 0; i < n_ghost_nodes; i++) {
      const idx_t idx = ghost_keys[i];
      int dest_rank =
          std::min(static_cast<int>(idx / n_lnodes_no_reminder), size - 1);

      assert(dest_rank < size);
      assert(dest_rank >= 0);

      const idx_t offset = send_displs[dest_rank] + send_count[dest_rank]++;
      exchange_buff[offset] = idx;
      recv_idx[offset] = i;
    }
  }

  SMESH_MPI_CATCH(
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm));

  recv_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  recv_key_buff =
      (idx_t *)realloc(recv_key_buff, recv_displs[size] * sizeof(idx_t));

  SMESH_MPI_CATCH(MPI_Alltoallv(
      exchange_buff, send_count, send_displs, smesh::mpi_type<idx_t>(),
      recv_key_buff, recv_count, recv_displs, smesh::mpi_type<idx_t>(), comm));

  // Query mapping
  for (int r = 0; r < size; r++) {
    int proc_begin = recv_displs[r];
    int proc_extent = recv_count[r];

    idx_t *keys = &recv_key_buff[proc_begin];

    for (int k = 0; k < proc_extent; k++) {
      idx_t iii = keys[k] - begin;

      if (iii >= n_lnodes_temp) {
        printf("[%d] %ld < %d < %ld\n", rank, begin, keys[k],
               begin + n_lnodes_temp);
      }

      assert(iii < n_lnodes_temp);
      assert(iii >= 0);
      assert(mapping[iii] >= 0);
      keys[k] = mapping[iii];
    }
  }

  // Send back
  SMESH_MPI_CATCH(MPI_Alltoallv(
      recv_key_buff, recv_count, recv_displs, smesh::mpi_type<idx_t>(),
      exchange_buff, send_count, send_displs, smesh::mpi_type<idx_t>(), comm));

  idx_t *ghosts = (idx_t *)malloc(n_ghost_nodes * sizeof(idx_t));
  for (ptrdiff_t i = 0; i < n_ghost_nodes; i++) {
    ghosts[recv_idx[i]] = exchange_buff[i];
  }

  free(recv_idx);
  free(exchange_buff);

  /////////////////////////////////////////////////

  *node_offsets_out = node_offsets;
  *ghosts_out = ghosts;

  // Clean-up
  free(send_count);
  free(recv_count);
  free(send_displs);
  free(recv_displs);
  free(recv_key_buff);
  free(recv_ids_buff);
  free(mapping);
  free(ghost_ids);

  return SMESH_SUCCESS;
}

// ------------------

template <typename idx_t>
int read_mapped_field(MPI_Comm comm, const char *input_path,
                      const ptrdiff_t n_local, const ptrdiff_t n_global,
                      const idx_t *SMESH_RESTRICT const mapping,
                      MPI_Datatype data_type,
                      void *SMESH_RESTRICT const data_out) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  uint8_t *const out = (uint8_t *const)data_out;

  int type_size = 0;
  SMESH_MPI_CATCH(MPI_Type_size(data_type, &type_size));

  const ptrdiff_t local_size_no_remainder = n_global / size;
  const ptrdiff_t begin = (n_global / size) * rank;

  ptrdiff_t local_size = local_size_no_remainder;
  if (rank == size - 1) {
    local_size = n_global - begin;
  }

  // Read this rank's contiguous chunk of the global field
  uint8_t *local_chunk =
      (uint8_t *)malloc((size_t)local_size * (size_t)type_size);
  if (!local_chunk)
    return SMESH_FAILURE;

  int err = array_read(comm, input_path, data_type, (void *)local_chunk,
                       local_size, n_global);
  if (err) {
    free(local_chunk);
    return err;
  }

  // Build request counts: how many global indices we need from each rank
  int *req_count = (int *)calloc((size_t)size, sizeof(int));
  if (!req_count) {
    free(local_chunk);
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];

    int src_rank = 0;
    if (local_size_no_remainder > 0) {
      src_rank = std::min(size - 1, (int)(idx / local_size_no_remainder));
    } else {
      // Degenerate distribution (n_global < size): everything lives on last
      // rank.
      src_rank = size - 1;
    }

    req_count[src_rank]++;
  }

  int *recv_req_count = (int *)malloc((size_t)size * sizeof(int));
  if (!recv_req_count) {
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  SMESH_MPI_CATCH(MPI_Alltoall(req_count, 1, smesh::mpi_type<int>(),
                               recv_req_count, 1, smesh::mpi_type<int>(),
                               comm));

  int *req_displs = (int *)malloc((size_t)size * sizeof(int));
  int *recv_displs = (int *)malloc((size_t)size * sizeof(int));
  if (!req_displs || !recv_displs) {
    free(req_displs);
    free(recv_displs);
    free(recv_req_count);
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  req_displs[0] = 0;
  recv_displs[0] = 0;
  for (int i = 0; i < size - 1; ++i) {
    req_displs[i + 1] = req_displs[i] + req_count[i];
    recv_displs[i + 1] = recv_displs[i] + recv_req_count[i];
  }

  const ptrdiff_t total_recv_req =
      (ptrdiff_t)recv_displs[size - 1] + (ptrdiff_t)recv_req_count[size - 1];

  // Pack requests: global indices + local positions (for unpack)
  idx_t *req_list = (idx_t *)malloc((size_t)n_local * sizeof(idx_t));
  ptrdiff_t *local_pos =
      (ptrdiff_t *)malloc((size_t)n_local * sizeof(ptrdiff_t));
  ptrdiff_t *book_keeping =
      (ptrdiff_t *)calloc((size_t)size, sizeof(ptrdiff_t));
  if (!req_list || !local_pos || !book_keeping) {
    free(book_keeping);
    free(local_pos);
    free(req_list);
    free(recv_displs);
    free(req_displs);
    free(recv_req_count);
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];

    int src_rank = 0;
    if (local_size_no_remainder > 0) {
      src_rank = std::min(size - 1, (int)(idx / local_size_no_remainder));
    } else {
      src_rank = size - 1;
    }

    const ptrdiff_t off =
        (ptrdiff_t)req_displs[src_rank] + (ptrdiff_t)book_keeping[src_rank];
    req_list[off] = idx;
    local_pos[off] = i;
    book_keeping[src_rank]++;
  }

  idx_t *recv_req_list =
      (idx_t *)malloc((size_t)total_recv_req * sizeof(idx_t));
  if (!recv_req_list) {
    free(book_keeping);
    free(local_pos);
    free(req_list);
    free(recv_displs);
    free(req_displs);
    free(recv_req_count);
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  // Exchange requested indices
  SMESH_MPI_CATCH(MPI_Alltoallv(
      req_list, req_count, req_displs, smesh::mpi_type<idx_t>(), recv_req_list,
      recv_req_count, recv_displs, smesh::mpi_type<idx_t>(), comm));

  // Build response buffer for received requests (same ordering as
  // recv_req_list)
  uint8_t *send_resp =
      (uint8_t *)malloc((size_t)total_recv_req * (size_t)type_size);
  if (!send_resp) {
    free(recv_req_list);
    free(book_keeping);
    free(local_pos);
    free(req_list);
    free(recv_displs);
    free(req_displs);
    free(recv_req_count);
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < total_recv_req; ++i) {
    const idx_t idx = recv_req_list[i];
    const ptrdiff_t loc = (ptrdiff_t)idx - begin;
    assert(loc >= 0);
    assert(loc < local_size);
    memcpy((void *)(send_resp + i * type_size),
           (const void *)(local_chunk + loc * type_size), (size_t)type_size);
  }

  // Exchange response data back to requesters
  uint8_t *recv_resp = (uint8_t *)malloc((size_t)n_local * (size_t)type_size);
  if (!recv_resp) {
    free(send_resp);
    free(recv_req_list);
    free(book_keeping);
    free(local_pos);
    free(req_list);
    free(recv_displs);
    free(req_displs);
    free(recv_req_count);
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  SMESH_MPI_CATCH(MPI_Alltoallv(send_resp, recv_req_count, recv_displs,
                                data_type, recv_resp, req_count, req_displs,
                                data_type, comm));

  // Unpack into local ordering
  for (ptrdiff_t off = 0; off < n_local; ++off) {
    const ptrdiff_t i = local_pos[off];
    memcpy((void *)(out + i * type_size),
           (const void *)(recv_resp + off * type_size), (size_t)type_size);
  }

  free(recv_resp);
  free(send_resp);
  free(recv_req_list);
  free(book_keeping);
  free(local_pos);
  free(req_list);
  free(recv_displs);
  free(req_displs);
  free(recv_req_count);
  free(req_count);
  free(local_chunk);
  return SMESH_SUCCESS;
}

template <typename idx_t>
int mesh_block_from_folder(MPI_Comm comm, const Path &folder,
                           int *nnodesxelem_out, idx_t ***const elems,
                           ptrdiff_t *const n_local_elements_out,
                           ptrdiff_t *const n_global_elements_out) {

  std::vector<Path> i_files =
      detect_files(folder / "i*.*", {".raw", ".int16", ".int32", ".int64"});

  int nnodesxelem = i_files.size();
  *elems = (idx_t **)malloc(sizeof(idx_t *) * nnodesxelem);
  for (int d = 0; d < nnodesxelem; ++d) {
    (*elems)[d] = nullptr;
  }

  ptrdiff_t n_local_elements = 0;
  ptrdiff_t n_global_elements = 0;
  int ret = SMESH_SUCCESS;
  {
    ptrdiff_t n_elements0 = 0;
    ptrdiff_t n_global_elements0 = 0;
    for (int d = 0; d < nnodesxelem; ++d) {
      Path i_path = i_files[d];
      std::string filename = i_path.file_name();
      int ii = std::stoi(filename.substr(1, filename.find_last_of('.')));

      idx_t *idx = 0;
      if (array_create_from_file_convert_from_extension<idx_t>(
              comm, i_path, &idx, &n_local_elements, &n_global_elements) !=
          SMESH_SUCCESS) {
        SMESH_ERROR("Failed to read index file %s\n", i_path.c_str());
        ret = SMESH_FAILURE;
      }
      // End of Selection
      (*elems)[ii] = idx;

      if (d == 0) {
        n_elements0 = n_local_elements;
        n_global_elements0 = n_global_elements;
      } else {
        assert(n_elements0 == n_local_elements);
        assert(n_global_elements0 == n_global_elements);

        if (n_elements0 != n_local_elements ||
            n_global_elements0 != n_global_elements) {
          SMESH_ERROR(
              "Inconsistent lenghts in input %ld != %ld or %ld != %ld\n",
              (long)n_elements0, (long)n_local_elements,
              (long)n_global_elements0, (long)n_global_elements);
          ret = SMESH_FAILURE;
          break;
        }
      }
    }
  }

  if (ret == SMESH_FAILURE) {
    *nnodesxelem_out = 0;
    *n_local_elements_out = 0;
    *n_global_elements_out = 0;
    *elems = nullptr;
    for (int d = 0; d < nnodesxelem; ++d) {
      free((*elems)[d]);
    }
    free(*elems);
    *elems = nullptr;
    return SMESH_FAILURE;
  } else {
    *nnodesxelem_out = nnodesxelem;
    *n_local_elements_out = n_local_elements;
    *n_global_elements_out = n_global_elements;
    return SMESH_SUCCESS;
  }
}

template <typename geom_t>
int mesh_coordinates_from_folder(MPI_Comm comm, const Path &folder,
                                 int *spatial_dim_out, geom_t ***points_out,
                                 ptrdiff_t *n_local_nodes_out,
                                 ptrdiff_t *n_global_nodes_out) {
  std::vector<Path> x_file = detect_files(
      folder / "x.*", {".raw", ".float16", ".float32", ".float64"});
  std::vector<Path> y_file = detect_files(
      folder / "y.*", {".raw", ".float16", ".float32", ".float64"});
  std::vector<Path> z_file = detect_files(
      folder / "z.*", {".raw", ".float16", ".float32", ".float64"});

  if (x_file.empty()) {
    x_file = detect_files(folder / "x0.*",
                          {".raw", ".float16", ".float32", ".float64"});
  }

  if (y_file.empty()) {
    y_file = detect_files(folder / "x1.*",
                          {".raw", ".float16", ".float32", ".float64"});
  }

  if (z_file.empty()) {
    z_file = detect_files(folder / "x2.*",
                          {".raw", ".float16", ".float32", ".float64"});
  }

  int ndims = x_file.empty() ? 0 : 1; // x only
  ndims += y_file.empty() ? 0 : 1;    // x and y
  ndims += z_file.empty() ? 0 : 1;    // x, y and z

  std::vector<Path> points_paths = x_file;
  if (!y_file.empty()) {
    points_paths.push_back(y_file[0]);
  }
  if (!z_file.empty()) {
    points_paths.push_back(z_file[0]);
  }

  geom_t **points = (geom_t **)malloc(sizeof(geom_t *) * ndims);

  ptrdiff_t n_local_nodes = 0;
  ptrdiff_t n_global_nodes = 0;
  ptrdiff_t n_local_nodes0 = 0;
  ptrdiff_t n_global_nodes0 = 0;
  int ret = SMESH_SUCCESS;
  for (int d = 0; d < ndims; ++d) {
    geom_t *points_d = 0;
    if (array_create_from_file_convert_from_extension(
            comm, points_paths[d], &points_d, &n_local_nodes,
            &n_global_nodes) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
    points[d] = points_d;

    if (d == 0) {
      n_local_nodes0 = n_local_nodes;
      n_global_nodes0 = n_global_nodes;
    } else {
      assert(n_local_nodes0 == n_local_nodes);
      assert(n_global_nodes0 == n_global_nodes);
    }

    if (n_local_nodes0 != n_local_nodes || n_global_nodes0 != n_global_nodes) {
      SMESH_ERROR("Inconsistent lenghts in input %ld != %ld or %ld != %ld\n",
                  (long)n_local_nodes0, (long)n_local_nodes,
                  (long)n_global_nodes0, (long)n_global_nodes);
      ret = SMESH_FAILURE;
      break;
    }
  }

  if (ret == SMESH_FAILURE) {
    *spatial_dim_out = 0;
    *points_out = nullptr;
    *n_local_nodes_out = 0;
    *n_global_nodes_out = 0;
    return SMESH_FAILURE;
  } else {
    *spatial_dim_out = ndims;
    *points_out = points;
    *n_local_nodes_out = n_local_nodes;
    *n_global_nodes_out = n_global_nodes;
    return SMESH_SUCCESS;
  }
}

template <typename idx_t, typename geom_t, typename element_idx_t>
int mesh_from_folder(
    const MPI_Comm comm, const Path &folder, int *nnodesxelem_out,
    ptrdiff_t *nelements_out, idx_t ***elements_out, int *spatial_dim_out,
    ptrdiff_t *nnodes_out, geom_t ***points_out, ptrdiff_t *n_owned_nodes_out,
    ptrdiff_t *n_owned_elements_out, element_idx_t **element_mapping_out,
    idx_t **node_mapping_out, int **node_owner_out, idx_t **node_offsets_out,
    idx_t **ghosts_out, ptrdiff_t *n_owned_nodes_with_ghosts_out,
    ptrdiff_t *n_shared_elements_out,
    ptrdiff_t *n_owned_elements_with_ghosts_out) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  static const int remap_elements = 1;

  if (size > 1) {
    double tick = MPI_Wtime();

    ////////////////////////////////////////////////////////////////////////////////
    // Read elements
    ////////////////////////////////////////////////////////////////////////////////

    int nnodesxelem = 0;
    ptrdiff_t n_local_elements = 0, n_elements = 0;
    idx_t **elems = nullptr;
    mesh_block_from_folder(comm, folder, &nnodesxelem, &elems,
                           &n_local_elements, &n_elements);

    idx_t *unique_idx =
        (idx_t *)malloc(sizeof(idx_t) * n_local_elements * nnodesxelem);
    for (int d = 0; d < nnodesxelem; ++d) {
      memcpy(&unique_idx[d * n_local_elements], elems[d],
             sizeof(idx_t) * n_local_elements);
    }

    ptrdiff_t n_unique =
        sort_and_unique<idx_t>(unique_idx, n_local_elements * nnodesxelem);

    int ndims = 0;
    ptrdiff_t n_local_nodes = 0, n_nodes = 0;
    geom_t **xyz = nullptr;
    mesh_coordinates_from_folder(comm, folder, &ndims, &xyz, &n_local_nodes,
                                 &n_nodes);

    ////////////////////////////////////////////////////////////////////////////////

    idx_t *input_node_partitions = (idx_t *)malloc(sizeof(idx_t) * (size + 1));
    memset(input_node_partitions, 0, sizeof(idx_t) * (size + 1));
    input_node_partitions[rank + 1] = n_local_nodes;

    SMESH_MPI_CATCH(MPI_Allreduce(MPI_IN_PLACE, &input_node_partitions[1], size,
                                  smesh::mpi_type<idx_t>(), MPI_SUM, comm));

    for (int r = 0; r < size; ++r) {
      input_node_partitions[r + 1] += input_node_partitions[r];
    }

    int *gather_node_count = (int *)malloc(size * sizeof(int));
    memset(gather_node_count, 0, size * sizeof(int));

    for (ptrdiff_t i = 0; i < n_unique; ++i) {
      idx_t idx = unique_idx[i];
      const int owner =
          find_owner_rank(idx, n_local_nodes, size, input_node_partitions);

      assert(owner < size);
      assert(owner >= 0);

      assert(idx >= input_node_partitions[owner]);
      assert(idx < input_node_partitions[owner + 1]);

      gather_node_count[owner]++;
    }

    int *scatter_node_count = (int *)malloc(size * sizeof(int));
    memset(scatter_node_count, 0, size * sizeof(int));

    SMESH_MPI_CATCH(MPI_Alltoall(gather_node_count, 1, smesh::mpi_type<idx_t>(),
                                 scatter_node_count, 1,
                                 smesh::mpi_type<idx_t>(), comm));

    int *gather_node_displs = (int *)malloc((size + 1) * sizeof(int));
    int *scatter_node_displs = (int *)malloc((size + 1) * sizeof(int));

    gather_node_displs[0] = 0;
    scatter_node_displs[0] = 0;

    for (int i = 0; i < size; i++) {
      gather_node_displs[i + 1] = gather_node_displs[i] + gather_node_count[i];
    }

    for (int i = 0; i < size; i++) {
      scatter_node_displs[i + 1] =
          scatter_node_displs[i] + scatter_node_count[i];
    }

    ptrdiff_t size_send_list = scatter_node_displs[size];

    idx_t *send_list = (idx_t *)malloc(sizeof(idx_t) * size_send_list);
    memset(send_list, 0, sizeof(idx_t) * size_send_list);

    SMESH_MPI_CATCH(
        MPI_Alltoallv(unique_idx, gather_node_count, gather_node_displs,
                      smesh::mpi_type<idx_t>(), send_list, scatter_node_count,
                      scatter_node_displs, smesh::mpi_type<idx_t>(), comm));

    ///////////////////////////////////////////////////////////////////////

    // Remove offset
    for (ptrdiff_t i = 0; i < size_send_list; ++i) {
      send_list[i] -= input_node_partitions[rank];
    }

    ///////////////////////////////////////////////////////////////////////
    // Exchange points

    geom_t *sendx = (geom_t *)malloc(size_send_list * sizeof(geom_t));
    geom_t **part_xyz = (geom_t **)malloc(sizeof(geom_t *) * ndims);

    for (int d = 0; d < ndims; ++d) {
      // Fill buffer
      for (ptrdiff_t i = 0; i < size_send_list; ++i) {
        sendx[i] = xyz[d][send_list[i]];
      }

      geom_t *recvx = (geom_t *)malloc(n_unique * sizeof(geom_t));
      SMESH_MPI_CATCH(
          MPI_Alltoallv(sendx, scatter_node_count, scatter_node_displs,
                        smesh::mpi_type<geom_t>(), recvx, gather_node_count,
                        gather_node_displs, smesh::mpi_type<geom_t>(), comm));
      part_xyz[d] = recvx;

      // Free space
      free(xyz[d]);
    }

    free(xyz);

    ///////////////////////////////////////////////////////////////////////
    // Determine owners
    int *node_owner = (int *)malloc(n_unique * sizeof(int));
    int *node_share_count = (int *)calloc(n_unique, sizeof(int));

    {
      int *decide_node_owner = (int *)malloc(n_local_nodes * sizeof(int));
      int *send_node_owner = (int *)malloc(size_send_list * sizeof(int));

      int *decide_share_count = (int *)calloc(n_local_nodes, sizeof(int));
      int *send_share_count = (int *)malloc(size_send_list * sizeof(int));

      for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
        decide_node_owner[i] = size;
      }

      for (int r = 0; r < size; ++r) {
        int begin = scatter_node_displs[r];
        int end = scatter_node_displs[r + 1];

        for (int i = begin; i < end; ++i) {
          decide_node_owner[send_list[i]] =
              std::min(decide_node_owner[send_list[i]], r);
        }

        for (int i = begin; i < end; ++i) {
          decide_share_count[send_list[i]]++;
        }
      }

      for (int r = 0; r < size; ++r) {
        int begin = scatter_node_displs[r];
        int end = scatter_node_displs[r + 1];

        for (int i = begin; i < end; ++i) {
          send_node_owner[i] = decide_node_owner[send_list[i]];
        }

        for (int i = begin; i < end; ++i) {
          send_share_count[i] = decide_share_count[send_list[i]];
        }
      }

      SMESH_MPI_CATCH(MPI_Alltoallv(
          send_node_owner, scatter_node_count, scatter_node_displs, MPI_INT,
          node_owner, gather_node_count, gather_node_displs, MPI_INT, comm));

      SMESH_MPI_CATCH(MPI_Alltoallv(send_share_count, scatter_node_count,
                                    scatter_node_displs, MPI_INT,
                                    node_share_count, gather_node_count,
                                    gather_node_displs, MPI_INT, comm));

      free(decide_node_owner);
      free(send_node_owner);

      free(decide_share_count);
      free(send_share_count);
    }

    ///////////////////////////////////////////////////////////////////////
    // Localize element index
    for (ptrdiff_t d = 0; d < nnodesxelem; ++d) {
      for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
        auto low =
            std::lower_bound(unique_idx, unique_idx + n_unique, elems[d][e]);
        elems[d][e] = std::distance(unique_idx, low);
      }
    }

    ////////////////////////////////////////////////////////////
    // Remap node index
    // We reorder nodes with the following order
    // 1) locally owned
    // 2) locally owned and shared by a remote process
    // 3) Shared, hence owned by a remote process

    idx_t *proc_ptr = (idx_t *)calloc((size + 1), sizeof(idx_t));
    idx_t *offset = (idx_t *)calloc((size), sizeof(idx_t));

    for (ptrdiff_t node = 0; node < n_unique; ++node) {
      proc_ptr[node_owner[node] + 1] += 1;
    }

    const ptrdiff_t n_owned_nodes = proc_ptr[rank + 1];
    // Remove offset
    proc_ptr[rank + 1] = 0;

    proc_ptr[0] = n_owned_nodes;
    for (int r = 0; r < size; ++r) {
      proc_ptr[r + 1] += proc_ptr[r];
    }

    // This rank comes first
    proc_ptr[rank] = 0;

    // Build local remap index
    idx_t *local_remap = (idx_t *)malloc((n_unique) * sizeof(idx_t));
    for (ptrdiff_t node = 0; node < n_unique; ++node) {
      int owner = node_owner[node];
      local_remap[node] = proc_ptr[owner] + offset[owner]++;
    }

    if (1) {
      // Remap based on shared
      ptrdiff_t owned_shared_count[2] = {0, 0};

      for (ptrdiff_t node = 0; node < n_unique; ++node) {
        int owner = node_owner[node];
        if (owner != rank)
          continue;
        owned_shared_count[node_share_count[node] > 1]++;
      }

      ptrdiff_t owned_shared_offset[2] = {0, owned_shared_count[0]};
      for (ptrdiff_t node = 0; node < n_unique; ++node) {
        int owner = node_owner[node];
        if (owner != rank)
          continue;

        int owned_and_shared = node_share_count[node] > 1;
        local_remap[node] =
            proc_ptr[owner] + owned_shared_offset[owned_and_shared]++;
      }

      *n_owned_nodes_with_ghosts_out = owned_shared_count[1];
    }

    const size_t max_sz =
        std::max(sizeof(int), std::max(sizeof(idx_t), sizeof(geom_t)));
    void *temp_buff = malloc(n_unique * max_sz);

    for (int d = 0; d < ndims; ++d) {
      geom_t *x = part_xyz[d];
      array_remap_scatter(n_unique, local_remap, x, temp_buff);
    }

    array_remap_scatter(n_unique, local_remap, node_owner, temp_buff);
    array_remap_scatter(n_unique, local_remap, unique_idx, temp_buff);
    array_remap_scatter(n_unique, local_remap, node_share_count, temp_buff);

    for (int d = 0; d < nnodesxelem; ++d) {
      for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
        elems[d][e] = local_remap[elems[d][e]];
      }
    }

    free(local_remap);
    free(proc_ptr);
    free(temp_buff);
    free(offset);

    ////////////////////////////////////////////////////////////
    // Remap element index

    // Reorder elements with the following order
    // 1) Locally owned element
    // 2) Elements that are locally owned but have node shared by a remote
    // process

    if (remap_elements) {
      idx_t *temp_buff = (idx_t *)malloc(n_local_elements * sizeof(idx_t));
      uint8_t *is_local = (uint8_t *)calloc(n_local_elements, sizeof(uint8_t));

      for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
        int is_local_e = 1;
        int is_owned_and_shared = 0;

        for (int d = 0; d < nnodesxelem; ++d) {
          const idx_t idx = elems[d][e];

          if (node_owner[idx] != rank) {
            is_local_e = 0;
          }

          is_owned_and_shared += node_share_count[idx] > 1;
        }

        is_local[e] = is_local_e;
        if (is_local_e && is_owned_and_shared) {
          is_local[e] += 1;
        }
      }

      // FIXME?
      idx_t *element_mapping =
          (idx_t *)malloc(n_local_elements * sizeof(idx_t));

      ptrdiff_t counter = 0;
      for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
        if (is_local[e] == 1) {
          element_mapping[counter++] = e;
        }
      }

      *n_owned_elements_with_ghosts_out = 0;
      for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
        if (is_local[e] == 2) {
          element_mapping[counter++] = e;
          (*n_owned_elements_with_ghosts_out)++;
        }
      }

      *n_shared_elements_out = n_local_elements - counter;
      for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
        if (!is_local[e]) {
          element_mapping[counter++] = e;
        }
      }

      for (int d = 0; d < nnodesxelem; ++d) {
        array_remap_gather(n_local_elements, element_mapping, elems[d],
                                  temp_buff);
      }

      *element_mapping_out = element_mapping;

      free(temp_buff);
      free(is_local);
    } else {
      *element_mapping_out = 0;
      *n_shared_elements_out = 0;
      *n_owned_elements_with_ghosts_out = 0;
    }

    ///////////////////////////////////////////////////////////////////////
    // Free space
    free(sendx);
    free(send_list);
    free(input_node_partitions);
    free(node_share_count);

    free(scatter_node_count);
    free(gather_node_count);
    free(gather_node_displs);
    free(scatter_node_displs);

    ///////////////////////////////////////////////////////////////////////
    *spatial_dim_out = ndims;
    *nnodesxelem_out = nnodesxelem;

    *nelements_out = n_local_elements;
    *nnodes_out = n_unique;

    *elements_out = elems;
    *points_out = part_xyz;

    *n_owned_nodes_out = n_owned_nodes;
    *n_owned_elements_out = n_local_elements - *n_shared_elements_out;

    *node_mapping_out = unique_idx;
    *node_owner_out = node_owner;

    *ghosts_out = 0;
    *node_offsets_out = 0;

    mesh_build_global_ids(comm, n_nodes, n_local_nodes, n_local_nodes,
                          *node_mapping_out,
                          /**node_owner_out,*/ node_offsets_out,
                          ghosts_out /*,n_owned_nodes_with_ghosts_out*/);

    // MPI_Barrier(comm);
    double tock = MPI_Wtime();
    if (!rank) {
      printf("read_mesh.c: read_mesh\t%g seconds\n", tock - tick);
    }

    return SMESH_SUCCESS;
  } else {
    SMESH_ERROR("Serial mesh reading not implemented\n");
    return SMESH_FAILURE;
  }
}

} // namespace smesh