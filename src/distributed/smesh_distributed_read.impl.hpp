#include "smesh_decompose.hpp"
#include "smesh_distributed_aura.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_graph.hpp"
#include "smesh_path.hpp"
#include "smesh_read.hpp"
#include "smesh_sort.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <assert.h>
#include <chrono>
#include <fstream>
#include <limits>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "matrixio_array.h"

#include "smesh_alltoallv.impl.hpp"
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
  // SMESH_TRACE_SCOPE("array_create_from_file_convert");

  FileType *temp = nullptr;
  if (array_create_from_file(comm, path.c_str(), smesh::mpi_type<FileType>(),
                             (void **)&temp, n_local_elements,
                             n_global_elements) != SMESH_SUCCESS) {
    *data = nullptr;
    *n_local_elements = 0;
    *n_global_elements = 0;
    SMESH_ERROR("Failed to create array from file %s\n", path.c_str());
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
  // SMESH_TRACE_SCOPE("array_create_from_file_convert_from_extension");
  auto ext = path.extension();
  if (ext == "raw") {
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
  i64 *req_count = (i64 *)calloc((size_t)size, sizeof(i64));
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

  i64 *recv_req_count = (i64 *)malloc((size_t)size * sizeof(i64));
  if (!recv_req_count) {
    free(req_count);
    free(local_chunk);
    return SMESH_FAILURE;
  }

  SMESH_MPI_CATCH(MPI_Alltoall(req_count, 1, smesh::mpi_type<i64>(),
                               recv_req_count, 1, smesh::mpi_type<i64>(),
                               comm));

  i64 *req_displs = (i64 *)malloc(((size_t)size + 1) * sizeof(i64));
  i64 *recv_displs = (i64 *)malloc(((size_t)size + 1) * sizeof(i64));
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
  for (int i = 0; i < size; ++i) {
    req_displs[i + 1] = req_displs[i] + req_count[i];
    recv_displs[i + 1] = recv_displs[i] + recv_req_count[i];
  }

  const i64 total_recv_req = recv_displs[size];

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
  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64(req_list, req_count, req_displs, recv_req_list,
                                 recv_req_count, recv_displs, comm,
                                 max_chunk_size));

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

  for (i64 i = 0; i < total_recv_req; ++i) {
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

  SMESH_MPI_CATCH(all_to_allv_64_b(send_resp, recv_req_count, recv_displs,
                                   data_type, recv_resp, req_count, req_displs,
                                   data_type, comm, max_chunk_size));

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
  // SMESH_TRACE_SCOPE("mesh_block_from_folder");

  std::vector<Path> i_files =
      detect_files(folder / "i*.*", {"raw", "int16", "int32", "int64"});

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

      if (ii >= nnodesxelem) {
        SMESH_ERROR("Index out of range: %d >= %d\n", ii, nnodesxelem);
        ret = SMESH_FAILURE;
        break;
      }

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
  std::vector<Path> x_file =
      detect_files(folder / "x.*", {"raw", "float16", "float32", "float64"});
  std::vector<Path> y_file =
      detect_files(folder / "y.*", {"raw", "float16", "float32", "float64"});
  std::vector<Path> z_file =
      detect_files(folder / "z.*", {"raw", "float16", "float32", "float64"});

  if (x_file.empty()) {
    x_file =
        detect_files(folder / "x0.*", {"raw", "float16", "float32", "float64"});
  }

  if (y_file.empty()) {
    y_file =
        detect_files(folder / "x1.*", {"raw", "float16", "float32", "float64"});
  }

  if (z_file.empty()) {
    z_file =
        detect_files(folder / "x2.*", {"raw", "float16", "float32", "float64"});
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

  geom_t **points = (geom_t **)calloc(ndims, sizeof(geom_t *));

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
    for (int d = 0; d < ndims; ++d) {
      free(points[d]);
    }
    free(points);
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

template <typename idx_t, typename geom_t, typename large_idx_t>
int mesh_from_folder(const MPI_Comm comm, const Path &folder,
                     // Elements
                     int *nnodesxelem_out, ptrdiff_t *n_global_elements_out,
                     ptrdiff_t *n_owned_elements_out,
                     ptrdiff_t *n_shared_elements_out,
                     ptrdiff_t *n_ghost_elements_out,
                     large_idx_t **element_mapping_out, idx_t ***elements_out,
                     // Nodes
                     int *spatial_dim_out, ptrdiff_t *n_global_nodes_out,
                     ptrdiff_t *n_owned_nodes_out, 
                     ptrdiff_t *n_shared_nodes_out,
                     ptrdiff_t *n_ghost_nodes_out,
                     large_idx_t **node_mapping_out, geom_t ***points_out,
                     // Distributed connectivities
                     int **node_owner_out, ptrdiff_t **node_offsets_out,
                     idx_t **ghosts_out) {

  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  int nnodesxelem;
  idx_t **elems;
  ptrdiff_t n_local_elements;
  ptrdiff_t n_global_elements;
  mesh_block_from_folder<idx_t>(comm, folder, &nnodesxelem, &elems,
                                &n_local_elements, &n_global_elements);

  int spatial_dim;
  geom_t **points;
  ptrdiff_t n_local2global;
  ptrdiff_t n_global_nodes;
  mesh_coordinates_from_folder(comm, folder, &spatial_dim, &points,
                               &n_local2global, &n_global_nodes);

  count_t *n2eptr;
  element_idx_t *n2e_idx;
  create_n2e<idx_t, count_t, element_idx_t>(
      comm, n_local_elements, n_global_elements, n_local2global, n_global_nodes,
      nnodesxelem, elems, &n2eptr, &n2e_idx);

  // Ensure it is always the same
  sort_n2e<count_t, element_idx_t>(n_local2global, n2eptr, n2e_idx);

  ptrdiff_t local2global_size = 0;
  // idx_t *local2global = nullptr;
  large_idx_t *local2global = nullptr;
  count_t *local_n2e_ptr = nullptr;
  element_idx_t *local_n2e_idx = nullptr;
  redistribute_n2e(comm, comm_size, comm_rank, n_local2global, n_global_nodes,
                   n_global_elements, n2eptr, n2e_idx, &local2global_size,
                   &local2global, &local_n2e_ptr, &local_n2e_idx);

  // We do not need them anymore
  free(n2eptr);
  free(n2e_idx);

  idx_t **local_elements = (idx_t **)malloc(nnodesxelem * sizeof(idx_t *));
  for (int d = 0; d < nnodesxelem; ++d) {
    local_elements[d] = (idx_t *)malloc(n_local_elements * sizeof(idx_t));
  }

  localize_element_indices(comm_size, comm_rank, n_global_elements,
                           n_local_elements, nnodesxelem, elems,
                           local2global_size, local_n2e_ptr, local_n2e_idx,
                           local2global, local_elements);
  for (int d = 0; d < nnodesxelem; ++d) {
    free(elems[d]);
  }
  free(elems);

  ptrdiff_t n_owned = 0;
  ptrdiff_t n_shared = 0;
  ptrdiff_t n_ghosts = 0;
  rearrange_local_nodes(comm_size, comm_rank, n_global_elements,
                        n_local_elements, nnodesxelem, local2global_size,
                        local_n2e_ptr, local_n2e_idx, local2global,
                        local_elements, &n_owned, &n_shared, &n_ghosts);

  large_idx_t *element_mapping =
      (large_idx_t *)malloc(n_local_elements * sizeof(large_idx_t));
  ptrdiff_t n_owned_not_shared = 0;
  rearrange_local_elements(comm_size, comm_rank, n_global_elements,
                           n_local_elements, nnodesxelem, local2global_size,
                           local_n2e_ptr, local_n2e_idx, local_elements,
                           n_owned, &n_owned_not_shared, element_mapping);

  idx_t *aura_elements = nullptr;
  idx_t **aura_element_nodes = (idx_t **)malloc(nnodesxelem * sizeof(idx_t *));
  for (int d = 0; d < nnodesxelem; ++d) {
    aura_element_nodes[d] = nullptr;
  }
  ptrdiff_t n_aura_elements = 0;
  expand_aura_elements_inconsistent(
      comm, n_global_elements, n_local_elements, nnodesxelem, local_n2e_ptr,
      local_n2e_idx, local2global, local_elements, element_mapping, n_owned,
      n_ghosts, &aura_elements, aura_element_nodes, &n_aura_elements);
  free(local_n2e_ptr);
  free(local_n2e_idx);

  long long owned_nodes_start_ll = 0;
  long long n_owned_ll = (long long)n_owned;
  SMESH_MPI_CATCH(MPI_Exscan(&n_owned_ll, &owned_nodes_start_ll, 1,
                             MPI_LONG_LONG, MPI_SUM, comm));
  if (!comm_rank) {
    owned_nodes_start_ll = 0;
  }
  const ptrdiff_t owned_nodes_start =
      static_cast<ptrdiff_t>(owned_nodes_start_ll);

  idx_t *global2owned = (idx_t *)calloc(
      rank_split(n_global_nodes, comm_size, comm_rank), sizeof(idx_t));
  prepare_node_renumbering(comm, n_global_nodes, owned_nodes_start, n_owned,
                           local2global, global2owned);

  ptrdiff_t *owned_node_ranges =
      (ptrdiff_t *)malloc((comm_size + 1) * sizeof(ptrdiff_t));
  node_ownership_ranges(comm, n_owned, owned_node_ranges);

  large_idx_t *local2global_with_aura = nullptr;

  // FIXME append aura nodes to element_mapping
  ptrdiff_t n_aura_nodes = 0;
  stitch_aura_elements(comm, n_owned, n_shared, n_ghosts, local2global,
                       nnodesxelem, n_aura_elements, aura_element_nodes,
                       n_local_elements, local_elements,
                       &local2global_with_aura, &n_aura_nodes);
  free(aura_elements);
  for (int d = 0; d < nnodesxelem; ++d) {
    free(aura_element_nodes[d]);
  }
  free(aura_element_nodes);
  free(local2global);
  local2global = local2global_with_aura;
  local2global_size = n_owned + n_ghosts + n_aura_nodes;

  SMESH_ASSERT(n_ghosts + n_aura_nodes > 0 || comm_size == 1);
  idx_t *ghost_and_aura_to_owned =
      (idx_t *)malloc((n_ghosts + n_aura_nodes) * sizeof(idx_t));
  collect_ghost_and_aura_import_indices(
      comm, n_owned, n_ghosts, n_aura_nodes, n_global_nodes, local2global,
      global2owned, owned_node_ranges, ghost_and_aura_to_owned);

  node_ownership_ranges(comm, n_owned, owned_node_ranges);
  int *owner = (int *)malloc((n_owned + n_ghosts + n_aura_nodes) * sizeof(int));
  determine_ownership(comm_size, comm_rank, n_owned, n_ghosts, n_aura_nodes,
                      ghost_and_aura_to_owned, owned_node_ranges, owner);

  group_ghost_and_aura_by_rank(comm_size, n_owned, n_ghosts, n_aura_nodes,
                               local2global, ghost_and_aura_to_owned, owner,
                               nnodesxelem, n_local_elements, n_aura_elements,
                               local_elements);

  const ptrdiff_t n_local_nodes = n_owned + n_ghosts + n_aura_nodes;
  geom_t **local_points = (geom_t **)malloc(spatial_dim * sizeof(geom_t *));
  for (int d = 0; d < spatial_dim; ++d) {
    local_points[d] = (geom_t *)malloc(n_local_nodes * sizeof(geom_t));
    gather_mapped_field(comm, n_local_nodes, n_global_nodes, local2global,
                        smesh::mpi_type<geom_t>(), points[d], local_points[d]);
    free(points[d]);
  }
  free(points);

  ptrdiff_t n_shared_elements = 0;
  for(ptrdiff_t i = 0; i < n_local_elements; ++i) {
    for(int d = 0; d < nnodesxelem; ++d) {
         idx_t node = local_elements[d][i];
         if(owner[node] != comm_rank) {
          n_shared_elements++;
            break;
      }
    }
  } 

  // Elements
  *nnodesxelem_out = nnodesxelem;
  *n_global_elements_out = n_global_elements;
  *n_owned_elements_out = n_local_elements;
  *n_shared_elements_out = n_shared_elements;
  *n_ghost_elements_out = n_aura_elements;
  *element_mapping_out = element_mapping;
  *elements_out = local_elements;

  // Nodes
  *spatial_dim_out = spatial_dim;
  *n_global_nodes_out = n_global_nodes;
  *n_owned_nodes_out = n_owned;
  *n_shared_nodes_out = n_shared;
  *n_ghost_nodes_out = n_ghosts;
  *node_mapping_out = local2global;
  *points_out = local_points;

  // Distributed connectivities
  *node_owner_out = owner;
  *node_offsets_out = owned_node_ranges;
  *ghosts_out = ghost_and_aura_to_owned;

  // Free memory that is not passed out
  free(global2owned);
  return SMESH_SUCCESS;
}

} // namespace smesh
