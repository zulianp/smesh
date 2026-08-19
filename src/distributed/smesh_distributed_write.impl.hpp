#define SMESH_DISTRIBUTED_WRITE_IMPL_HPP

#include "matrixio_array.h"
#include "smesh_alloc.hpp"
#include "smesh_alltoallv.impl.hpp"
#include "smesh_base.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <limits>

namespace smesh {

// Map sparse unique GIDs onto [0, n_global) once. Same rank_owner Alltoallv
// as unique_by_id / parallel read: occupancy lives in a GID-space slice
// (n_space / P), not a replicated n_global buffer. No comparison sort.
// Identity when mapping already lies in [0, n_global).
template <typename large_idx_t>
static int densify_file_indices_rank_owner(MPI_Comm           comm,
                                           const ptrdiff_t    n_local,
                                           const ptrdiff_t    n_global,
                                           const large_idx_t *mapping,
                                           large_idx_t       *file_idx) {
    int rank = 0;
    int size = 1;
    SMESH_MPI_CATCH(MPI_Comm_rank(comm, &rank));
    SMESH_MPI_CATCH(MPI_Comm_size(comm, &size));
    if (n_global < 0 || n_local < 0) {
        return SMESH_FAILURE;
    }
    if (n_global == 0) {
        return n_local == 0 ? SMESH_SUCCESS : SMESH_FAILURE;
    }
    if (n_local > 0 && (!mapping || !file_idx)) {
        return SMESH_FAILURE;
    }

    large_idx_t local_max = (large_idx_t)-1;
    for (ptrdiff_t i = 0; i < n_local; ++i) {
        if (mapping[i] > local_max) {
            local_max = mapping[i];
        }
    }
    large_idx_t gmax = (large_idx_t)-1;
    SMESH_MPI_CATCH(MPI_Allreduce(&local_max, &gmax, 1, mpi_type<large_idx_t>(), MPI_MAX, comm));
    if (gmax < 0) {
        return SMESH_FAILURE;
    }
    if (gmax < (large_idx_t)n_global) {
        if (n_local > 0) {
            memcpy(file_idx, mapping, (size_t)n_local * sizeof(large_idx_t));
        }
        return SMESH_SUCCESS;
    }

    const ptrdiff_t n_space     = (ptrdiff_t)gmax + 1;
    const ptrdiff_t n_space_pad = n_space >= (ptrdiff_t)size ? n_space : (ptrdiff_t)size;
    const ptrdiff_t ids_start   = rank_start(n_space_pad, size, rank);
    const ptrdiff_t n_owned_ids = rank_split(n_space_pad, size, rank);

    constexpr i64 k_chunk = static_cast<i64>(1) << 20;
    constexpr int k_nreq  = 2;
    constexpr int k_nrep  = 2;

    i64 *send_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
    i64 *send_count  = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
    i64 *recv_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
    i64 *recv_count  = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
    if (!send_displs || !send_count || !recv_displs || !recv_count) {
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        return SMESH_FAILURE;
    }

    for (ptrdiff_t i = 0; i < n_local; ++i) {
        const int p = rank_owner(n_space_pad, (ptrdiff_t)mapping[i], size);
        send_displs[p + 1]++;
    }
    SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count, 1, mpi_type<i64>(), comm));
    recv_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        recv_displs[r + 1] = recv_displs[r] + recv_count[r];
    }
    send_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] += send_displs[r];
    }

    const ptrdiff_t n_send = (ptrdiff_t)send_displs[size];
    const ptrdiff_t n_recv = (ptrdiff_t)recv_displs[size];
    large_idx_t *send_pack = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_send, 1) * k_nreq * sizeof(large_idx_t));
    large_idx_t *recv_pack = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_recv, 1) * k_nreq * sizeof(large_idx_t));
    if (!send_pack || !recv_pack) {
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        SMESH_FREE(send_pack);
        SMESH_FREE(recv_pack);
        return SMESH_FAILURE;
    }

    for (ptrdiff_t i = 0; i < n_local; ++i) {
        const int p    = rank_owner(n_space_pad, (ptrdiff_t)mapping[i], size);
        const i64 slot = send_displs[p] + send_count[p]++;
        large_idx_t *row = send_pack + (ptrdiff_t)slot * k_nreq;
        row[0] = mapping[i];
        row[1] = (large_idx_t)i;
    }

    if (all_to_allv_64v(send_pack,
                        send_count,
                        send_displs,
                        recv_pack,
                        recv_count,
                        recv_displs,
                        k_nreq,
                        comm,
                        k_chunk) != SMESH_SUCCESS) {
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        SMESH_FREE(send_pack);
        SMESH_FREE(recv_pack);
        return SMESH_FAILURE;
    }

    uint8_t     *occ   = (uint8_t *)SMESH_CALLOC((size_t)n_owned_ids, sizeof(uint8_t));
    large_idx_t *dense = (large_idx_t *)SMESH_ALLOC((size_t)n_owned_ids * sizeof(large_idx_t));
    if (!occ || !dense) {
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        SMESH_FREE(send_pack);
        SMESH_FREE(recv_pack);
        SMESH_FREE(occ);
        SMESH_FREE(dense);
        return SMESH_FAILURE;
    }

    int dup = 0;
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const ptrdiff_t off = (ptrdiff_t)recv_pack[i * k_nreq] - ids_start;
        if (off < 0 || off >= n_owned_ids) {
            dup = 1;
            break;
        }
        if (occ[off]) {
            dup = 1;
            break;
        }
        occ[off] = 1;
    }
    int dup_any = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&dup, &dup_any, 1, MPI_INT, MPI_MAX, comm));
    if (dup_any) {
        SMESH_ERROR("write_mapped_field: mapping is not unique\n");
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        SMESH_FREE(send_pack);
        SMESH_FREE(recv_pack);
        SMESH_FREE(occ);
        SMESH_FREE(dense);
        return SMESH_FAILURE;
    }

    ptrdiff_t n_occ = 0;
    for (ptrdiff_t off = 0; off < n_owned_ids; ++off) {
        n_occ += (ptrdiff_t)occ[off];
    }
    ptrdiff_t base = 0;
    SMESH_MPI_CATCH(MPI_Exscan(&n_occ, &base, 1, mpi_type<ptrdiff_t>(), MPI_SUM, comm));
    if (rank == 0) {
        base = 0;
    }
    ptrdiff_t n_occ_sum = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&n_occ, &n_occ_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, comm));
    if (n_occ_sum != n_global) {
        SMESH_ERROR("write_mapped_field: unique gids %ld != n_global=%ld\n",
                    (long)n_occ_sum,
                    (long)n_global);
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        SMESH_FREE(send_pack);
        SMESH_FREE(recv_pack);
        SMESH_FREE(occ);
        SMESH_FREE(dense);
        return SMESH_FAILURE;
    }

    ptrdiff_t k = 0;
    for (ptrdiff_t off = 0; off < n_owned_ids; ++off) {
        if (occ[off]) {
            dense[off] = (large_idx_t)(base + k);
            ++k;
        }
    }

    memset(send_displs, 0, ((size_t)size + 1) * sizeof(i64));
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] = recv_count[r];
    }
    send_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] += send_displs[r];
    }
    memcpy(send_count, recv_count, (size_t)size * sizeof(i64));

    i64 *reply_recv_count  = (i64 *)SMESH_ALLOC((size_t)size * sizeof(i64));
    i64 *reply_recv_displs = (i64 *)SMESH_ALLOC(((size_t)size + 1) * sizeof(i64));
    SMESH_MPI_CATCH(MPI_Alltoall(send_count, 1, mpi_type<i64>(), reply_recv_count, 1, mpi_type<i64>(), comm));
    reply_recv_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        reply_recv_displs[r + 1] = reply_recv_displs[r] + reply_recv_count[r];
    }

    const ptrdiff_t n_reply_send = (ptrdiff_t)send_displs[size];
    const ptrdiff_t n_reply_recv = (ptrdiff_t)reply_recv_displs[size];
    large_idx_t *reply_send = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_reply_send, 1) * k_nrep * sizeof(large_idx_t));
    large_idx_t *reply_recv = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_reply_recv, 1) * k_nrep * sizeof(large_idx_t));

    memset(send_count, 0, (size_t)size * sizeof(i64));
    for (int r = 0; r < size; ++r) {
        for (i64 j = 0; j < recv_count[r]; ++j) {
            const ptrdiff_t i    = (ptrdiff_t)recv_displs[r] + (ptrdiff_t)j;
            const ptrdiff_t off  = (ptrdiff_t)recv_pack[i * k_nreq] - ids_start;
            const i64       slot = send_displs[r] + send_count[r]++;
            large_idx_t    *row  = reply_send + (ptrdiff_t)slot * k_nrep;
            row[0] = recv_pack[i * k_nreq + 1];
            row[1] = dense[off];
        }
    }

    const int reply_err = all_to_allv_64v(reply_send,
                                          send_count,
                                          send_displs,
                                          reply_recv,
                                          reply_recv_count,
                                          reply_recv_displs,
                                          k_nrep,
                                          comm,
                                          k_chunk);

    SMESH_FREE(send_displs);
    SMESH_FREE(send_count);
    SMESH_FREE(recv_displs);
    SMESH_FREE(recv_count);
    SMESH_FREE(send_pack);
    SMESH_FREE(recv_pack);
    SMESH_FREE(occ);
    SMESH_FREE(dense);
    SMESH_FREE(reply_send);

    if (reply_err != SMESH_SUCCESS) {
        SMESH_FREE(reply_recv_count);
        SMESH_FREE(reply_recv_displs);
        SMESH_FREE(reply_recv);
        return SMESH_FAILURE;
    }

    for (ptrdiff_t i = 0; i < n_reply_recv; ++i) {
        const ptrdiff_t local_i = (ptrdiff_t)reply_recv[i * k_nrep];
        if (local_i < 0 || local_i >= n_local) {
            SMESH_ERROR("write_mapped_field: densify reply index out of range\n");
            SMESH_FREE(reply_recv_count);
            SMESH_FREE(reply_recv_displs);
            SMESH_FREE(reply_recv);
            return SMESH_FAILURE;
        }
        file_idx[local_i] = reply_recv[i * k_nrep + 1];
    }

    SMESH_FREE(reply_recv_count);
    SMESH_FREE(reply_recv_displs);
    SMESH_FREE(reply_recv);
    return SMESH_SUCCESS;
}

template <typename FileType, typename T>
int array_write_convert(MPI_Comm comm, const Path &path,
                        const T *const SMESH_RESTRICT data,
                        const ptrdiff_t n_local_elements,
                        const ptrdiff_t n_global_elements) {
  if (std::is_same_v<FileType, T>) {
    return array_write(comm, path.c_str(), smesh::mpi_type<T>(), data,
                       n_local_elements, n_global_elements);
  }

  FileType *buffer = (FileType *)SMESH_ALLOC(n_local_elements * sizeof(FileType));
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    buffer[i] = static_cast<FileType>(data[i]);
  }
  int ret = array_write(comm, path.c_str(), smesh::mpi_type<FileType>(), buffer,
                        n_local_elements, n_global_elements);

  SMESH_FREE(buffer);
  return ret;
}

template <typename T>
int array_write_convert_from_extension(MPI_Comm comm, const Path &path,
                                       const T *const SMESH_RESTRICT data,
                                       const ptrdiff_t n_local_elements,
                                       const ptrdiff_t n_global_elements) {
  auto ext = path.extension();
  if (ext == "raw") {
    return array_write(comm, path.c_str(), smesh::mpi_type<T>(), data,
                       n_local_elements, n_global_elements);
  } else if (ext == "float16") {
    return array_write_convert<f16, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "float32") {
    return array_write_convert<f32, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "float64") {
    return array_write_convert<f64, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "int16") {
    return array_write_convert<i16, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "int32") {
    return array_write_convert<i32, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "int64") {
    return array_write_convert<i64, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  }
  // else  if (ext == "txt") {
  //   return SMESH_SUCCESS;
  // }
  else {
    SMESH_ERROR("Unsupported file extension %s for file %s\n", ext.c_str(),
                path.c_str());
    return SMESH_FAILURE;
  }
}

template <typename large_idx_t>
int write_mapped_field(MPI_Comm comm, const Path &output_path,
                       const ptrdiff_t n_local, const ptrdiff_t n_global,
                       const large_idx_t *const mapping, MPI_Datatype data_type,
                       const void *const data_in) {
  using byte_t = uint8_t;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  const byte_t *const data = (const byte_t *const)data_in;

  int type_size;
  SMESH_MPI_CATCH(MPI_Type_size(data_type, &type_size));

  const ptrdiff_t local_output_size_no_remainder = n_global / size;
  const ptrdiff_t begin = (n_global / size) * rank;

  ptrdiff_t local_output_size = local_output_size_no_remainder;
  if (rank == size - 1) {
    local_output_size = n_global - begin;
  }

  i64 *send_count = (i64 *)SMESH_ALLOC((size) * sizeof(i64));
  memset(send_count, 0, (size) * sizeof(i64));

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const large_idx_t idx = mapping[i];
    int dest_rank = size - 1;
    if (local_output_size_no_remainder > 0) {
      dest_rank = std::min(size - 1, (int)(idx / (large_idx_t)local_output_size_no_remainder));
    }
    send_count[dest_rank]++;
  }

  i64 *recv_count = (i64 *)SMESH_ALLOC((size) * sizeof(i64));
  SMESH_MPI_CATCH(MPI_Alltoall(send_count, 1, mpi_type<i64>(), recv_count, 1,
                               mpi_type<i64>(), comm));

  i64 *send_displs = (i64 *)SMESH_ALLOC(size * sizeof(i64));
  i64 *recv_displs = (i64 *)SMESH_ALLOC(size * sizeof(i64));
  i64 *book_keeping = (i64 *)SMESH_CALLOC(size, sizeof(i64));

  send_displs[0] = 0;
  recv_displs[0] = 0;

  // Create data displacements for sending
  for (int i = 0; i < size - 1; ++i) {
    send_displs[i + 1] = send_displs[i] + send_count[i];
  }

  // Create data displacements for receiving
  for (int i = 0; i < size - 1; ++i) {
    recv_displs[i + 1] = recv_displs[i] + recv_count[i];
  }

  const i64 n_send = send_displs[size - 1] + send_count[size - 1];
  const i64 n_recv = recv_displs[size - 1] + recv_count[size - 1];

  large_idx_t *send_list = (large_idx_t *)SMESH_ALLOC((size_t)std::max<i64>(n_send, 1) * sizeof(large_idx_t));

  ptrdiff_t n_buff = std::max((ptrdiff_t)n_send, local_output_size);
  n_buff = std::max<ptrdiff_t>(n_buff, 1);
  uint8_t *send_data_and_final_storage = (uint8_t *)SMESH_ALLOC((size_t)n_buff * (size_t)type_size);

  // Pack data and indices
  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const large_idx_t idx = mapping[i];
    int dest_rank = size - 1;
    if (local_output_size_no_remainder > 0) {
      dest_rank = std::min(size - 1, (int)(idx / (large_idx_t)local_output_size_no_remainder));
    }
    SMESH_ASSERT(dest_rank < size);

    // Put index and data into buffers
    const ptrdiff_t offset = send_displs[dest_rank] + book_keeping[dest_rank];
    send_list[offset] = idx;
    memcpy((void *)&send_data_and_final_storage[offset * type_size],
           (void *)&data[i * type_size], type_size);

    book_keeping[dest_rank]++;
  }

  large_idx_t *recv_list = (large_idx_t *)SMESH_ALLOC((size_t)std::max<i64>(n_recv, 1) * sizeof(large_idx_t));
  uint8_t *recv_data = (uint8_t *)SMESH_ALLOC((size_t)std::max<i64>(n_recv, 1) * (size_t)type_size);

  ///////////////////////////////////
  // Send indices
  ///////////////////////////////////

  ptrdiff_t max_chunk_size = std::numeric_limits<ptrdiff_t>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64_b(send_list, send_count, send_displs,
                                mpi_type<large_idx_t>(), recv_list, recv_count,
                                recv_displs, mpi_type<large_idx_t>(), comm, max_chunk_size));

  ///////////////////////////////////
  // Send data
  ///////////////////////////////////

  SMESH_MPI_CATCH(all_to_allv_64_b(send_data_and_final_storage, send_count,
                                send_displs, data_type, recv_data, recv_count,
                                recv_displs, data_type, comm, max_chunk_size));

  ///////////////////////////////////
  // Unpack indexed data
  ///////////////////////////////////

  if (local_output_size > 0) {
    memset(send_data_and_final_storage, 0, (size_t)local_output_size * (size_t)type_size);
  }
  for (i64 i = 0; i < n_recv; ++i) {
    const ptrdiff_t dest = (ptrdiff_t)recv_list[i] - begin;
    if (dest < 0 || dest >= local_output_size) {
      SMESH_ERROR("write_mapped_field: dest %ld out of range [0, %ld)\n",
                  (long)dest,
                  (long)local_output_size);
      SMESH_FREE(send_count);
      SMESH_FREE(send_displs);
      SMESH_FREE(recv_count);
      SMESH_FREE(recv_displs);
      SMESH_FREE(book_keeping);
      SMESH_FREE(send_list);
      SMESH_FREE(recv_list);
      SMESH_FREE(recv_data);
      SMESH_FREE(send_data_and_final_storage);
      return SMESH_FAILURE;
    }
    memcpy((void *)&send_data_and_final_storage[dest * type_size],
           (void *)&recv_data[i * type_size], type_size);
  }

  array_write(comm, output_path.c_str(), data_type, (void *)send_data_and_final_storage,
              local_output_size, n_global);

  ///////////////////////////////////
  // Clean-up
  ///////////////////////////////////
  SMESH_FREE(send_count);
  SMESH_FREE(send_displs);
  SMESH_FREE(recv_count);
  SMESH_FREE(recv_displs);
  SMESH_FREE(book_keeping);
  SMESH_FREE(send_list);
  SMESH_FREE(recv_list);
  SMESH_FREE(recv_data);
  SMESH_FREE(send_data_and_final_storage);
  return 0;
}

int write_distributed_mesh_coordinates(
    MPI_Comm comm, const Path &path, int spatial_dim, ptrdiff_t n_global_nodes,
    ptrdiff_t n_owned_nodes, const large_idx_t *node_mapping,
    geom_t **local_points) {
  SMESH_TRACE_SCOPE("write_distributed_mesh_coordinates");

  int err = SMESH_SUCCESS;
  static constexpr char xyz[3] = {'x', 'y', 'z'};
  for (int d = 0; d < spatial_dim; ++d) {
    std::string fname =
        std::string(1, xyz[d]) + "." + std::string(TypeToString<geom_t>::value());
    Path coord_path = path / fname;
    err |= write_mapped_field(comm, coord_path, n_owned_nodes, n_global_nodes,
                              node_mapping, smesh::mpi_type<geom_t>(),
                              local_points[d]);
  }
  return err == SMESH_SUCCESS ? SMESH_SUCCESS : SMESH_FAILURE;
}

int write_distributed_block_connectivity(
    MPI_Comm comm, const Path &path, ptrdiff_t n_global_elements,
    ptrdiff_t n_owned_elements, const large_idx_t *element_mapping,
    int nnodesxelem, idx_t **local_elements, const large_idx_t *node_mapping) {
  SMESH_TRACE_SCOPE("write_distributed_block_connectivity");

  large_idx_t *file_idx = nullptr;
  if (n_owned_elements > 0) {
    file_idx = (large_idx_t *)SMESH_ALLOC((size_t)n_owned_elements * sizeof(large_idx_t));
    if (!file_idx) {
      return SMESH_FAILURE;
    }
  }
  if (densify_file_indices_rank_owner(comm,
                                      n_owned_elements,
                                      n_global_elements,
                                      element_mapping,
                                      file_idx) != SMESH_SUCCESS) {
    SMESH_FREE(file_idx);
    return SMESH_FAILURE;
  }

  int err = SMESH_SUCCESS;
  for (int v = 0; v < nnodesxelem; ++v) {
    std::string fname =
        "i" + std::to_string(v) + "." + std::string(TypeToString<idx_t>::value());
    Path conn_path = path / fname;

    idx_t *buffer =
        (idx_t *)SMESH_ALLOC((size_t)n_owned_elements * sizeof(idx_t));
    if (!buffer) {
      SMESH_FREE(file_idx);
      return SMESH_FAILURE;
    }

    idx_t *local_col = local_elements[v];
    for (ptrdiff_t e = 0; e < n_owned_elements; ++e) {
      const idx_t local_node = local_col[e];
      buffer[e] = static_cast<idx_t>(node_mapping[local_node]);
    }

    err |= write_mapped_field(comm, conn_path, n_owned_elements,
                              n_global_elements, file_idx,
                              smesh::mpi_type<idx_t>(), buffer);

    SMESH_FREE(buffer);
  }

  SMESH_FREE(file_idx);
  return err == SMESH_SUCCESS ? SMESH_SUCCESS : SMESH_FAILURE;
}

int write_distributed_mesh_topology(
    MPI_Comm comm, const Path &path, enum ElemType /*element_type*/,
    int spatial_dim, ptrdiff_t n_global_elements, ptrdiff_t n_owned_elements,
    const large_idx_t *element_mapping, int nnodesxelem, idx_t **local_elements,
    ptrdiff_t n_global_nodes, ptrdiff_t n_owned_nodes,
    const large_idx_t *node_mapping, geom_t **local_points) {
  SMESH_TRACE_SCOPE("write_distributed_mesh_topology");

  int err = write_distributed_mesh_coordinates(
      comm, path, spatial_dim, n_global_nodes, n_owned_nodes, node_mapping,
      local_points);
  if (err != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  err = write_distributed_block_connectivity(
      comm, path, n_global_elements, n_owned_elements, element_mapping,
      nnodesxelem, local_elements, node_mapping);

  return err == SMESH_SUCCESS ? SMESH_SUCCESS : SMESH_FAILURE;
}

} // namespace smesh
