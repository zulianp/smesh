#ifndef SMESH_DISTRIBUTED_AURA_IMPL_HPP
#define SMESH_DISTRIBUTED_AURA_IMPL_HPP

#include "smesh_distributed_aura.hpp"

#include "smesh_alltoallv.impl.hpp"
#include "smesh_alloc.hpp"
#include "smesh_decompose.hpp"

#include <cstring>
#include <limits>

namespace smesh {

template <typename idx_t>
int exchange_create(MPI_Comm comm, const ptrdiff_t n_local_nodes,
                    const ptrdiff_t n_owned_nodes,
                    const int *const SMESH_RESTRICT node_owner,
                    const ptrdiff_t *const SMESH_RESTRICT node_offsets,
                    const idx_t *const SMESH_RESTRICT ghosts,
                    i64 *const SMESH_RESTRICT send_count,
                    i64 *const SMESH_RESTRICT send_displs,
                    i64 *const SMESH_RESTRICT recv_count,
                    i64 *const SMESH_RESTRICT recv_displs,
                    idx_t **const SMESH_RESTRICT out_sparse_idx,
                    idx_t **const SMESH_RESTRICT out_import_idx) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  memset(send_count, 0, (size_t)size * sizeof(i64));
  memset(send_displs, 0, ((size_t)size + 1) * sizeof(i64));
  memset(recv_count, 0, (size_t)size * sizeof(i64));
  memset(recv_displs, 0, ((size_t)size + 1) * sizeof(i64));

  //   Loop over all ghost nodes
  for (ptrdiff_t i = n_owned_nodes; i < n_local_nodes; i++) {
    SMESH_ASSERT(node_owner[i] >= 0);
    SMESH_ASSERT(node_owner[i] < size);

    send_count[node_owner[i]]++;
  }

  send_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    send_displs[r + 1] = send_displs[r] + send_count[r];
  }

  SMESH_MPI_CATCH(
      MPI_Alltoall(send_count, 1, mpi_type<i64>(), recv_count, 1,
                   mpi_type<i64>(), comm));

  recv_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  const ptrdiff_t n_import = n_local_nodes - n_owned_nodes;
  idx_t *ordered_import_idx = (idx_t *)SMESH_ALLOC((size_t)n_import * sizeof(idx_t));
  idx_t *packed_ghosts = (idx_t *)SMESH_ALLOC((size_t)n_import * sizeof(idx_t));
  i64 *cursor = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const int owner = node_owner[n_owned_nodes + i];
    const i64 slot = send_displs[owner] + cursor[owner]++;
    ordered_import_idx[slot] = static_cast<idx_t>(n_owned_nodes + i);
    packed_ghosts[slot] = ghosts[i];
  }
  SMESH_FREE(cursor);

  idx_t *remote_ghosts =
      (idx_t *)SMESH_ALLOC((size_t)recv_displs[size] * sizeof(idx_t));

  // send slave nodes to process with master nodes
  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64(packed_ghosts, send_count, send_displs,
                                 remote_ghosts, recv_count, recv_displs, comm,
                                 max_chunk_size));
  SMESH_FREE(packed_ghosts);

  // Replace global indexing with local indexing for identifying master nodes
  for (int r = 0; r < size; r++) {
    const int begin = recv_displs[r];
    const int extent = recv_count[r];
    idx_t *nodes = &remote_ghosts[begin];
    for (int k = 0; k < extent; k++) {
      idx_t global_idx = nodes[k];
      idx_t local_idx = global_idx - node_offsets[rank];
      nodes[k] = local_idx;
    }
  }

  *out_sparse_idx = remote_ghosts;
  *out_import_idx = ordered_import_idx;
  return SMESH_SUCCESS;
}

template <typename idx_t>
int exchange_create_ghosts(MPI_Comm comm, const ptrdiff_t n_local_nodes,
                           const ptrdiff_t n_owned_nodes,
                           const int *const SMESH_RESTRICT node_owner,
                           const ptrdiff_t *const SMESH_RESTRICT node_offsets,
                           const idx_t *const SMESH_RESTRICT ghosts,
                           i64 *const SMESH_RESTRICT send_count,
                           i64 *const SMESH_RESTRICT send_displs,
                           i64 *const SMESH_RESTRICT recv_count,
                           i64 *const SMESH_RESTRICT recv_displs,
                           idx_t **const SMESH_RESTRICT out_sparse_idx) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  memset(send_count, 0, (size_t)size * sizeof(i64));
  memset(send_displs, 0, ((size_t)size + 1) * sizeof(i64));
  memset(recv_count, 0, (size_t)size * sizeof(i64));
  memset(recv_displs, 0, ((size_t)size + 1) * sizeof(i64));

#ifndef NDEBUG
  for (ptrdiff_t i = n_owned_nodes; i < n_local_nodes - 1; i++) {
    SMESH_ASSERT(node_owner[i] <= node_owner[i + 1]);
  }
#endif

  for (ptrdiff_t i = n_owned_nodes; i < n_local_nodes; i++) {
    SMESH_ASSERT(node_owner[i] >= 0);
    SMESH_ASSERT(node_owner[i] < size);
    send_count[node_owner[i]]++;
  }

  send_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    send_displs[r + 1] = send_displs[r] + send_count[r];
  }

  SMESH_MPI_CATCH(
      MPI_Alltoall(send_count, 1, mpi_type<i64>(), recv_count, 1,
                   mpi_type<i64>(), comm));

  recv_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  idx_t *remote_ghosts =
      (idx_t *)SMESH_ALLOC((size_t)recv_displs[size] * sizeof(idx_t));

  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64(ghosts, send_count, send_displs, remote_ghosts,
                                 recv_count, recv_displs, comm, max_chunk_size));

  for (int r = 0; r < size; r++) {
    const int begin = recv_displs[r];
    const int extent = recv_count[r];
    idx_t *nodes = &remote_ghosts[begin];
    for (int k = 0; k < extent; k++) {
      const idx_t global_idx = nodes[k];
      nodes[k] = global_idx - node_offsets[rank];
    }
  }

  *out_sparse_idx = remote_ghosts;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_scatter_add(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                         const i64 *const SMESH_RESTRICT send_count,
                         const i64 *const SMESH_RESTRICT send_displs,
                         const i64 *const SMESH_RESTRICT recv_count,
                         const i64 *const SMESH_RESTRICT recv_displs,
                         const idx_t *const SMESH_RESTRICT scatter_idx,
                         const idx_t *const SMESH_RESTRICT import_idx,
                         T *const SMESH_RESTRICT inout,
                         T *const SMESH_RESTRICT send_buffer,
                         T *const SMESH_RESTRICT recv_buffer,
                         const ptrdiff_t block_size) {
  SMESH_UNUSED(n_owned_nodes);
  // Exchange ghosts
  int size;
  MPI_Comm_size(comm, &size);
  const ptrdiff_t n_import = send_displs[size];
  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const ptrdiff_t src = (ptrdiff_t)import_idx[i] * block_size;
    const ptrdiff_t dst = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      send_buffer[dst + c] = inout[src + c];
    }
  }
  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64v(send_buffer, send_count, send_displs,
                                  recv_buffer, recv_count, recv_displs,
                                  block_size, comm, max_chunk_size));

  const ptrdiff_t count = recv_displs[size];
  for (ptrdiff_t i = 0; i < count; i++) {
    const ptrdiff_t dst = (ptrdiff_t)scatter_idx[i] * block_size;
    const ptrdiff_t src = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      SMESH_ASSERT(recv_buffer[src + c] == recv_buffer[src + c]);
      inout[dst + c] += recv_buffer[src + c];
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_scatter_add_ghosts(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                                const i64 *const SMESH_RESTRICT send_count,
                                const i64 *const SMESH_RESTRICT send_displs,
                                const i64 *const SMESH_RESTRICT recv_count,
                                const i64 *const SMESH_RESTRICT recv_displs,
                                const idx_t *const SMESH_RESTRICT scatter_idx,
                                T *const SMESH_RESTRICT inout,
                                T *const SMESH_RESTRICT recv_buffer,
                                const ptrdiff_t block_size) {
  int size;
  MPI_Comm_size(comm, &size);
  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64v(&inout[n_owned_nodes * block_size],
                                  send_count, send_displs, recv_buffer,
                                  recv_count, recv_displs, block_size, comm,
                                  max_chunk_size));

  const ptrdiff_t count = recv_displs[size];
  for (ptrdiff_t i = 0; i < count; i++) {
    const ptrdiff_t dst = (ptrdiff_t)scatter_idx[i] * block_size;
    const ptrdiff_t src = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      SMESH_ASSERT(recv_buffer[src + c] == recv_buffer[src + c]);
      inout[dst + c] += recv_buffer[src + c];
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_gather(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                    const i64 *const SMESH_RESTRICT send_count,
                    const i64 *const SMESH_RESTRICT send_displs,
                    const i64 *const SMESH_RESTRICT recv_count,
                    const i64 *const SMESH_RESTRICT recv_displs,
                    const idx_t *const SMESH_RESTRICT gather_idx,
                    const idx_t *const SMESH_RESTRICT import_idx,
                    T *const SMESH_RESTRICT inout,
                    T *const SMESH_RESTRICT send_buffer,
                    T *const SMESH_RESTRICT recv_buffer,
                    const ptrdiff_t block_size) {
  SMESH_UNUSED(n_owned_nodes);
  int size;
  MPI_Comm_size(comm, &size);
  const ptrdiff_t total_send = recv_displs[size - 1] + recv_count[size - 1];
  for (ptrdiff_t i = 0; i < total_send; i++) {
    const ptrdiff_t src = (ptrdiff_t)gather_idx[i] * block_size;
    const ptrdiff_t dst = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      send_buffer[dst + c] = inout[src + c];
    }
  }

  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64v(send_buffer, recv_count, recv_displs,
                                  recv_buffer, send_count, send_displs,
                                  block_size, comm, max_chunk_size));
  const ptrdiff_t n_import = send_displs[size];
  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const ptrdiff_t dst = (ptrdiff_t)import_idx[i] * block_size;
    const ptrdiff_t src = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      inout[dst + c] = recv_buffer[src + c];
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_gather_ghosts(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                           const i64 *const SMESH_RESTRICT send_count,
                           const i64 *const SMESH_RESTRICT send_displs,
                           const i64 *const SMESH_RESTRICT recv_count,
                           const i64 *const SMESH_RESTRICT recv_displs,
                           const idx_t *const SMESH_RESTRICT gather_idx,
                           T *const SMESH_RESTRICT inout,
                           T *const SMESH_RESTRICT send_buffer,
                           const ptrdiff_t block_size) {
  int size;
  MPI_Comm_size(comm, &size);
  const ptrdiff_t total_send = recv_displs[size - 1] + recv_count[size - 1];
  for (ptrdiff_t i = 0; i < total_send; i++) {
    const ptrdiff_t src = (ptrdiff_t)gather_idx[i] * block_size;
    const ptrdiff_t dst = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      send_buffer[dst + c] = inout[src + c];
    }
  }

  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64v(send_buffer, recv_count, recv_displs,
                                  &inout[n_owned_nodes * block_size],
                                  send_count, send_displs, block_size, comm,
                                  max_chunk_size));

  return SMESH_SUCCESS;
}

// template <typename idx_t, typename T>
//  int exchange_gather(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
//                     const int *const SMESH_RESTRICT send_count,
//                     const int *const SMESH_RESTRICT send_displs,
//                     const int *const SMESH_RESTRICT recv_count,
//                     const int *const SMESH_RESTRICT recv_displs,
//                     const idx_t *const SMESH_RESTRICT gather_idx,
//                     T *const SMESH_RESTRICT inout,
//                     T *const SMESH_RESTRICT aux_gather_buffer) {
//   int size;
//   MPI_Comm_size(comm, &size);

//   ptrdiff_t total_send = 0;
// #ifndef NDEBUG
//   int expected_displ = 0;
// #endif
//   for (int r = 0; r < size; ++r) {
//     SMESH_ASSERT(send_count[r] >= 0);
// #ifndef NDEBUG
//     SMESH_ASSERT(send_displs[r] == expected_displ);
//     expected_displ += send_count[r];
// #endif
//     total_send += (ptrdiff_t)send_count[r];
//   }

//   for (ptrdiff_t i = 0; i < total_send; ++i) {
//     aux_gather_buffer[i] = inout[gather_idx[i]];
//   }

//   SMESH_MPI_CATCH(MPI_Alltoallv(aux_gather_buffer, send_count, send_displs,
//                                 mpi_type<T>(), &inout[n_owned_nodes],
//                                 recv_count, recv_displs, mpi_type<T>(),
//                                 comm));

//   return SMESH_SUCCESS;
// }

template <typename idx_t>
int gather_mapped_field(MPI_Comm comm, const ptrdiff_t n_local,
                        const ptrdiff_t n_global,
                        const idx_t *SMESH_RESTRICT const mapping,
                        MPI_Datatype data_type,
                        const void *SMESH_RESTRICT const data_in,
                        void *SMESH_RESTRICT const data_out) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  uint8_t *const in = (uint8_t *const)data_in;
  uint8_t *const out = (uint8_t *const)data_out;

  int type_size = 0;
  SMESH_MPI_CATCH(MPI_Type_size(data_type, &type_size));

  const ptrdiff_t begin = rank_start(n_global, size, rank);

  // Build request counts: how many global indices we need from each rank
  i64 *req_count = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];
    req_count[rank_owner(n_global, idx, size)]++;
  }

  i64 *recv_req_count = (i64 *)SMESH_ALLOC((size_t)size * sizeof(i64));
  SMESH_MPI_CATCH(MPI_Alltoall(req_count, 1, smesh::mpi_type<i64>(),
                               recv_req_count, 1, smesh::mpi_type<i64>(),
                               comm));

  i64 *req_displs = (i64 *)SMESH_ALLOC(((size_t)size + 1) * sizeof(i64));
  i64 *recv_displs = (i64 *)SMESH_ALLOC(((size_t)size + 1) * sizeof(i64));

  req_displs[0] = 0;
  recv_displs[0] = 0;
  for (int i = 0; i < size; ++i) {
    req_displs[i + 1] = req_displs[i] + req_count[i];
    recv_displs[i + 1] = recv_displs[i] + recv_req_count[i];
  }

  const i64 total_recv_req = recv_displs[size];

  // Pack requests: global indices + local positions (for unpack)
  idx_t *req_list = (idx_t *)SMESH_ALLOC((size_t)n_local * sizeof(idx_t));
  ptrdiff_t *local_pos =
      (ptrdiff_t *)SMESH_ALLOC((size_t)n_local * sizeof(ptrdiff_t));
  ptrdiff_t *book_keeping =
      (ptrdiff_t *)SMESH_CALLOC((size_t)size, sizeof(ptrdiff_t));

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];
    int src_rank = rank_owner(n_global, idx, size);

    const ptrdiff_t off =
        (ptrdiff_t)req_displs[src_rank] + (ptrdiff_t)book_keeping[src_rank];
    req_list[off] = idx;
    local_pos[off] = i;
    book_keeping[src_rank]++;
  }

  idx_t *recv_req_list =
      (idx_t *)SMESH_ALLOC((size_t)total_recv_req * sizeof(idx_t));

  // Exchange requested indices
  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64(req_list, req_count, req_displs, recv_req_list,
                                 recv_req_count, recv_displs, comm,
                                 max_chunk_size));

  // Build response buffer for received requests (same ordering as
  // recv_req_list)
  uint8_t *send_resp =
      (uint8_t *)SMESH_ALLOC((size_t)total_recv_req * (size_t)type_size);

  for (i64 i = 0; i < total_recv_req; ++i) {
    const idx_t idx = recv_req_list[i];
    const ptrdiff_t loc = (ptrdiff_t)idx - begin;
    assert(loc >= 0);
    // assert(loc < local_size);
    memcpy((void *)(send_resp + i * type_size),
           (const void *)(in + loc * type_size), (size_t)type_size);
  }

  // Exchange response data back to requesters
  uint8_t *recv_resp = (uint8_t *)SMESH_ALLOC((size_t)n_local * (size_t)type_size);

  SMESH_MPI_CATCH(all_to_allv_64_b(send_resp, recv_req_count, recv_displs,
                                data_type, recv_resp, req_count, req_displs,
                                data_type, comm, max_chunk_size));

  // Unpack into local ordering
  for (ptrdiff_t off = 0; off < n_local; ++off) {
    const ptrdiff_t i = local_pos[off];
    memcpy((void *)(out + i * type_size),
           (const void *)(recv_resp + off * type_size), (size_t)type_size);
  }

  SMESH_FREE(recv_resp);
  SMESH_FREE(send_resp);
  SMESH_FREE(recv_req_list);
  SMESH_FREE(book_keeping);
  SMESH_FREE(local_pos);
  SMESH_FREE(req_list);
  SMESH_FREE(recv_displs);
  SMESH_FREE(req_displs);
  SMESH_FREE(recv_req_count);
  SMESH_FREE(req_count);
  return SMESH_SUCCESS;
}
} // namespace smesh

#endif // SMESH_DISTRIBUTED_AURA_IMPL_HPP
