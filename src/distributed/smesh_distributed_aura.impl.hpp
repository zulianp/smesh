#ifndef SMESH_DISTRIBUTED_AURA_IMPL_HPP
#define SMESH_DISTRIBUTED_AURA_IMPL_HPP

#include "smesh_distributed_aura.hpp"

#include "smesh_decompose.hpp"

#include <cstring>

namespace smesh {

template <typename idx_t>
int exchange_create(MPI_Comm comm, const ptrdiff_t n_local_nodes,
                    const ptrdiff_t n_owned_nodes,
                    const int *const SMESH_RESTRICT node_owner,
                    const ptrdiff_t *const SMESH_RESTRICT node_offsets,
                    const idx_t *const SMESH_RESTRICT ghosts,
                    int *const SMESH_RESTRICT send_count,
                    int *const SMESH_RESTRICT send_displs,
                    int *const SMESH_RESTRICT recv_count,
                    int *const SMESH_RESTRICT recv_displs,
                    idx_t **const SMESH_RESTRICT out_sparse_idx) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  memset(send_count, 0, size * sizeof(int));
  memset(send_displs, 0, (size + 1) * sizeof(int));
  memset(recv_count, 0, size * sizeof(int));
  memset(recv_displs, 0, (size + 1) * sizeof(int));

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
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm));

  recv_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  idx_t *remote_ghosts = (idx_t *)malloc(recv_displs[size] * sizeof(idx_t));

  // send slave nodes to process with master nodes
  SMESH_MPI_CATCH(MPI_Alltoallv(ghosts, send_count, send_displs,
                                mpi_type<idx_t>(), remote_ghosts, recv_count,
                                recv_displs, mpi_type<idx_t>(), comm));

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
  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_scatter_add(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                         const int *const SMESH_RESTRICT send_count,
                         const int *const SMESH_RESTRICT send_displs,
                         const int *const SMESH_RESTRICT recv_count,
                         const int *const SMESH_RESTRICT recv_displs,
                         const idx_t *const SMESH_RESTRICT scatter_idx,
                         T *const SMESH_RESTRICT inout,
                         T *const SMESH_RESTRICT added_buffer) {
  // Exchange ghosts
  SMESH_MPI_CATCH(MPI_Alltoallv(&inout[n_owned_nodes], send_count, send_displs,
                                mpi_type<T>(), added_buffer, recv_count,
                                recv_displs, mpi_type<T>(), comm));

  int size;
  MPI_Comm_size(comm, &size);
  ptrdiff_t count = recv_count[size - 1] + recv_displs[size - 1];
  for (ptrdiff_t i = 0; i < count; i++) {
    SMESH_ASSERT(added_buffer[i] == added_buffer[i]);
    inout[scatter_idx[i]] += added_buffer[i];
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_gather(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                    const int *const SMESH_RESTRICT send_count,
                    const int *const SMESH_RESTRICT send_displs,
                    const int *const SMESH_RESTRICT recv_count,
                    const int *const SMESH_RESTRICT recv_displs,
                    const idx_t *const SMESH_RESTRICT gather_idx,
                    T *const SMESH_RESTRICT inout,
                    T *const SMESH_RESTRICT aux_gather_buffer) {
  int size;
  MPI_Comm_size(comm, &size);
  for (ptrdiff_t i = 0; i < send_displs[size - 1] + send_count[size - 1]; i++) {
    aux_gather_buffer[i] = inout[gather_idx[i]];
  }

  SMESH_MPI_CATCH(MPI_Alltoallv(aux_gather_buffer, send_count, send_displs,
                                mpi_type<T>(), &inout[n_owned_nodes],
                                recv_count, recv_displs, mpi_type<T>(), comm));

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
  int *req_count = (int *)calloc((size_t)size, sizeof(int));
  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];
    req_count[rank_owner(n_global, idx, size)]++;
  }

  int *recv_req_count = (int *)malloc((size_t)size * sizeof(int));
  SMESH_MPI_CATCH(MPI_Alltoall(req_count, 1, smesh::mpi_type<int>(),
                               recv_req_count, 1, smesh::mpi_type<int>(),
                               comm));

  int *req_displs = (int *)malloc((size_t)size * sizeof(int));
  int *recv_displs = (int *)malloc((size_t)size * sizeof(int));

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
      (idx_t *)malloc((size_t)total_recv_req * sizeof(idx_t));

  // Exchange requested indices
  SMESH_MPI_CATCH(MPI_Alltoallv(
      req_list, req_count, req_displs, smesh::mpi_type<idx_t>(), recv_req_list,
      recv_req_count, recv_displs, smesh::mpi_type<idx_t>(), comm));

  // Build response buffer for received requests (same ordering as
  // recv_req_list)
  uint8_t *send_resp =
      (uint8_t *)malloc((size_t)total_recv_req * (size_t)type_size);

  for (ptrdiff_t i = 0; i < total_recv_req; ++i) {
    const idx_t idx = recv_req_list[i];
    const ptrdiff_t loc = (ptrdiff_t)idx - begin;
    assert(loc >= 0);
    // assert(loc < local_size);
    memcpy((void *)(send_resp + i * type_size),
           (const void *)(in + loc * type_size), (size_t)type_size);
  }

  // Exchange response data back to requesters
  uint8_t *recv_resp = (uint8_t *)malloc((size_t)n_local * (size_t)type_size);

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
  return SMESH_SUCCESS;
}
} // namespace smesh

#endif // SMESH_DISTRIBUTED_AURA_IMPL_HPP
