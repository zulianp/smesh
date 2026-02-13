#ifndef SMESH_DISTRIBUTED_AURA_IMPL_HPP
#define SMESH_DISTRIBUTED_AURA_IMPL_HPP

#include "smesh_distributed_aura.hpp"

#include <cstring>

namespace smesh {

template <typename idx_t>
int exchange_create(MPI_Comm comm, const ptrdiff_t nnodes,
                    const ptrdiff_t n_owned_nodes, int *const node_owner,
                    const idx_t *const node_offsets, const idx_t *const ghosts,
                    int *const SMESH_RESTRICT send_count,
                    int *const SMESH_RESTRICT send_displs,
                    int *const SMESH_RESTRICT recv_count,
                    int *const SMESH_RESTRICT recv_displs,
                    idx_t **const SMESH_RESTRICT sparse_idx) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  memset(send_count, 0, size * sizeof(int));
  memset(send_displs, 0, (size + 1) * sizeof(int));
  memset(recv_count, 0, size * sizeof(int));
  memset(recv_displs, 0, (size + 1) * sizeof(int));

  //   Loop over all ghost nodes
  for (ptrdiff_t i = n_owned_nodes; i < nnodes; i++) {
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

  idx_t *slave_nodes = (idx_t *)malloc(recv_displs[size] * sizeof(idx_t));

  // send slave nodes to process with master nodes
  SMESH_MPI_CATCH(MPI_Alltoallv(ghosts, send_count, send_displs,
                                mpi_type<idx_t>(), slave_nodes, recv_count,
                                recv_displs, mpi_type<idx_t>(), comm));

  // Replace global indexing with local indexing for identifying master nodes
  for (int r = 0; r < size; r++) {
    const int begin = recv_displs[r];
    const int extent = recv_count[r];
    idx_t *nodes = &slave_nodes[begin];
    for (int k = 0; k < extent; k++) {
      idx_t global_idx = nodes[k];
      idx_t local_idx = global_idx - node_offsets[rank];
      nodes[k] = local_idx;
    }
  }

  *sparse_idx = slave_nodes;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int exchange_scatter_add(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                         const int *const SMESH_RESTRICT send_count,
                         const int *const SMESH_RESTRICT send_displs,
                         const int *const SMESH_RESTRICT recv_count,
                         const int *const SMESH_RESTRICT recv_displs,
                         const idx_t *const SMESH_RESTRICT sparse_idx,
                         T *const SMESH_RESTRICT inout,
                         T *const SMESH_RESTRICT real_buffer) {
  // Exchange ghosts
  SMESH_MPI_CATCH(MPI_Alltoallv(&inout[n_owned_nodes], send_count, send_displs,
                                mpi_type<T>(), real_buffer, recv_count,
                                recv_displs, mpi_type<T>(), comm));

  int size;
  MPI_Comm_size(comm, &size);
  ptrdiff_t count = recv_count[size - 1] + recv_displs[size - 1];
  for (ptrdiff_t i = 0; i < count; i++) {
    SMESH_ASSERT(real_buffer[i] == real_buffer[i]);
    inout[sparse_idx[i]] += real_buffer[i];
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_DISTRIBUTED_AURA_IMPL_HPP