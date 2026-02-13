#ifndef SMESH_DISTRIBUTED_AURA_HPP
#define SMESH_DISTRIBUTED_AURA_HPP

#include "smesh_distributed_base.hpp"

namespace smesh {

template <typename idx_t>
int exchange_create(MPI_Comm comm, const ptrdiff_t nnodes,
                     const ptrdiff_t n_owned_nodes, int *const node_owner,
                     const idx_t *const node_offsets, const idx_t *const ghosts,
                     int *const SMESH_RESTRICT send_count,
                     int *const SMESH_RESTRICT send_displs,
                     int *const SMESH_RESTRICT recv_count,
                     int *const SMESH_RESTRICT recv_displs,
                     idx_t **const SMESH_RESTRICT sparse_idx);

template <typename idx_t, typename T>
int exchange_scatter_add(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                          const int *const SMESH_RESTRICT send_count,
                          const int *const SMESH_RESTRICT send_displs,
                          const int *const SMESH_RESTRICT recv_count,
                          const int *const SMESH_RESTRICT recv_displs,
                          const idx_t *const SMESH_RESTRICT sparse_idx,
                          T *const SMESH_RESTRICT inout,
                          T *const SMESH_RESTRICT real_buffer);

} // namespace smesh

#endif // SMESH_DISTRIBUTED_AURA_HPP