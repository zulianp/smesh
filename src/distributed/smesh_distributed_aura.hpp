#ifndef SMESH_DISTRIBUTED_AURA_HPP
#define SMESH_DISTRIBUTED_AURA_HPP

#include "smesh_distributed_base.hpp"

#include <cstddef>
#include <mpi.h>

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
                    idx_t **const SMESH_RESTRICT out_import_idx);

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
                           idx_t **const SMESH_RESTRICT out_sparse_idx);

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
                         const ptrdiff_t block_size);

template <typename idx_t, typename T>
inline int exchange_scatter_add(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                                const i64 *const SMESH_RESTRICT send_count,
                                const i64 *const SMESH_RESTRICT send_displs,
                                const i64 *const SMESH_RESTRICT recv_count,
                                const i64 *const SMESH_RESTRICT recv_displs,
                                const idx_t *const SMESH_RESTRICT scatter_idx,
                                const idx_t *const SMESH_RESTRICT import_idx,
                                T *const SMESH_RESTRICT inout,
                                T *const SMESH_RESTRICT send_buffer,
                                T *const SMESH_RESTRICT recv_buffer);

template <typename idx_t>
int gather_mapped_field(MPI_Comm comm, const ptrdiff_t n_local,
                        const ptrdiff_t n_global,
                        const idx_t *SMESH_RESTRICT const mapping,
                        MPI_Datatype data_type,
                        const void *SMESH_RESTRICT const data_in,
                        void *SMESH_RESTRICT const data_out);

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
                    const ptrdiff_t block_size);

template <typename idx_t, typename T>
inline int exchange_gather(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                           const i64 *const SMESH_RESTRICT send_count,
                           const i64 *const SMESH_RESTRICT send_displs,
                           const i64 *const SMESH_RESTRICT recv_count,
                           const i64 *const SMESH_RESTRICT recv_displs,
                           const idx_t *const SMESH_RESTRICT gather_idx,
                           const idx_t *const SMESH_RESTRICT import_idx,
                           T *const SMESH_RESTRICT inout,
                           T *const SMESH_RESTRICT send_buffer,
                           T *const SMESH_RESTRICT recv_buffer);

template <typename idx_t, typename T>
int exchange_scatter_add_ghosts(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                                const i64 *const SMESH_RESTRICT send_count,
                                const i64 *const SMESH_RESTRICT send_displs,
                                const i64 *const SMESH_RESTRICT recv_count,
                                const i64 *const SMESH_RESTRICT recv_displs,
                                const idx_t *const SMESH_RESTRICT scatter_idx,
                                T *const SMESH_RESTRICT inout,
                                T *const SMESH_RESTRICT recv_buffer,
                                const ptrdiff_t block_size);

template <typename idx_t, typename T>
inline int exchange_scatter_add_ghosts(MPI_Comm comm,
                                       const ptrdiff_t n_owned_nodes,
                                       const i64 *const SMESH_RESTRICT send_count,
                                       const i64 *const SMESH_RESTRICT send_displs,
                                       const i64 *const SMESH_RESTRICT recv_count,
                                       const i64 *const SMESH_RESTRICT recv_displs,
                                       const idx_t *const SMESH_RESTRICT scatter_idx,
                                       T *const SMESH_RESTRICT inout,
                                       T *const SMESH_RESTRICT recv_buffer);

template <typename idx_t, typename T>
int exchange_gather_ghosts(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                           const i64 *const SMESH_RESTRICT send_count,
                           const i64 *const SMESH_RESTRICT send_displs,
                           const i64 *const SMESH_RESTRICT recv_count,
                           const i64 *const SMESH_RESTRICT recv_displs,
                           const idx_t *const SMESH_RESTRICT gather_idx,
                           T *const SMESH_RESTRICT inout,
                           T *const SMESH_RESTRICT send_buffer,
                           const ptrdiff_t block_size);

template <typename idx_t, typename T>
inline int exchange_gather_ghosts(MPI_Comm comm,
                                  const ptrdiff_t n_owned_nodes,
                                  const i64 *const SMESH_RESTRICT send_count,
                                  const i64 *const SMESH_RESTRICT send_displs,
                                  const i64 *const SMESH_RESTRICT recv_count,
                                  const i64 *const SMESH_RESTRICT recv_displs,
                                  const idx_t *const SMESH_RESTRICT gather_idx,
                                  T *const SMESH_RESTRICT inout,
                                  T *const SMESH_RESTRICT send_buffer);

template <typename idx_t, typename T>
inline int exchange_scatter_add(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                                const i64 *const SMESH_RESTRICT send_count,
                                const i64 *const SMESH_RESTRICT send_displs,
                                const i64 *const SMESH_RESTRICT recv_count,
                                const i64 *const SMESH_RESTRICT recv_displs,
                                const idx_t *const SMESH_RESTRICT scatter_idx,
                                const idx_t *const SMESH_RESTRICT import_idx,
                                T *const SMESH_RESTRICT inout,
                                T *const SMESH_RESTRICT send_buffer,
                                T *const SMESH_RESTRICT recv_buffer) {
  return exchange_scatter_add(comm, n_owned_nodes, send_count, send_displs,
                              recv_count, recv_displs, scatter_idx, import_idx,
                              inout, send_buffer, recv_buffer, 1);
}

template <typename idx_t, typename T>
inline int exchange_gather(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                           const i64 *const SMESH_RESTRICT send_count,
                           const i64 *const SMESH_RESTRICT send_displs,
                           const i64 *const SMESH_RESTRICT recv_count,
                           const i64 *const SMESH_RESTRICT recv_displs,
                           const idx_t *const SMESH_RESTRICT gather_idx,
                           const idx_t *const SMESH_RESTRICT import_idx,
                           T *const SMESH_RESTRICT inout,
                           T *const SMESH_RESTRICT send_buffer,
                           T *const SMESH_RESTRICT recv_buffer) {
  return exchange_gather(comm, n_owned_nodes, send_count, send_displs,
                         recv_count, recv_displs, gather_idx, import_idx, inout,
                         send_buffer, recv_buffer, 1);
}

template <typename idx_t, typename T>
inline int exchange_scatter_add_ghosts(MPI_Comm comm,
                                       const ptrdiff_t n_owned_nodes,
                                       const i64 *const SMESH_RESTRICT send_count,
                                       const i64 *const SMESH_RESTRICT send_displs,
                                       const i64 *const SMESH_RESTRICT recv_count,
                                       const i64 *const SMESH_RESTRICT recv_displs,
                                       const idx_t *const SMESH_RESTRICT scatter_idx,
                                       T *const SMESH_RESTRICT inout,
                                       T *const SMESH_RESTRICT recv_buffer) {
  return exchange_scatter_add_ghosts(comm, n_owned_nodes, send_count,
                                     send_displs, recv_count, recv_displs,
                                     scatter_idx, inout, recv_buffer, 1);
}

template <typename idx_t, typename T>
inline int exchange_gather_ghosts(MPI_Comm comm,
                                  const ptrdiff_t n_owned_nodes,
                                  const i64 *const SMESH_RESTRICT send_count,
                                  const i64 *const SMESH_RESTRICT send_displs,
                                  const i64 *const SMESH_RESTRICT recv_count,
                                  const i64 *const SMESH_RESTRICT recv_displs,
                                  const idx_t *const SMESH_RESTRICT gather_idx,
                                  T *const SMESH_RESTRICT inout,
                                  T *const SMESH_RESTRICT send_buffer) {
  return exchange_gather_ghosts(comm, n_owned_nodes, send_count, send_displs,
                                recv_count, recv_displs, gather_idx, inout,
                                send_buffer, 1);
}
} // namespace smesh

#endif // SMESH_DISTRIBUTED_AURA_HPP
