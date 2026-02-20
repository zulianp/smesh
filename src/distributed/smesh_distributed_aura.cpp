#include "smesh_distributed_aura.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_INSTANTIATE_EXCHANGE_CREATE(IDX_T)                               \
  template int exchange_create<IDX_T>(                                         \
      MPI_Comm comm, const ptrdiff_t n_local_nodes,                            \
      const ptrdiff_t n_owned_nodes,                                           \
      const int *const SMESH_RESTRICT node_owner,                              \
      const ptrdiff_t *const SMESH_RESTRICT node_offsets,                      \
      const IDX_T *const SMESH_RESTRICT ghosts,                                \
      int *const SMESH_RESTRICT send_count,                                    \
      int *const SMESH_RESTRICT send_displs,                                   \
      int *const SMESH_RESTRICT recv_count,                                    \
      int *const SMESH_RESTRICT recv_displs,                                   \
      IDX_T **const SMESH_RESTRICT out_sparse_idx);                            \
  template int gather_mapped_field<IDX_T>(MPI_Comm, ptrdiff_t, ptrdiff_t,      \
                                          const IDX_T *, MPI_Datatype,         \
                                          const void *, void *)

#define SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, T)                       \
  template int exchange_scatter_add<IDX_T, T>(                                 \
      MPI_Comm, ptrdiff_t, const int *, const int *, const int *, const int *, \
      const IDX_T *, T *, T *);                                                \
  template int exchange_gather<IDX_T, T>(                                      \
      MPI_Comm, const ptrdiff_t, const int *const SMESH_RESTRICT,              \
      const int *const SMESH_RESTRICT, const int *const SMESH_RESTRICT,        \
      const int *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,      \
      T *const SMESH_RESTRICT, T *const SMESH_RESTRICT)

#define SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX(IDX_T)                  \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, i32);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, i64);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, f32);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, f64)

namespace smesh {

SMESH_INSTANTIATE_EXCHANGE_CREATE(i32);
SMESH_INSTANTIATE_EXCHANGE_CREATE(i64);

SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX(i32);
SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX(i64);

} // namespace smesh

#undef SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX
#undef SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD
#undef SMESH_INSTANTIATE_EXCHANGE_CREATE