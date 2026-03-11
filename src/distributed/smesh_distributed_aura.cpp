#include "smesh_distributed_aura.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_INSTANTIATE_EXCHANGE_CREATE(IDX_T)                               \
  template int exchange_create<IDX_T>(                                         \
      MPI_Comm comm, const ptrdiff_t n_local_nodes,                            \
      const ptrdiff_t n_owned_nodes,                                           \
      const int *const SMESH_RESTRICT node_owner,                              \
      const ptrdiff_t *const SMESH_RESTRICT node_offsets,                      \
      const IDX_T *const SMESH_RESTRICT ghosts,                                \
      i64 *const SMESH_RESTRICT send_count,                                    \
      i64 *const SMESH_RESTRICT send_displs,                                   \
      i64 *const SMESH_RESTRICT recv_count,                                    \
      i64 *const SMESH_RESTRICT recv_displs,                                   \
      IDX_T **const SMESH_RESTRICT out_sparse_idx);                            \
  template int gather_mapped_field<IDX_T>(MPI_Comm, ptrdiff_t, ptrdiff_t,      \
                                          const IDX_T *, MPI_Datatype,         \
                                          const void *, void *)

#define SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, T)                       \
  template int exchange_scatter_add<IDX_T, T>(                                 \
      MPI_Comm, ptrdiff_t, const i64 *, const i64 *, const i64 *, const i64 *, \
      const IDX_T *, T *, T *);                                                \
  template int exchange_gather<IDX_T, T>(                                      \
      MPI_Comm, const ptrdiff_t, const i64 *const SMESH_RESTRICT,              \
      const i64 *const SMESH_RESTRICT, const i64 *const SMESH_RESTRICT,        \
      const i64 *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,      \
      T *const SMESH_RESTRICT, T *const SMESH_RESTRICT)

#define SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX(IDX_T)                  \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, f16);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, i32);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, i8);                           \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, char);                         \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, i64);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, i16);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, f32);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, f64);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, u8);                           \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, u16);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, u32);                          \
  SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, u64)

namespace smesh {

SMESH_INSTANTIATE_EXCHANGE_CREATE(i32);
SMESH_INSTANTIATE_EXCHANGE_CREATE(i64);

SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX(i32);
SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX(i64);

#if defined(__clang__)
SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(i32, long);
SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(i64, long);
#endif

} // namespace smesh

#undef SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD_FOR_IDX
#undef SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD
#undef SMESH_INSTANTIATE_EXCHANGE_CREATE
