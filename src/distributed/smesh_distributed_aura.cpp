#include "smesh_distributed_aura.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_INSTANTIATE_EXCHANGE_CREATE(IDX_T)                               \
  template int exchange_create<IDX_T>(MPI_Comm, ptrdiff_t, ptrdiff_t, int *,  \
                                       const IDX_T *, const IDX_T *, int *,    \
                                       int *, int *, int *, IDX_T **)

#define SMESH_INSTANTIATE_EXCHANGE_SCATTER_ADD(IDX_T, T)                       \
  template int exchange_scatter_add<IDX_T, T>(                                \
      MPI_Comm, ptrdiff_t, const int *, const int *, const int *, const int *, \
      const IDX_T *, T *, T *)

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