#include "smesh_distributed_reorder.impl.hpp"

#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_ORDERING(GEOM_T)                \
  template struct Hilbert3ElementOrdering<GEOM_T>

#define SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(IDX_T, GEOM_T)          \
  template int distributed_reorder_elements<IDX_T, GEOM_T,                     \
                                            Hilbert3ElementOrdering<GEOM_T>>(  \
      MPI_Comm, const int, const ptrdiff_t, const ptrdiff_t,                   \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const ptrdiff_t,      \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                      \
      Hilbert3ElementOrdering<GEOM_T>)

SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_ORDERING(f32);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_ORDERING(f64);

SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i64, f64);

#undef SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER
#undef SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_ORDERING

} // namespace smesh
