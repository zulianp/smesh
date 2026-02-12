#include "smesh_refine.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_REFINEMENT(IDX_T, COUNT_T,   \
                                              GEOM_T)                          \
  template int mesh_refine<IDX_T, COUNT_T, GEOM_T>(             \
      const enum ElemType, const ptrdiff_t,                                    \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const int,      \
      const ptrdiff_t,                                                         \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                \
      const COUNT_T *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,  \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_REFINEMENT(i32, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_REFINEMENT(i64, i64, f32);
} // namespace smesh
