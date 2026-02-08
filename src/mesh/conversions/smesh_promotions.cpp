#include "smesh_promotions.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(IDX_T, COUNT_T, ELEM_IDX_T, GEOM_T) \
  template void mesh_tet4_to_tet15<IDX_T, COUNT_T, ELEM_IDX_T>(                 \
      const ptrdiff_t, const ptrdiff_t,                                         \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                  \
      const COUNT_T *const SMESH_RESTRICT,                                      \
      const IDX_T *const SMESH_RESTRICT,                                        \
      const ELEM_IDX_T *const SMESH_RESTRICT,                                   \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                        \
      ptrdiff_t *const SMESH_RESTRICT);                                         \
  template void mesh_tet4_to_tet15_points<IDX_T, COUNT_T, GEOM_T>(              \
      const ptrdiff_t, const ptrdiff_t,                                         \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const COUNT_T *const SMESH_RESTRICT,                                      \
      const IDX_T *const SMESH_RESTRICT,                                        \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                        \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT)

namespace smesh {
// Match the common explicit-instantiation convention used across the repo.
SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(i16, i32, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(i32, i32, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(i64, i32, i32, f32);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS