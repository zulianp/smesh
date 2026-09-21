#include "smesh_smooth.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SMOOTH(IDX_T, COUNT_T, GEOM_T)                                    \
    template int mesh_smooth_feature<IDX_T, COUNT_T, GEOM_T>(                                        \
            const int,                                                                               \
            const ptrdiff_t,                                                                         \
            GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                      \
            const COUNT_T *const SMESH_RESTRICT,                                                     \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const int,                                                                               \
            const GEOM_T,                                                                            \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            const GEOM_T,                                                                            \
            const GEOM_T,                                                                            \
            const GEOM_T *const SMESH_RESTRICT,                                                      \
            const GEOM_T *const SMESH_RESTRICT,                                                      \
            const GEOM_T *const SMESH_RESTRICT);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_SMOOTH(i32, i32, f32)
SMESH_EXPLICIT_INSTANTIATE_SMOOTH(i64, i32, f32)
SMESH_EXPLICIT_INSTANTIATE_SMOOTH(i64, i64, f32)
SMESH_EXPLICIT_INSTANTIATE_SMOOTH(i32, i32, f64)
SMESH_EXPLICIT_INSTANTIATE_SMOOTH(i64, i32, f64)
SMESH_EXPLICIT_INSTANTIATE_SMOOTH(i64, i64, f64)
}  // namespace smesh
