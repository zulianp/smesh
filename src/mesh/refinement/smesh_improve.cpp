#include "smesh_improve.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_IMPROVE(IDX_T, COUNT_T, GEOM_T)                                   \
    template int mesh_improve<IDX_T, COUNT_T, GEOM_T>(                                               \
            const enum ElemType,                                                                     \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                 \
            const int,                                                                               \
            const ptrdiff_t,                                                                         \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const GEOM_T,                                                                            \
            const GEOM_T,                                                                            \
            const GEOM_T,                                                                            \
            const int,                                                                               \
            const int,                                                                               \
            const int,                                                                               \
            const int,                                                                               \
            ptrdiff_t *,                                                                             \
            IDX_T ***,                                                                               \
            ptrdiff_t *,                                                                             \
            GEOM_T ***,                                                                              \
            uint8_t **,                                                                              \
            uint8_t **,                                                                              \
            GEOM_T ***,                                                                              \
            ptrdiff_t *);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_IMPROVE(i32, i32, f32)
SMESH_EXPLICIT_INSTANTIATE_IMPROVE(i32, i32, f64)
}  // namespace smesh
