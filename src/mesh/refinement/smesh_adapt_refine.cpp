#include "smesh_adapt_refine.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_ADAPT(IDX_T, COUNT_T, GEOM_T)                                     \
    template int mesh_adapt_refine<IDX_T, COUNT_T, GEOM_T>(                                          \
            const enum ElemType,                                                                     \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                 \
            const int,                                                                               \
            const ptrdiff_t,                                                                         \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            const GEOM_T *const SMESH_RESTRICT,                                                      \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            const GEOM_T,                                                                            \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            const int,                                                                               \
            ptrdiff_t *,                                                                             \
            IDX_T ***,                                                                               \
            ptrdiff_t *,                                                                             \
            GEOM_T ***,                                                                              \
            ptrdiff_t **,                                                                            \
            COUNT_T **,                                                                              \
            IDX_T **,                                                                                \
            IDX_T **,                                                                                \
            IDX_T **,                                                                                \
            int (*)(void *, GEOM_T, GEOM_T, GEOM_T, IDX_T, IDX_T, const IDX_T *, const IDX_T *, ptrdiff_t, GEOM_T *, GEOM_T *, GEOM_T *), \
            void (*)(void *, ptrdiff_t, IDX_T **, ptrdiff_t, GEOM_T **, const IDX_T *, const IDX_T *), \
            void *);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_ADAPT(i32, i32, f32)
SMESH_EXPLICIT_INSTANTIATE_ADAPT(i64, i32, f32)
SMESH_EXPLICIT_INSTANTIATE_ADAPT(i64, i64, f32)
SMESH_EXPLICIT_INSTANTIATE_ADAPT(i32, i32, f64)
SMESH_EXPLICIT_INSTANTIATE_ADAPT(i64, i32, f64)
SMESH_EXPLICIT_INSTANTIATE_ADAPT(i64, i64, f64)
}  // namespace smesh
