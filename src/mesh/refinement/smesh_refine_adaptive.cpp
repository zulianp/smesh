#include "smesh_refine_adaptive.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(IDX_T, COUNT_T, GEOM_T)                                             \
    template int mesh_refine_adaptive<IDX_T, COUNT_T, GEOM_T>(                                                         \
            const ptrdiff_t, const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const int, const ptrdiff_t,     \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,              \
            const ptrdiff_t, const int, IDX_T ***, ptrdiff_t *, GEOM_T ***, ptrdiff_t *, COUNT_T **, IDX_T **,        \
            IDX_T **, COUNT_T **, IDX_T **, IDX_T **, IDX_T **, IDX_T **, IDX_T **, ptrdiff_t *, IDX_T **,           \
            const int *);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(i32, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(i64, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(i64, i64, f32);
SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(i32, i32, f64);
SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(i64, i32, f64);
SMESH_EXPLICIT_INSTANTIATE_REFINE_ADAPTIVE(i64, i64, f64);
}  // namespace smesh
