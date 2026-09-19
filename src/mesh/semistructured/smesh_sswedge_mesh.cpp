#include <cstddef>

#include "smesh_sswedge_mesh.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(IDX_T, GEOM_T)            \
    template int sswedge_fill_points<IDX_T, GEOM_T>(                      \
            const int,                                                    \
            const ptrdiff_t,                                              \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,      \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,     \
            GEOM_T *SMESH_RESTRICT *const SMESH_RESTRICT)

namespace smesh {
template int sswedge_to_standard_wedge6_mesh<i16>(
        const int,
        const ptrdiff_t,
        const i16 *const SMESH_RESTRICT *const SMESH_RESTRICT,
        i16 *SMESH_RESTRICT *const SMESH_RESTRICT);
template int sswedge_to_standard_wedge6_mesh<i32>(
        const int,
        const ptrdiff_t,
        const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
        i32 *SMESH_RESTRICT *const SMESH_RESTRICT);
template int sswedge_to_standard_wedge6_mesh<i64>(
        const int,
        const ptrdiff_t,
        const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
        i64 *SMESH_RESTRICT *const SMESH_RESTRICT);

SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i64, f64);
}  // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH
