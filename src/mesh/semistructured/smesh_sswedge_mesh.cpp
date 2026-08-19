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
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH(i64, f64);
}  // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSWEDGE_MESH
