#include "smesh_build.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(IDX_T, GEOM_T)                   \
  template void mesh_fill_hex8_cube<IDX_T, GEOM_T>(                            \
      int, int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T,           \
      IDX_T *SMESH_RESTRICT *const SMESH_RESTRICT,                             \
      GEOM_T *SMESH_RESTRICT *const SMESH_RESTRICT);                           \
  template void mesh_fill_tri3_square<IDX_T, GEOM_T>(                          \
      int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T,                                \
      IDX_T *SMESH_RESTRICT *const SMESH_RESTRICT,                             \
      GEOM_T *SMESH_RESTRICT *const SMESH_RESTRICT);                           \
  template void mesh_fill_quad4_square<IDX_T, GEOM_T>(                         \
      int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T,                                \
      IDX_T *SMESH_RESTRICT *const SMESH_RESTRICT,                             \
      GEOM_T *SMESH_RESTRICT *const SMESH_RESTRICT)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(i64, f32);

} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD
