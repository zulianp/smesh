#include "smesh_build.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(IDX_T, GEOM_T)                   \
  template void mesh_fill_hex8_reference_cube<IDX_T, GEOM_T>(                  \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,              \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT points);              \
  template void mesh_fill_hex8_cube<IDX_T, GEOM_T>(                            \
      int, int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T,           \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                     \
  template void mesh_fill_tet4_cube<IDX_T, GEOM_T>(                            \
      int, int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T,           \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                     \
  template void mesh_fill_tri3_square<IDX_T, GEOM_T>(                          \
      int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T,                                \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                     \
  template void mesh_fill_quad4_square<IDX_T, GEOM_T>(                         \
      int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T,                                \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                     \
  template void mesh_fill_hex8_bidomain_cube<IDX_T, GEOM_T>(                   \
      int, int, int, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T, GEOM_T, int,      \
      IDX_T, IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                     \
  template void mesh_fill_hex8_checkerboard_cube<IDX_T, GEOM_T>(               \
      const int nx, const int ny, const int nz, const GEOM_T xmin,             \
      const GEOM_T ymin, const GEOM_T zmin, const GEOM_T xmax,                 \
      const GEOM_T ymax, const GEOM_T zmax,                                    \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT white_elements,        \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT black_elements,        \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT points);              \
  template void mesh_fill_quad4_ring<IDX_T, GEOM_T>(                           \
      const GEOM_T inner_radius, const GEOM_T outer_radius,                    \
      const ptrdiff_t nlayers, const ptrdiff_t nelements,                      \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,              \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT points)

namespace smesh {
// SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD(i64, f32);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_MESH_BUILD
