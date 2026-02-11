#include <cstddef>

#include "smesh_sshex8_mesh.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_MESH(IDX_T, GEOM_T, REF_T)                  \
  template int sshex8_to_standard_hex8_mesh<IDX_T>(                            \
      const int, const ptrdiff_t,                                              \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      IDX_T *SMESH_RESTRICT *const SMESH_RESTRICT);                            \
  template int sshex8_fill_points<IDX_T, GEOM_T>(                              \
      const int, const ptrdiff_t,                                              \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                \
      GEOM_T *SMESH_RESTRICT *const SMESH_RESTRICT);                           \
  template int sshex8_fill_points_1D_map<IDX_T, GEOM_T, REF_T>(                \
      const int, const ptrdiff_t,                                              \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                \
      const REF_T *const SMESH_RESTRICT,                                       \
      GEOM_T *SMESH_RESTRICT *const SMESH_RESTRICT)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_MESH(i16, f32, real_t);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_MESH(i32, f32, real_t);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_MESH(i64, f32, real_t);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_MESH
