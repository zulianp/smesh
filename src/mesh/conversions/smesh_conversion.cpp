#include "smesh_conversion.impl.hpp"

#include "smesh_types.hpp"
// TODO: explicit instantiation
#define SMESH_EXPLICIT_INSTANTIATE_CONVERSIONS(IDX_T)                          \
  template void mesh_hex8_to_6x_tet4<IDX_T>(                                   \
      const ptrdiff_t,                                                         \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                      \
  template void mesh_tet15_to_4x_hex8<IDX_T>(                                  \
      const ptrdiff_t,                                                         \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                      \
  template void mesh_wedge6_to_3x_tet4<IDX_T>(                                 \
      const ptrdiff_t,                                                         \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_CONVERSIONS(i16);
SMESH_EXPLICIT_INSTANTIATE_CONVERSIONS(i32);
SMESH_EXPLICIT_INSTANTIATE_CONVERSIONS(i64);
} // namespace smesh
