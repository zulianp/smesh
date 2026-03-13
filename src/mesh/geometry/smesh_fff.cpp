#include "smesh_fff.impl.hpp"

#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_FFF(FFFType)                                \
  template int tet4_fff_fill<FFFType>(                                         \
      const ptrdiff_t nelements,                                               \
      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const ptrdiff_t stride,                                                  \
      FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);                \
  template int hex8_fff_fill<FFFType>(                                         \
      const ptrdiff_t nelements,                                               \
      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const geom_t qx, const geom_t qy, const geom_t qz,                       \
      const ptrdiff_t stride,                                                  \
      FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff)

SMESH_EXPLICIT_INSTANTIATE_FFF(f32);
SMESH_EXPLICIT_INSTANTIATE_FFF(f16);

#undef SMESH_EXPLICIT_INSTANTIATE_FFF

} // namespace smesh
