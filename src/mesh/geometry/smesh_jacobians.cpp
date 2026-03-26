#include "smesh_jacobians.impl.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_JACOBIANS(AdjugateType, DeterminantType)    \
  template int tet4_adjugate_fill<AdjugateType, DeterminantType>(              \
      const ptrdiff_t nelements,                                               \
      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const ptrdiff_t stride,                                                  \
      AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,       \
      DeterminantType *const SMESH_RESTRICT determinant);                      \
  template int hex8_adjugate_fill<AdjugateType, DeterminantType>(              \
      const ptrdiff_t nelements,                                               \
      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const geom_t qx, const geom_t qy, const geom_t qz,                       \
      const ptrdiff_t stride,                                                  \
      AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,       \
      DeterminantType *const SMESH_RESTRICT determinant);                      \
  template int sshex8_macro_adjugate_fill<AdjugateType, DeterminantType>(      \
      const int level,                                                         \
      const ptrdiff_t nelements,                                               \
      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const ptrdiff_t stride,                                                  \
      AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,       \
      DeterminantType *const SMESH_RESTRICT determinant);

SMESH_EXPLICIT_INSTANTIATE_JACOBIANS(f32, f32);
SMESH_EXPLICIT_INSTANTIATE_JACOBIANS(f16, f32);

#undef SMESH_EXPLICIT_INSTANTIATE_JACOBIANS

} // namespace smesh
