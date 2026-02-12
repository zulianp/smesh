#include "smesh_sfc.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SFC(GEOM_T)                                 \
  template int encode_morton3<GEOM_T>(                                         \
      const ptrdiff_t, const GEOM_T *const SMESH_RESTRICT,                     \
      const GEOM_T *const SMESH_RESTRICT, const GEOM_T *const SMESH_RESTRICT,  \
      u32 *const SMESH_RESTRICT);                                              \
  template int encode_hilbert3<GEOM_T>(                                        \
      const ptrdiff_t, const GEOM_T *const SMESH_RESTRICT,                     \
      const GEOM_T *const SMESH_RESTRICT, const GEOM_T *const SMESH_RESTRICT,  \
      u32 *const SMESH_RESTRICT);                                              \
  template int encode_cartesian3<GEOM_T>(                                      \
      const ptrdiff_t, const GEOM_T *const SMESH_RESTRICT,                     \
      const GEOM_T *const SMESH_RESTRICT, const GEOM_T *const SMESH_RESTRICT,  \
      int, int, int, u32 *const SMESH_RESTRICT)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_SFC(f32);
SMESH_EXPLICIT_INSTANTIATE_SFC(f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SFC
