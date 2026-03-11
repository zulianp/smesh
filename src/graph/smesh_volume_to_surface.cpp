#include "smesh_volume_to_surface.impl.hpp"

// TODO: explicit instantiation
#define SMESH_EXPLICIT_INSTANTIATE_VOLUME_TO_SURFACE(T, C, E)                    \
  template int extract_skin_sideset_from_n2e<T, C, E>(                           \
      const ptrdiff_t, const ptrdiff_t, const enum ElemType,                      \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT,                        \
      const C *const SMESH_RESTRICT,                                              \
      const E *const SMESH_RESTRICT,                                              \
      ptrdiff_t *, E **, i16 **);                                                 \

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_VOLUME_TO_SURFACE(i32, i32, i32);
SMESH_EXPLICIT_INSTANTIATE_VOLUME_TO_SURFACE(i64, i64, i64);
} // namespace smesh
