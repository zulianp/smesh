#include "smesh_restriction.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_RESTRICTION(IDX_T, T)                       \
  template int hierarchical_restriction<IDX_T, T>(                             \
      const enum ElemType, const enum ElemType, const ptrdiff_t,               \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const u16 *const SMESH_RESTRICT e2n_count, const int vec_size,           \
      const T *const SMESH_RESTRICT from, T *const SMESH_RESTRICT to)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_RESTRICTION(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_RESTRICTION(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_RESTRICTION(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_RESTRICTION(i64, f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_RESTRICTION
