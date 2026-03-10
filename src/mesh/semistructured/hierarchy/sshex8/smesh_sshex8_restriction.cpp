#include "smesh_sshex8_restriction.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(IDX_T, REAL_T)           \
  template int sshex8_hierarchical_restriction<IDX_T, REAL_T>(                 \
      int level, const ptrdiff_t nelements,                                    \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const u16 *const SMESH_RESTRICT element_to_node_incidence_count,         \
      const int vec_size, const REAL_T *const SMESH_RESTRICT from,             \
      REAL_T *const SMESH_RESTRICT to);

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICT(IDX_T, REAL_T)           \
  template int sshex8_restrict<IDX_T, REAL_T>(                              \
      const ptrdiff_t nelements, const int from_level, const int from_level_stride, \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements, const int to_level, const int to_level_stride, \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements, const int vec_size, const REAL_T *const SMESH_RESTRICT from, \
      REAL_T *const SMESH_RESTRICT to);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i64, f32);

SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i64, f64);

SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICT(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICT(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICT(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICT(i64, f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION
#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICT