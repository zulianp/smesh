#include "smesh_sshex8_restriction.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(IDX_T, REAL_T)           \
  template int sshex8_hierarchical_restriction<IDX_T, REAL_T>(                 \
      int level, const ptrdiff_t nelements,                                    \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const u16 *const SMESH_RESTRICT element_to_node_incidence_count,         \
      const int vec_size, const REAL_T *const SMESH_RESTRICT from,             \
      REAL_T *const SMESH_RESTRICT to);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i64, f32);

SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION(i64, f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_RESTRICTION
