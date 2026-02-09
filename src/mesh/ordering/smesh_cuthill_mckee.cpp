#include "smesh_cuthill_mckee.impl.hpp"
#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_CUTHILL_MCKEE(COUNT_T, IDX_T)                \
  template int eccentricity<COUNT_T, IDX_T>(                                    \
      const ptrdiff_t, const COUNT_T *const SMESH_RESTRICT,                     \
      const IDX_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT);          \
  template int cuthill_mckee<COUNT_T, IDX_T>(                                   \
      const ptrdiff_t, const COUNT_T *const SMESH_RESTRICT,                     \
      const IDX_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT)

SMESH_EXPLICIT_INSTANTIATE_CUTHILL_MCKEE(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_CUTHILL_MCKEE(i64, i32);
SMESH_EXPLICIT_INSTANTIATE_CUTHILL_MCKEE(i64, i64);

#undef SMESH_EXPLICIT_INSTANTIATE_CUTHILL_MCKEE

} // namespace smesh
