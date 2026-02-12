#include "smesh_ops.impl.hpp"
#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_BARYCENTERS(IDX_T, GEOM_T)                   \
  template int barycenters<IDX_T, GEOM_T>(                                      \
      const int, const ptrdiff_t,                                               \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                  \
      const int,                                                                \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT)

SMESH_EXPLICIT_INSTANTIATE_BARYCENTERS(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_BARYCENTERS(i64, f32);

#undef SMESH_EXPLICIT_INSTANTIATE_BARYCENTERS

} // namespace smesh

