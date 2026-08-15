#include "smesh_types.hpp"
#include "smesh_ssquad4_mesh.impl.hpp"

namespace smesh {
template int ssquad4_to_standard_quad4_mesh<i32>(
    int, ptrdiff_t,
    const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    i32 *SMESH_RESTRICT *const SMESH_RESTRICT);
template int ssquad4_to_standard_quad4_mesh<i64>(
    int, ptrdiff_t,
    const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    i64 *SMESH_RESTRICT *const SMESH_RESTRICT);

template int ssquad4_fill_points<i32, geom_t>(const int,
                                              const ptrdiff_t,
                                              const int,
                                              const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                              const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                              geom_t *SMESH_RESTRICT *const SMESH_RESTRICT);
template int ssquad4_fill_points<i64, geom_t>(const int,
                                              const ptrdiff_t,
                                              const int,
                                              const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                              const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                              geom_t *SMESH_RESTRICT *const SMESH_RESTRICT);
} // namespace smesh

