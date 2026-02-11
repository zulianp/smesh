#include "smesh_types.hpp"
#include "smesh_ssquad4_mesh.impl.hpp"

namespace smesh {
template int ssquad4_to_standard_quad4_mesh<i32>(
    int, ptrdiff_t,
    const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    i32 *SMESH_RESTRICT *const SMESH_RESTRICT);
} // namespace smesh

