#ifndef SMESH_SSQUAD4_MESH_HPP
#define SMESH_SSQUAD4_MESH_HPP

#include "smesh_base.hpp"

namespace smesh {

template <typename idx_t>
int ssquad4_to_standard_quad4_mesh(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements);

} // namespace smesh

#endif // SMESH_SSQUAD4_MESH_HPP
