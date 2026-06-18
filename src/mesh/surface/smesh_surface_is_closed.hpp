#ifndef SMESH_SURFACE_IS_CLOSED_HPP
#define SMESH_SURFACE_IS_CLOSED_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {

bool surface_is_closed(
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_nodes, const count_t *const SMESH_RESTRICT n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx);

} // namespace smesh

#endif // SMESH_SURFACE_IS_CLOSED_HPP
