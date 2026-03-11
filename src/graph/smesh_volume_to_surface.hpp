#ifndef SMESH_VOLUME_TO_SURFACE_HPP
#define SMESH_VOLUME_TO_SURFACE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int extract_skin_sideset_from_n2e(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx,
    ptrdiff_t *SMESH_RESTRICT n_surf_elements,
    element_idx_t **SMESH_RESTRICT parent_element,
    i16 **SMESH_RESTRICT side_idx);

} // namespace smesh

#endif // SMESH_VOLUME_TO_SURFACE_HPP
