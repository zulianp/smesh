#ifndef SMESH_SIDESETS_HPP
#define SMESH_SIDESETS_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename element_idx_t>
int extract_nodeset_from_sideset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT side_idx, ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out);

template <typename idx_t, typename element_idx_t>
int extract_nodeset_from_sidesets(
    ptrdiff_t n_sidesets, const enum ElemType element_type[],
    idx_t **const SMESH_RESTRICT elems[], const ptrdiff_t n_surf_elements[],
    const element_idx_t *const SMESH_RESTRICT parent_element[],
    const int16_t *const SMESH_RESTRICT side_idx[], ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out);

template <typename idx_t, typename element_idx_t>
int extract_surface_from_sideset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT side_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT sides);
} // namespace smesh

#endif // SMESH_SIDESETS_HPP
