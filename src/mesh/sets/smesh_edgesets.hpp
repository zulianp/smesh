#ifndef SMESH_EDGESETS_HPP
#define SMESH_EDGESETS_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename element_idx_t>
int extract_edges_from_edgeset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_edges,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT edge_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT edges);

template <typename idx_t, typename element_idx_t>
int extract_nodeset_from_edgeset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_edges,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT edge_idx, ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out);

}  // namespace smesh

#endif  // SMESH_EDGESETS_HPP
