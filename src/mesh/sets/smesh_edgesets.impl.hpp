#ifndef SMESH_EDGESETS_IMPL_HPP
#define SMESH_EDGESETS_IMPL_HPP

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_edgesets.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_sort.hpp"

namespace smesh {

template <typename idx_t, typename element_idx_t>
int extract_edges_from_edgeset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_edges,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT edge_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT edges) {
    LocalEdgeTable let;
    if (let.fill(element_type) != SMESH_SUCCESS) {
        LocalEdgeTable::report_unsupported("extract_edges_from_edgeset", element_type);
        return SMESH_FAILURE;
    }

    for (ptrdiff_t i = 0; i < n_edges; ++i) {
        const ptrdiff_t e  = parent_element[i];
        const int       le = edge_idx[i];
        const int       nn = let.nnxe_edge[le];
        for (int n = 0; n < nn; ++n) {
            edges[n][i] = elems[let(le, n)][e];
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int extract_nodeset_from_edgeset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_edges,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT edge_idx, ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out) {
    LocalEdgeTable let;
    if (let.fill(element_type) != SMESH_SUCCESS) {
        LocalEdgeTable::report_unsupported("extract_nodeset_from_edgeset", element_type);
        return SMESH_FAILURE;
    }

    ptrdiff_t n = 0;
    for (ptrdiff_t i = 0; i < n_edges; ++i) {
        n += let.nnxe_edge[edge_idx[i]];
    }
    idx_t *nodes = (idx_t *)SMESH_ALLOC((size_t)n * sizeof(idx_t));

    ptrdiff_t write = 0;
    for (ptrdiff_t i = 0; i < n_edges; ++i) {
        const ptrdiff_t e  = parent_element[i];
        const int       le = edge_idx[i];
        const int       nn = let.nnxe_edge[le];
        for (int k = 0; k < nn; ++k) {
            nodes[write++] = elems[let(le, k)][e];
        }
    }

    *n_nodes_out = (ptrdiff_t)sort_and_unique(nodes, (size_t)n);
    *nodes_out   = (idx_t *)SMESH_REALLOC(nodes, (size_t)(*n_nodes_out) * sizeof(idx_t));
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif  // SMESH_EDGESETS_IMPL_HPP
