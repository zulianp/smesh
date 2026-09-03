#ifndef SMESH_REFINE_HPP
#define SMESH_REFINE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

#include <stddef.h>
#include <stdio.h>

namespace smesh {

/// Unstructured `refine()` contract: HEX8, TET4, TRI3/TRISHELL3, QUAD4/QUADSHELL4,
/// WEDGE6, EDGE2/EDGESHELL2 (see REFINE.md).
inline bool refine_is_tri_family(const enum ElemType t) {
    return t == TRI3 || t == TRISHELL3;
}

inline bool refine_is_edge_family(const enum ElemType t) {
    return t == EDGE2 || t == EDGESHELL2;
}

inline bool refine_type_supported(const enum ElemType t) {
    return t == HEX8 || t == TET4 || refine_is_tri_family(t) || t == QUAD4 || t == QUADSHELL4 ||
           t == WEDGE6 || refine_is_edge_family(t);
}

/// Child factor for one edge-midpoint `mesh_refine` / `refine_edges_once` step.
inline int refine_edge_midpoint_factor(const enum ElemType t) {
    if (refine_is_tri_family(t)) {
        return 4;
    }
    if (refine_is_edge_family(t)) {
        return 2;
    }
    return 8;
}

inline const char *refine_unsupported_workaround(const enum ElemType t) {
    if (is_semistructured_type(t)) {
        return "Explode to linear (sshex_to_hex8 / sstet_to_tet4 / ssquad_to_quad4 / sswedge_to_wedge6) then refine, "
               "or regenerate SS from the linear parent at a higher L";
    }
    switch (t) {
        case PYRAMID5:
            return "Use to_semistructured (unstructured PYRAMID refine is not implemented)";
        case TET10:
        case TET15:
        case TET20:
        case TRI6:
        case TRI10:
        case TRISHELL6:
        case HEX27:
        case QUAD9:
        case QUADSHELL9:
            return "refine() is P1-only; demote is not implemented";
        default:
            return "Supported types: HEX8, TET4, TRI3, TRISHELL3, QUAD4, QUADSHELL4, WEDGE6, EDGE2, EDGESHELL2";
    }
}

inline void refine_print_unsupported(const enum ElemType t) {
    fprintf(stderr,
            "refine: unsupported element type %s. %s\n",
            type_to_string(t),
            refine_unsupported_workaround(t));
}

inline void refine_print_mixed_types(const enum ElemType a, const size_t b, const enum ElemType tb) {
    fprintf(stderr,
            "refine: mixed block types are not supported (block 0 is %s, block %zu is %s). "
            "Convert or split first.\n",
            type_to_string(a),
            b,
            type_to_string(tb));
}

template <typename idx_t, typename count_t,
          typename geom_t>
int mesh_refine(
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT coarse_elements,
    const int spatial_dim, const ptrdiff_t n_coarse_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT coarse_points,
    const count_t *const SMESH_RESTRICT n2n_ptr,
    const idx_t *const SMESH_RESTRICT n2n_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT refined_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT refined_points);

} // namespace smesh

#endif // SMESH_REFINE_HPP
