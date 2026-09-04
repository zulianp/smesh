#ifndef SMESH_REFINE_HPP
#define SMESH_REFINE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

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

inline bool refine_is_quad_family(const enum ElemType t) {
    return t == QUAD4 || t == QUADSHELL4;
}

inline bool refine_type_supported(const enum ElemType t) {
    return t == HEX8 || t == TET4 || refine_is_tri_family(t) || refine_is_quad_family(t) ||
           t == WEDGE6 || refine_is_edge_family(t) || t == PYRAMID5;
}

struct RefineTypeSet {
    bool          hex{false};
    bool          tet{false};
    bool          tri{false};
    bool          quad{false};
    bool          wedge{false};
    bool          edge{false};
    bool          pyramid{false};
    bool          other{false};
    bool          all_same{true};
    bool          all_supported{true};
    enum ElemType et0{INVALID};
    enum ElemType first_unsupported{INVALID};
    size_t        mixed_block{0};
    enum ElemType mixed_type{INVALID};
};

inline void refine_note_type(RefineTypeSet &s, const size_t b, const enum ElemType t) {
    if (b == 0) {
        s.et0 = t;
    } else if (s.all_same && t != s.et0) {
        s.all_same     = false;
        s.mixed_block  = b;
        s.mixed_type   = t;
    }
    if (!refine_type_supported(t)) {
        s.all_supported = false;
        if (s.first_unsupported == INVALID) {
            s.first_unsupported = t;
        }
    }
    if (t == HEX8) {
        s.hex = true;
    } else if (t == TET4) {
        s.tet = true;
    } else if (refine_is_tri_family(t)) {
        s.tri = true;
    } else if (refine_is_quad_family(t)) {
        s.quad = true;
    } else if (t == WEDGE6) {
        s.wedge = true;
    } else if (refine_is_edge_family(t)) {
        s.edge = true;
    } else if (t == PYRAMID5) {
        s.pyramid = true;
    } else {
        s.other = true;
    }
}

template <typename MeshT>
inline RefineTypeSet refine_scan_mesh(const MeshT &mesh) {
    RefineTypeSet s;
    for (size_t b = 0; b < mesh.n_blocks(); ++b) {
        refine_note_type(s, b, mesh.element_type(static_cast<block_idx_t>(b)));
    }
    return s;
}

inline bool refine_hex_wedge_only(const RefineTypeSet &s) {
    return s.hex && s.wedge && !s.tet && !s.tri && !s.quad && !s.edge && !s.pyramid && !s.other;
}

inline bool refine_quad_family_only(const RefineTypeSet &s) {
    return s.quad && !s.hex && !s.tet && !s.tri && !s.wedge && !s.edge && !s.pyramid && !s.other;
}

/// True for any mix of HEX/TET/WEDGE/PYRAMID (volume types that go through the SS lattice).
/// No surface types (TRI/QUAD/EDGE) and no other types.
inline bool refine_mixed_volume_ss(const RefineTypeSet &s) {
    return !s.tri && !s.quad && !s.edge && !s.other &&
           (s.hex || s.tet || s.wedge || s.pyramid) &&
           !s.all_same;
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

/// Bey 8-split. 10-node numbering: corners 0–3, Exodus edge mids 4–9.
/// Shared by serial `mesh_refine`, MPI `refine_edges_once`, and TET set mappers.
inline constexpr int tet4_refine_pattern[8][4] = {
        {0, 4, 6, 7},
        {4, 1, 5, 8},
        {6, 5, 2, 9},
        {7, 8, 9, 3},
        {4, 5, 6, 8},
        {7, 4, 6, 8},
        {6, 5, 9, 8},
        {7, 6, 9, 8}};
inline constexpr int tet4_refine_edges[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
/// Child tets on coarse face `lfi` (order is the sideset contract).
inline constexpr int tet4_face_child[4][4] = {{0, 1, 3, 5}, {1, 2, 3, 6}, {0, 2, 3, 7}, {0, 1, 2, 4}};
/// Child tets on coarse edge `lei` (order is the edgeset contract).
inline constexpr int tet4_edge_child[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

/// 4-split. 6-node numbering: corners 0–2, Exodus edge mids 3–5.
inline constexpr int tri3_refine_pattern[4][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};
inline constexpr int tri3_refine_edges[3][2]   = {{0, 1}, {1, 2}, {0, 2}};
/// Child tris on coarse edge `lfi`/`lei` (order is the set-mapper contract).
inline constexpr int tri3_edge_child[3][2] = {{0, 1}, {1, 2}, {2, 0}};

/// 2-split. 3-node numbering: endpoints 0–1, midpoint 2.
inline constexpr int edge2_refine_pattern[2][2] = {{0, 2}, {2, 1}};
inline constexpr int edge2_refine_edges[1][2]   = {{0, 1}};

inline bool refine_is_higher_order(const enum ElemType t) {
    switch (t) {
        case TET10:
        case TET15:
        case TET20:
        case TRI6:
        case TRI10:
        case TRISHELL6:
        case TRISHELL10:
        case HEX27:
        case QUAD9:
        case QUADSHELL9:
        case EDGE3:
        case EDGESHELL3:
            return true;
        default:
            return false;
    }
}

inline const char *refine_unsupported_workaround(const enum ElemType t) {
    if (is_semistructured_type(t)) {
        return "refine() does not h-refine SS lattices or return PROTEUS_* types. "
               "Explode to linear (sshex_to_hex8 / sstet_to_tet4 / ssquad_to_quad4 / sswedge_to_wedge6) then refine, "
               "or regenerate SS from the linear parent at a higher L";
    }
    if (t == PYRAMID5) {
        return "A pyramid does not tile into pyramids under uniform edge bisection (tets appear). "
               "Use to_semistructured for an SS lattice. convert_to(TET4) splits each pyramid into 2 tets. "
               "Mixed PYRAMID+TET explode is not a refine() path";
    }
    if (refine_is_higher_order(t)) {
        return "refine() is P1-only; dropping mid-nodes then h-refine is lossy and is not a refine() path. "
               "Use the linear parent, or promote_to after refine()";
    }
    return "Supported types: HEX8, TET4, TRI3, TRISHELL3, QUAD4, QUADSHELL4, WEDGE6, EDGE2, EDGESHELL2";
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
            "HEX+WEDGE and QUAD4+QUADSHELL4 are supported; HEX+TET is not "
            "(Bey vs Kuhn). PYRAMID stays SS-only (to_semistructured). "
            "Split other combinations.\n",
            type_to_string(a),
            b,
            type_to_string(tb));
}

inline void refine_print_mixed_hex_tet() {
    fprintf(stderr,
            "refine: mixed HEX+TET is not supported (TET refine is Bey 8-split, not Kuhn SS explode). "
            "Workaround: to_semistructured then sstet_to_tet4 (changes TET child layout), "
            "or split blocks and refine separately.\n");
}

inline void refine_print_mixed_hex_quad() {
    fprintf(stderr,
            "refine: mixed HEX/WEDGE+QUAD is not supported (no shared SS lattice for volume+surface). "
            "Split by family and refine separately.\n");
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
