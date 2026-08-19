#ifndef SMESH_SSMIXED_HEX_DOMINANT_IMPL_HPP
#define SMESH_SSMIXED_HEX_DOMINANT_IMPL_HPP

#include "smesh_sspyramid.hpp"
#include "smesh_sswedge.hpp"

#include <map>
#include <utility>

namespace smesh {

static int hexdom_wedge_conn[6][3] = {{1, 2, 3}, {0, 2, 4}, {0, 1, 5}, {0, 4, 5}, {1, 3, 5}, {2, 3, 4}};
static int hexdom_pyr_nneigh[5]    = {3, 3, 3, 3, 4};
static int hexdom_pyr_conn[5][4]   = {{1, 3, 4, -1}, {0, 2, 4, -1}, {1, 3, 4, -1}, {0, 2, 4, -1}, {0, 1, 2, 3}};
static const int hexdom_wedge_xyz[6][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
static const int hexdom_pyr_ijk[5][3]   = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}};

static int hexdom_n_macro(const enum ElemType f) {
    if (f == HEX8) {
        return 8;
    }
    if (f == TET4) {
        return 4;
    }
    if (f == WEDGE6) {
        return 6;
    }
    return 5;
}

static int hexdom_nsides(const enum ElemType f) {
    if (f == HEX8) {
        return 6;
    }
    if (f == TET4) {
        return 4;
    }
    return 5;
}

static int hexdom_side_nnxs(const enum ElemType f, const int s) {
    if (f == HEX8) {
        return 4;
    }
    if (f == TET4) {
        return 3;
    }
    if (f == WEDGE6) {
        return s < 3 ? 4 : 3;
    }
    return s < 4 ? 3 : 4;
}

static int hexdom_lidx(const enum ElemType f, const int L, const int x, const int y, const int z) {
    if (f == HEX8) {
        return sshex8_lidx(L, x, y, z);
    }
    if (f == TET4) {
        return sstet4_lidx(L, x, y, z);
    }
    if (f == WEDGE6) {
        return sswedge_lidx(L, x, y, z);
    }
    return sspyramid_lidx(L, x, y, z);
}

static int hexdom_pyr_corner_lidx(const int L, const int c) {
    if (c == 4) {
        return sspyramid_lidx(L, 0, 0, L);
    }
    return sspyramid_lidx(L, hexdom_pyr_ijk[c][0] * L, hexdom_pyr_ijk[c][1] * L, hexdom_pyr_ijk[c][2] * L);
}

template <typename idx_t, typename count_t, typename element_idx_t>
static int mixed_hex_dominant_build_edge_graph_from_n2e(
        const ptrdiff_t                                                                n_blocks,
        const enum ElemType *const                                                     block_types,
        const ptrdiff_t *const                                                         n_elements,
        const ptrdiff_t                                                                nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
        const count_t *const SMESH_RESTRICT                                            n2eptr,
        const element_idx_t *const SMESH_RESTRICT                                      elindex,
        const block_idx_t *const SMESH_RESTRICT                                        block_number,
        count_t                                                                      **out_rowptr,
        idx_t                                                                        **out_colidx) {
    count_t *rowptr = (count_t *)SMESH_ALLOC((nnodes + 1) * sizeof(count_t));
    rowptr[0]       = 0;
#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        idx_t   n2nbuff[2048];
        count_t ebegin  = n2eptr[node];
        count_t eend    = n2eptr[node + 1];
        idx_t   nneighs = 0;
        for (count_t e = ebegin; e < eend; ++e) {
            const block_idx_t   b    = block_number[e];
            const element_idx_t eidx = elindex[e];
            const enum ElemType et   = block_types[b];
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];
            SMESH_UNUSED(n_elements);
            SMESH_UNUSED(n_blocks);
            if (et == HEX8) {
                int lidx = -1;
                for (int i = 0; i < 8; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                SMESH_ASSERT(lidx != -1);
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[mixed_hex8_edge_connectivity[lidx][d]][eidx];
                }
            } else if (et == TET4) {
                int lidx = -1;
                for (int i = 0; i < 4; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                SMESH_ASSERT(lidx != -1);
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[mixed_tet4_edge_connectivity[lidx][d]][eidx];
                }
            } else if (et == WEDGE6) {
                int lidx = -1;
                for (int i = 0; i < 6; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                SMESH_ASSERT(lidx != -1);
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[hexdom_wedge_conn[lidx][d]][eidx];
                }
            } else {
                int lidx = -1;
                for (int i = 0; i < 5; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                SMESH_ASSERT(lidx != -1);
                const int nn = hexdom_pyr_nneigh[lidx];
                for (int d = 0; d < nn; d++) {
                    n2nbuff[nneighs++] = block_elems[hexdom_pyr_conn[lidx][d]][eidx];
                }
            }
        }
        nneighs          = sort_and_unique(n2nbuff, nneighs);
        rowptr[node + 1] = nneighs;
    }
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        rowptr[node + 1] += rowptr[node];
    }
    idx_t *colidx = (idx_t *)SMESH_ALLOC(rowptr[nnodes] * sizeof(idx_t));
#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        idx_t   n2nbuff[2048];
        count_t ebegin  = n2eptr[node];
        count_t eend    = n2eptr[node + 1];
        idx_t   nneighs = 0;
        for (count_t e = ebegin; e < eend; ++e) {
            const block_idx_t   b    = block_number[e];
            const element_idx_t eidx = elindex[e];
            const enum ElemType et   = block_types[b];
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];
            if (et == HEX8) {
                int lidx = -1;
                for (int i = 0; i < 8; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[mixed_hex8_edge_connectivity[lidx][d]][eidx];
                }
            } else if (et == TET4) {
                int lidx = -1;
                for (int i = 0; i < 4; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[mixed_tet4_edge_connectivity[lidx][d]][eidx];
                }
            } else if (et == WEDGE6) {
                int lidx = -1;
                for (int i = 0; i < 6; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[hexdom_wedge_conn[lidx][d]][eidx];
                }
            } else {
                int lidx = -1;
                for (int i = 0; i < 5; ++i) {
                    if (block_elems[i][eidx] == static_cast<idx_t>(node)) {
                        lidx = i;
                        break;
                    }
                }
                const int nn = hexdom_pyr_nneigh[lidx];
                for (int d = 0; d < nn; d++) {
                    n2nbuff[nneighs++] = block_elems[hexdom_pyr_conn[lidx][d]][eidx];
                }
            }
        }
        nneighs = sort_and_unique(n2nbuff, nneighs);
        for (idx_t i = 0; i < nneighs; ++i) {
            colidx[rowptr[node] + i] = n2nbuff[i];
        }
    }
    *out_rowptr = rowptr;
    *out_colidx = colidx;
    return 0;
}

template <typename idx_t>
static void hexdom_write_wedge_edges(const int                                               L,
                                     const ptrdiff_t                                         e,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                     const idx_t *const                                      rowptr,
                                     const idx_t *const                                      colidx,
                                     const idx_t *const                                      edge_idx,
                                     const int *const                                        lagr,
                                     int *const *const                                       coords,
                                     const idx_t                                             index_base,
                                     const int                                               nxedge,
                                     idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements) {
    idx_t nodes[6];
    for (int d = 0; d < 6; d++) {
        nodes[d] = m_elements[d][e];
    }
    for (int d1 = 0; d1 < 6; d1++) {
        const idx_t        node1     = nodes[d1];
        const idx_t *const columns   = &colidx[rowptr[node1]];
        const idx_t *const edge_view = &edge_idx[rowptr[node1]];
        idx_t              g_neigh[3];
        for (int k = 0; k < 3; k++) {
            g_neigh[k] = nodes[hexdom_wedge_conn[d1][k]];
        }
        idx_t offsets[3];
        find<3>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);
        for (int d2 = 0; d2 < 3; d2++) {
            const idx_t node2 = g_neigh[d2];
            if (node1 > node2) {
                continue;
            }
            const int lid1 = lagr[d1];
            const int lid2 = lagr[hexdom_wedge_conn[d1][d2]];
            int       P1[3], P2[3];
            for (int d = 0; d < 3; d++) {
                P1[d] = coords[d][lid1];
                P2[d] = coords[d][lid2];
            }
            const idx_t edge_start = index_base + edge_view[offsets[d2]] * static_cast<idx_t>(nxedge);
            for (int t = 1; t < L; ++t) {
                const int xi          = (P1[0] * (L - t) + P2[0] * t) / L;
                const int yi          = (P1[1] * (L - t) + P2[1] * t) / L;
                const int zi          = (P1[2] * (L - t) + P2[2] * t) / L;
                elements[sswedge_lidx(L, xi, yi, zi)][e] = edge_start + (t - 1);
            }
        }
    }
}

template <typename idx_t>
static void hexdom_write_pyr_edges(const int                                               L,
                                   const ptrdiff_t                                         e,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                   const idx_t *const                                      rowptr,
                                   const idx_t *const                                      colidx,
                                   const idx_t *const                                      edge_idx,
                                   const int *const                                        lagr,
                                   int *const *const                                       coords,
                                   const idx_t                                             index_base,
                                   const int                                               nxedge,
                                   idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements) {
    idx_t nodes[5];
    for (int d = 0; d < 5; d++) {
        nodes[d] = m_elements[d][e];
    }
    for (int d1 = 0; d1 < 5; d1++) {
        const idx_t        node1     = nodes[d1];
        const idx_t *const columns   = &colidx[rowptr[node1]];
        const idx_t *const edge_view = &edge_idx[rowptr[node1]];
        const int          nn        = hexdom_pyr_nneigh[d1];
        idx_t              g_neigh[4];
        idx_t              offsets[4];
        for (int k = 0; k < nn; k++) {
            g_neigh[k] = nodes[hexdom_pyr_conn[d1][k]];
        }
        if (nn == 3) {
            find<3>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);
        } else {
            find<4>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);
        }
        for (int d2 = 0; d2 < nn; d2++) {
            const idx_t node2 = g_neigh[d2];
            if (node1 > node2) {
                continue;
            }
            const int lid1 = lagr[d1];
            const int lid2 = lagr[hexdom_pyr_conn[d1][d2]];
            int       P1[3], P2[3];
            for (int d = 0; d < 3; d++) {
                P1[d] = coords[d][lid1];
                P2[d] = coords[d][lid2];
            }
            const idx_t edge_start = index_base + edge_view[offsets[d2]] * static_cast<idx_t>(nxedge);
            for (int t = 1; t < L; ++t) {
                const int xi = (P1[0] * (L - t) + P2[0] * t) / L;
                const int yi = (P1[1] * (L - t) + P2[1] * t) / L;
                const int zi = (P1[2] * (L - t) + P2[2] * t) / L;
                elements[sspyramid_lidx(L, xi, yi, zi)][e] = edge_start + (t - 1);
            }
        }
    }
}

template <typename idx_t>
static void hexdom_index_tri_face(const enum ElemType                                    family,
                                  const int                                               L,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                  const int *const                                        local_side_table,
                                  const int *const                                        lagr,
                                  int *const *const                                       coords,
                                  const idx_t                                             off,
                                  const ptrdiff_t                                         e,
                                  const int                                               f,
                                  idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements) {
    static const int nnxs   = LocalSideTable::MAX_NUM_NODES_PER_SIDE;
    int              argmin = 0;
    idx_t            valmin = m_elements[local_side_table[f * nnxs + 0]][e];
    for (int i = 1; i < 3; i++) {
        const idx_t temp = m_elements[local_side_table[f * nnxs + i]][e];
        if (temp < valmin) {
            argmin = i;
            valmin = temp;
        }
    }
    int lst_o = argmin;
    int lst_u = (argmin + 1) % 3;
    int lst_v = (argmin + 2) % 3;
    if (m_elements[local_side_table[f * nnxs + lst_u]][e] > m_elements[local_side_table[f * nnxs + lst_v]][e]) {
        const int temp = lst_v;
        lst_v          = lst_u;
        lst_u          = temp;
    }
    const int lidx_o = lagr[local_side_table[f * nnxs + lst_o]];
    const int lidx_u = lagr[local_side_table[f * nnxs + lst_u]];
    const int lidx_v = lagr[local_side_table[f * nnxs + lst_v]];
    int       Po[3], Pu[3], Pv[3];
    for (int d = 0; d < 3; d++) {
        Po[d] = coords[d][lidx_o];
        Pu[d] = coords[d][lidx_u];
        Pv[d] = coords[d][lidx_v];
    }
    int local_offset = 0;
    for (int t = 1; t <= L - 2; ++t) {
        for (int s = 1; s <= L - 1 - t; ++s) {
            const int w    = L - s - t;
            const int xi   = (Po[0] * w + Pu[0] * s + Pv[0] * t) / L;
            const int yi   = (Po[1] * w + Pu[1] * s + Pv[1] * t) / L;
            const int zi   = (Po[2] * w + Pu[2] * s + Pv[2] * t) / L;
            const int pidx = hexdom_lidx(family, L, xi, yi, zi);
            elements[pidx][e] = off + local_offset++;
        }
    }
}

template <typename idx_t>
static void hexdom_index_quad_face(const enum ElemType                                    family,
                                   const int                                               L,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                   const int *const                                        local_side_table,
                                   const int *const                                        lagr,
                                   int *const *const                                       coords,
                                   const idx_t                                             off,
                                   const ptrdiff_t                                         e,
                                   const int                                               f,
                                   idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements) {
    static const int nnxs   = LocalSideTable::MAX_NUM_NODES_PER_SIDE;
    int              argmin = 0;
    idx_t            valmin = m_elements[local_side_table[f * nnxs + 0]][e];
    for (int i = 1; i < 4; i++) {
        const idx_t temp = m_elements[local_side_table[f * nnxs + i]][e];
        if (temp < valmin) {
            argmin = i;
            valmin = temp;
        }
    }
    int lst_o = argmin;
    int lst_u = (lst_o + 1) % 4;
    int lst_v = (lst_o + 3) % 4;
    if (m_elements[local_side_table[f * nnxs + lst_u]][e] > m_elements[local_side_table[f * nnxs + lst_v]][e]) {
        const int temp = lst_v;
        lst_v          = lst_u;
        lst_u          = temp;
    }
    const int lidx_o = lagr[local_side_table[f * nnxs + lst_o]];
    const int lidx_u = lagr[local_side_table[f * nnxs + lst_u]];
    const int lidx_v = lagr[local_side_table[f * nnxs + lst_v]];
    int       Po[3], Pu[3], Pv[3];
    for (int d = 0; d < 3; d++) {
        Po[d] = coords[d][lidx_o];
        Pu[d] = coords[d][lidx_u];
        Pv[d] = coords[d][lidx_v];
    }
    int local_offset = 0;
    for (int t = 1; t < L; ++t) {
        for (int s = 1; s < L; ++s) {
            const int xi   = Po[0] + ((Pu[0] - Po[0]) * s + (Pv[0] - Po[0]) * t) / L;
            const int yi   = Po[1] + ((Pu[1] - Po[1]) * s + (Pv[1] - Po[1]) * t) / L;
            const int zi   = Po[2] + ((Pu[2] - Po[2]) * s + (Pv[2] - Po[2]) * t) / L;
            const int pidx = hexdom_lidx(family, L, xi, yi, zi);
            elements[pidx][e] = off + local_offset++;
        }
    }
}

struct HexDomFaceKey {
    idx_t n[4];
    bool  operator<(const HexDomFaceKey &o) const {
        for (int i = 0; i < 4; ++i) {
            if (n[i] != o.n[i]) {
                return n[i] < o.n[i];
            }
        }
        return false;
    }
};

template <typename idx_t>
static HexDomFaceKey hexdom_face_key(const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                     const LocalSideTable                                   &lst,
                                     const ptrdiff_t                                         e,
                                     const int                                               f,
                                     const int                                               nnxs) {
    HexDomFaceKey key;
    key.n[0] = key.n[1] = key.n[2] = key.n[3] = invalid_idx<idx_t>();
    idx_t tmp[4];
    for (int i = 0; i < nnxs; ++i) {
        tmp[i] = m_elements[lst(f, i)][e];
    }
    for (int i = 0; i < nnxs; ++i) {
        idx_t m = tmp[i];
        int   p = i;
        for (int j = i + 1; j < nnxs; ++j) {
            if (tmp[j] < m) {
                m = tmp[j];
                p = j;
            }
        }
        const idx_t swap = tmp[i];
        tmp[i]           = tmp[p];
        tmp[p]           = swap;
        key.n[i]         = tmp[i];
    }
    return key;
}

template <typename idx_t>
int ssmixed_hex_dominant_generate_elements_blocks(
        const int                                                                     L,
        const ptrdiff_t                                                               n_blocks,
        const enum ElemType *const                                                    block_types,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               m_nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        const element_idx_t *const *const                                             /*hft*/,
        const block_idx_t *const *const                                               /*hnbb*/,
        ptrdiff_t                                                                    *n_unique_nodes_out,
        ptrdiff_t                                                                    *interior_start_out) {
    SMESH_ASSERT(L >= 1);
    double tick = time_seconds();

    int hex_corners[8] = {sshex8_lidx(L, 0, 0, 0),
                          sshex8_lidx(L, L, 0, 0),
                          sshex8_lidx(L, L, L, 0),
                          sshex8_lidx(L, 0, L, 0),
                          sshex8_lidx(L, 0, 0, L),
                          sshex8_lidx(L, L, 0, L),
                          sshex8_lidx(L, L, L, L),
                          sshex8_lidx(L, 0, L, L)};
    int tet_corners[4] = {sstet4_lidx(L, 0, 0, 0),
                          sstet4_lidx(L, L, 0, 0),
                          sstet4_lidx(L, 0, L, 0),
                          sstet4_lidx(L, 0, 0, L)};
    int wedge_corners[6];
    int pyr_corners[5];
    for (int d = 0; d < 6; ++d) {
        wedge_corners[d] = sswedge_lidx(L, hexdom_wedge_xyz[d][0] * L, hexdom_wedge_xyz[d][1] * L, hexdom_wedge_xyz[d][2] * L);
    }
    for (int d = 0; d < 5; ++d) {
        pyr_corners[d] = hexdom_pyr_corner_lidx(L, d);
    }

    int *hex_coords[3], *tet_coords[3], *wedge_coords[3], *pyr_coords[3];
    for (int d = 0; d < 3; d++) {
        hex_coords[d]   = (int *)SMESH_ALLOC(static_cast<size_t>(sshex8_nxe(L)) * sizeof(int));
        tet_coords[d]   = (int *)SMESH_ALLOC(static_cast<size_t>(sstet4_nxe(L)) * sizeof(int));
        wedge_coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(sswedge_nxe(L)) * sizeof(int));
        pyr_coords[d]   = (int *)SMESH_ALLOC(static_cast<size_t>(sspyramid_nxe(L)) * sizeof(int));
    }
    for (int zi = 0; zi <= L; zi++) {
        for (int yi = 0; yi <= L; yi++) {
            for (int xi = 0; xi <= L; xi++) {
                const int lidx      = sshex8_lidx(L, xi, yi, zi);
                hex_coords[0][lidx] = xi;
                hex_coords[1][lidx] = yi;
                hex_coords[2][lidx] = zi;
            }
        }
    }
    for (int z = 0; z <= L; ++z) {
        for (int y = 0; y <= L - z; ++y) {
            for (int x = 0; x <= L - z - y; ++x) {
                const int lidx      = sstet4_lidx(L, x, y, z);
                tet_coords[0][lidx] = x;
                tet_coords[1][lidx] = y;
                tet_coords[2][lidx] = z;
            }
        }
    }
    for (int z = 0; z <= L; ++z) {
        for (int y = 0; y <= L; ++y) {
            for (int x = 0; x <= L - y; ++x) {
                const int lidx        = sswedge_lidx(L, x, y, z);
                wedge_coords[0][lidx] = x;
                wedge_coords[1][lidx] = y;
                wedge_coords[2][lidx] = z;
            }
        }
    }
    for (int k = 0; k <= L; ++k) {
        for (int j = 0; j <= L - k; ++j) {
            for (int i = 0; i <= L - k; ++i) {
                const int lidx      = sspyramid_lidx(L, i, j, k);
                pyr_coords[0][lidx] = i;
                pyr_coords[1][lidx] = j;
                pyr_coords[2][lidx] = k;
            }
        }
    }

    ptrdiff_t n_e_total = 0, n_e_hex = 0, n_e_tet = 0, n_e_wedge = 0, n_e_pyr = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        n_e_total += n_elements[b];
        const enum ElemType f = block_types[b];
        if (f == HEX8) {
            n_e_hex += n_elements[b];
            for (int d = 0; d < 8; d++) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    elements[b][hex_corners[d]][e] = m_elements[b][d][e];
                }
            }
        } else if (f == TET4) {
            n_e_tet += n_elements[b];
            for (int d = 0; d < 4; d++) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    elements[b][tet_corners[d]][e] = m_elements[b][d][e];
                }
            }
        } else if (f == WEDGE6) {
            n_e_wedge += n_elements[b];
            for (int d = 0; d < 6; d++) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    elements[b][wedge_corners[d]][e] = m_elements[b][d][e];
                }
            }
        } else if (f == PYRAMID5) {
            n_e_pyr += n_elements[b];
            for (int d = 0; d < 5; d++) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    elements[b][pyr_corners[d]][e] = m_elements[b][d][e];
                }
            }
        } else {
            SMESH_ERROR("ssmixed hex-dominant: unsupported block type %s\n", type_to_string(f));
            return SMESH_FAILURE;
        }
    }

    idx_t     index_base = static_cast<idx_t>(m_nnodes);
    const int nxedge     = (L > 1) ? (L - 1) : 0;
    if (nxedge) {
        idx_t       *n2eptr{nullptr};
        idx_t       *elindex{nullptr};
        block_idx_t *block_number{nullptr};
        create_multiblock_n2e<idx_t, idx_t, idx_t>(static_cast<block_idx_t>(n_blocks),
                                                   block_types,
                                                   n_elements,
                                                   m_elements,
                                                   m_nnodes,
                                                   &block_number,
                                                   &n2eptr,
                                                   &elindex);
        idx_t *rowptr{nullptr};
        idx_t *colidx{nullptr};
        mixed_hex_dominant_build_edge_graph_from_n2e<idx_t, idx_t, idx_t>(n_blocks,
                                                                          block_types,
                                                                          n_elements,
                                                                          m_nnodes,
                                                                          m_elements,
                                                                          n2eptr,
                                                                          elindex,
                                                                          block_number,
                                                                          &rowptr,
                                                                          &colidx);
        const ptrdiff_t nedges   = rowptr[m_nnodes] / 2;
        const ptrdiff_t nnz      = rowptr[m_nnodes];
        idx_t          *edge_idx = (idx_t *)SMESH_CALLOC(static_cast<size_t>(nnz), sizeof(idx_t));
        ptrdiff_t       edge_count = 0;
        idx_t           next_id    = 0;
        for (ptrdiff_t i = 0; i < m_nnodes; i++) {
            for (idx_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                if (i < colidx[k]) {
                    edge_count += 1;
                    edge_idx[k] = next_id++;
                }
            }
        }
        SMESH_ASSERT(edge_count == nedges);
        SMESH_UNUSED(edge_count);
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            if (block_types[b] == HEX8) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    mixed_write_hex_owned_edges(L, e, m_elements[b], rowptr, colidx, edge_idx, hex_corners, hex_coords, index_base, nxedge, elements[b]);
                }
            } else if (block_types[b] == TET4) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    mixed_write_tet_owned_edges(L, e, m_elements[b], rowptr, colidx, edge_idx, tet_corners, tet_coords, index_base, nxedge, elements[b]);
                }
            } else if (block_types[b] == WEDGE6) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    hexdom_write_wedge_edges(L, e, m_elements[b], rowptr, colidx, edge_idx, wedge_corners, wedge_coords, index_base, nxedge, elements[b]);
                }
            } else {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    hexdom_write_pyr_edges(L, e, m_elements[b], rowptr, colidx, edge_idx, pyr_corners, pyr_coords, index_base, nxedge, elements[b]);
                }
            }
        }
        SMESH_FREE(rowptr);
        SMESH_FREE(colidx);
        SMESH_FREE(edge_idx);
        SMESH_FREE(n2eptr);
        SMESH_FREE(elindex);
        SMESH_FREE(block_number);
        index_base += static_cast<idx_t>(nedges * nxedge);
    }

    LocalSideTable lst_hex, lst_tet, lst_wedge, lst_pyr;
    lst_hex.fill(HEX8);
    lst_tet.fill(TET4);
    lst_wedge.fill(WEDGE6);
    lst_pyr.fill(PYRAMID5);

    const int nxf_tri  = sstet4_nxface(L);
    const int nxf_quad = (L > 1) ? ((L - 1) * (L - 1)) : 0;

    auto family_lst    = [&](const enum ElemType f) -> LocalSideTable & {
        if (f == HEX8) {
            return lst_hex;
        }
        if (f == TET4) {
            return lst_tet;
        }
        if (f == WEDGE6) {
            return lst_wedge;
        }
        return lst_pyr;
    };
    auto family_lagr = [&](const enum ElemType f) -> int * {
        if (f == HEX8) {
            return hex_corners;
        }
        if (f == TET4) {
            return tet_corners;
        }
        if (f == WEDGE6) {
            return wedge_corners;
        }
        return pyr_corners;
    };
    auto family_coords = [&](const enum ElemType f) -> int ** {
        if (f == HEX8) {
            return hex_coords;
        }
        if (f == TET4) {
            return tet_coords;
        }
        if (f == WEDGE6) {
            return wedge_coords;
        }
        return pyr_coords;
    };

    if (nxf_tri) {
        std::map<HexDomFaceKey, idx_t> tri_off;
        idx_t                          n_tri = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            const enum ElemType f      = block_types[b];
            const int           nsides = hexdom_nsides(f);
            auto               &lst    = family_lst(f);
            int                *lagr   = family_lagr(f);
            int               **coords = family_coords(f);
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                for (int s = 0; s < nsides; ++s) {
                    if (hexdom_side_nnxs(f, s) != 3) {
                        continue;
                    }
                    const HexDomFaceKey key = hexdom_face_key(m_elements[b], lst, e, s, 3);
                    auto                it  = tri_off.find(key);
                    idx_t               off;
                    if (it == tri_off.end()) {
                        off            = index_base + n_tri * nxf_tri;
                        tri_off[key]   = off;
                        n_tri++;
                    } else {
                        off = it->second;
                    }
                    hexdom_index_tri_face(f, L, m_elements[b], lst.table, lagr, coords, off, e, s, elements[b]);
                }
            }
        }
        index_base += n_tri * nxf_tri;
    }

    if (nxf_quad) {
        std::map<HexDomFaceKey, idx_t> quad_off;
        idx_t                          n_quad = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            const enum ElemType f      = block_types[b];
            const int           nsides = hexdom_nsides(f);
            auto               &lst    = family_lst(f);
            int                *lagr   = family_lagr(f);
            int               **coords = family_coords(f);
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                for (int s = 0; s < nsides; ++s) {
                    if (hexdom_side_nnxs(f, s) != 4) {
                        continue;
                    }
                    const HexDomFaceKey key = hexdom_face_key(m_elements[b], lst, e, s, 4);
                    auto                it  = quad_off.find(key);
                    idx_t               off;
                    if (it == quad_off.end()) {
                        off             = index_base + n_quad * nxf_quad;
                        quad_off[key]   = off;
                        n_quad++;
                    } else {
                        off = it->second;
                    }
                    hexdom_index_quad_face(f, L, m_elements[b], lst.table, lagr, coords, off, e, s, elements[b]);
                }
            }
        }
        index_base += n_quad * nxf_quad;
    }

    const int       nxvol_hex   = (L - 1) * (L - 1) * (L - 1);
    const int       nxvol_tet   = sstet4_nxvol(L);
    const int       nxvol_wedge = sswedge_nxvol(L);
    const int       nxvol_pyr   = sspyramid_nxvol(L);
    const ptrdiff_t interior_start = index_base;

    if (nxvol_hex) {
        ptrdiff_t e_off = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            if (block_types[b] != HEX8) {
                continue;
            }
            for (int zi = 1; zi < L; zi++) {
                for (int yi = 1; yi < L; yi++) {
                    for (int xi = 1; xi < L; xi++) {
                        const int lidx_vol = sshex8_lidx(L, xi, yi, zi);
                        const int Lm1      = L - 1;
                        const int en       = (zi - 1) * Lm1 * Lm1 + (yi - 1) * Lm1 + xi - 1;
                        for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                            elements[b][lidx_vol][e] = index_base + (e_off + e) * nxvol_hex + en;
                        }
                    }
                }
            }
            e_off += n_elements[b];
        }
        index_base += static_cast<idx_t>(n_e_hex * nxvol_hex);
    }
    if (nxvol_tet) {
        ptrdiff_t e_off = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            if (block_types[b] != TET4) {
                continue;
            }
            for (int z = 1; z <= L - 3; ++z) {
                for (int y = 1; y <= L - 2 - z; ++y) {
                    for (int x = 1; x <= L - 1 - z - y; ++x) {
                        const int lidx_vol = sstet4_lidx(L, x, y, z);
                        const int en       = sstet4_lidx(L - 4, x - 1, y - 1, z - 1);
                        for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                            elements[b][lidx_vol][e] = index_base + (e_off + e) * nxvol_tet + en;
                        }
                    }
                }
            }
            e_off += n_elements[b];
        }
        index_base += static_cast<idx_t>(n_e_tet * nxvol_tet);
    }
    if (nxvol_wedge) {
        ptrdiff_t e_off = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            if (block_types[b] != WEDGE6) {
                continue;
            }
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                int en = 0;
                for (int z = 1; z <= L - 1; ++z) {
                    for (int y = 1; y <= L - 2; ++y) {
                        for (int x = 1; x <= L - 1 - y; ++x) {
                            elements[b][sswedge_lidx(L, x, y, z)][e] =
                                    index_base + (e_off + e) * nxvol_wedge + en++;
                        }
                    }
                }
            }
            e_off += n_elements[b];
        }
        index_base += static_cast<idx_t>(n_e_wedge * nxvol_wedge);
    }
    if (nxvol_pyr) {
        ptrdiff_t e_off = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            if (block_types[b] != PYRAMID5) {
                continue;
            }
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                int en = 0;
                for (int k = 1; k <= L - 2; ++k) {
                    for (int j = 1; j <= L - k - 1; ++j) {
                        for (int i = 1; i <= L - k - 1; ++i) {
                            elements[b][sspyramid_lidx(L, i, j, k)][e] =
                                    index_base + (e_off + e) * nxvol_pyr + en++;
                        }
                    }
                }
            }
            e_off += n_elements[b];
        }
        index_base += static_cast<idx_t>(n_e_pyr * nxvol_pyr);
    }

    for (int d = 0; d < 3; d++) {
        SMESH_FREE(hex_coords[d]);
        SMESH_FREE(tet_coords[d]);
        SMESH_FREE(wedge_coords[d]);
        SMESH_FREE(pyr_coords[d]);
    }
    *n_unique_nodes_out = index_base;
    *interior_start_out = interior_start;
    printf("Create idx (HEX-dominant mixed) took\t%g [s]\n", time_seconds() - tick);
    printf("#macroelements %ld, #macronodes %ld\n", n_e_total, m_nnodes);
    printf("#micronodes %ld\n", *n_unique_nodes_out);
    return SMESH_SUCCESS;
}

template <typename idx_t>
int ssmixed_hex_dominant_hierarchical_renumbering_blocks(
        const int                                                                     L,
        const int                                                                     nlevels,
        int *const                                                                    levels,
        const ptrdiff_t                                                               n_blocks,
        const enum ElemType *const                                                    block_types,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               nnodes,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        idx_t *const SMESH_RESTRICT                                                   node_mapping,
        const bool                                                                    preserve_corner_ordering) {
    static const int tet_corner_xyz[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < nnodes; i++) {
        node_mapping[i] = invalid_idx<idx_t>();
    }
    idx_t next_id = 0;

    auto visit_idx = [&](const idx_t idx, const bool assign_sequential) {
        SMESH_ASSERT(idx < nnodes);
        if (assign_sequential) {
            if (node_mapping[idx] == invalid_idx<idx_t>()) {
                node_mapping[idx] = next_id++;
            }
        } else {
            node_mapping[idx] = idx;
            next_id           = std::max(next_id, idx);
        }
    };

    auto visit_corners = [&](const bool assign_sequential) {
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            const enum ElemType f = block_types[b];
            if (f == HEX8) {
                for (int zi = 0; zi <= 1; zi++) {
                    for (int yi = 0; yi <= 1; yi++) {
                        for (int xi = 0; xi <= 1; xi++) {
                            const int v = sshex8_lidx(L, xi * L, yi * L, zi * L);
                            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                                visit_idx(elements[b][v][e], assign_sequential);
                            }
                        }
                    }
                }
            } else if (f == TET4) {
                for (int c = 0; c < 4; c++) {
                    const int v = sstet4_lidx(L, tet_corner_xyz[c][0] * L, tet_corner_xyz[c][1] * L, tet_corner_xyz[c][2] * L);
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        visit_idx(elements[b][v][e], assign_sequential);
                    }
                }
            } else if (f == WEDGE6) {
                for (int c = 0; c < 6; c++) {
                    const int v = sswedge_lidx(L, hexdom_wedge_xyz[c][0] * L, hexdom_wedge_xyz[c][1] * L, hexdom_wedge_xyz[c][2] * L);
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        visit_idx(elements[b][v][e], assign_sequential);
                    }
                }
            } else {
                for (int c = 0; c < 5; c++) {
                    const int v = hexdom_pyr_corner_lidx(L, c);
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        visit_idx(elements[b][v][e], assign_sequential);
                    }
                }
            }
        }
    };

    if (preserve_corner_ordering) {
        visit_corners(false);
        next_id++;
    } else {
        visit_corners(true);
    }

    for (int k = 1; k < nlevels; k++) {
        const int l           = levels[k];
        const int step_factor = L / l;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            const enum ElemType f = block_types[b];
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                if (f == HEX8) {
                    for (int zi = 0; zi <= l; zi++) {
                        for (int yi = 0; yi <= l; yi++) {
                            for (int xi = 0; xi <= l; xi++) {
                                visit_idx(elements[b][sshex8_lidx(L, xi * step_factor, yi * step_factor, zi * step_factor)][e], true);
                            }
                        }
                    }
                } else if (f == TET4) {
                    for (int z = 0; z <= l; ++z) {
                        for (int y = 0; y <= l - z; ++y) {
                            for (int x = 0; x <= l - z - y; ++x) {
                                visit_idx(elements[b][sstet4_lidx(L, x * step_factor, y * step_factor, z * step_factor)][e], true);
                            }
                        }
                    }
                } else if (f == WEDGE6) {
                    for (int z = 0; z <= l; ++z) {
                        for (int y = 0; y <= l; ++y) {
                            for (int x = 0; x <= l - y; ++x) {
                                visit_idx(elements[b][sswedge_lidx(L, x * step_factor, y * step_factor, z * step_factor)][e], true);
                            }
                        }
                    }
                } else {
                    for (int kk = 0; kk <= l; ++kk) {
                        for (int j = 0; j <= l - kk; ++j) {
                            for (int i = 0; i <= l - kk; ++i) {
                                visit_idx(elements[b][sspyramid_lidx(L, i * step_factor, j * step_factor, kk * step_factor)][e], true);
                            }
                        }
                    }
                }
            }
        }
    }

    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const enum ElemType f = block_types[b];
        const int           nxe =
                (f == HEX8)      ? sshex8_nxe(L)
                : (f == TET4)    ? sstet4_nxe(L)
                : (f == WEDGE6)  ? sswedge_nxe(L)
                                 : sspyramid_nxe(L);
        for (int v = 0; v < nxe; ++v) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                const idx_t idx = elements[b][v][e];
                if (node_mapping[idx] == invalid_idx<idx_t>()) {
                    SMESH_ERROR("Uninitialized node mapping\n");
                }
                elements[b][v][e] = node_mapping[idx];
            }
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
