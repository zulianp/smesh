#ifndef SMESH_SSPYRAMID_GRAPH_IMPL_HPP
#define SMESH_SSPYRAMID_GRAPH_IMPL_HPP

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_buffer.hpp"
#include "smesh_graph.hpp"
#include "smesh_multiblock_graph.hpp"
#include "smesh_search.hpp"
#include "smesh_sort.hpp"
#include "smesh_sspyramid.hpp"
#include "smesh_sspyramid_graph.hpp"

#include <algorithm>

namespace smesh {

static const int pyramid5_n_edge_neigh[5]    = {3, 3, 3, 3, 4};
static int       pyramid5_edge_connectivity[5][4] = {
        {1, 3, 4, -1}, {0, 2, 4, -1}, {1, 3, 4, -1}, {0, 2, 4, -1}, {0, 1, 2, 3}};

static const int sspyramid_corner_ijk[5][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}};

template <typename idx_t, typename count_t, typename element_idx_t>
static int pyramid5_build_edge_graph_from_n2e(const ptrdiff_t                                         nelements,
                                              const ptrdiff_t                                         nnodes,
                                              const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                              const count_t *const SMESH_RESTRICT                     n2eptr,
                                              const element_idx_t *const SMESH_RESTRICT               elindex,
                                              count_t                                               **out_rowptr,
                                              idx_t                                                 **out_colidx) {
    SMESH_UNUSED(nelements);
    count_t *rowptr = (count_t *)SMESH_ALLOC((nnodes + 1) * sizeof(count_t));
    idx_t   *colidx = 0;
    rowptr[0]       = 0;
#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        idx_t   n2nbuff[2048];
        count_t ebegin  = n2eptr[node];
        count_t eend    = n2eptr[node + 1];
        idx_t   nneighs = 0;
        for (count_t e = ebegin; e < eend; ++e) {
            element_idx_t eidx = elindex[e];
            int           lidx = -1;
            for (int edof_i = 0; edof_i < 5; ++edof_i) {
                if (elems[edof_i][eidx] == node) {
                    lidx = edof_i;
                    break;
                }
            }
            SMESH_ASSERT(lidx != -1);
            const int nn = pyramid5_n_edge_neigh[lidx];
            for (int d = 0; d < nn; d++) {
                SMESH_ASSERT(nneighs < 2048);
                n2nbuff[nneighs++] = elems[pyramid5_edge_connectivity[lidx][d]][eidx];
            }
        }
        nneighs          = sort_and_unique(n2nbuff, nneighs);
        rowptr[node + 1] = nneighs;
    }
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        rowptr[node + 1] += rowptr[node];
    }
    colidx = (idx_t *)SMESH_ALLOC(rowptr[nnodes] * sizeof(idx_t));
#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        idx_t   n2nbuff[2048];
        count_t ebegin  = n2eptr[node];
        count_t eend    = n2eptr[node + 1];
        idx_t   nneighs = 0;
        for (count_t e = ebegin; e < eend; ++e) {
            element_idx_t eidx = elindex[e];
            int           lidx = -1;
            for (int edof_i = 0; edof_i < 5; ++edof_i) {
                if (elems[edof_i][eidx] == node) {
                    lidx = edof_i;
                    break;
                }
            }
            const int nn = pyramid5_n_edge_neigh[lidx];
            for (int d = 0; d < nn; d++) {
                n2nbuff[nneighs++] = elems[pyramid5_edge_connectivity[lidx][d]][eidx];
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

template <typename idx_t, typename count_t, typename element_idx_t>
static int pyramid5_build_edge_graph(const ptrdiff_t                                         nelements,
                                     const ptrdiff_t                                         nnodes,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                     count_t                                               **out_rowptr,
                                     idx_t                                                 **out_colidx) {
    count_t       *n2eptr;
    element_idx_t *elindex;
    create_n2e<idx_t, count_t, element_idx_t>(nelements, nnodes, 5, elems, &n2eptr, &elindex);
    int err = pyramid5_build_edge_graph_from_n2e(nelements, nnodes, elems, n2eptr, elindex, out_rowptr, out_colidx);
    SMESH_FREE(n2eptr);
    SMESH_FREE(elindex);
    return err;
}

template <typename idx_t, typename count_t, typename element_idx_t>
static int pyramid5_build_edge_graph_from_multiblock_n2e(
        const ptrdiff_t                                                                n_blocks,
        const ptrdiff_t *const                                                         n_elements,
        const ptrdiff_t                                                                nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
        const count_t *const SMESH_RESTRICT                                            n2eptr,
        const element_idx_t *const SMESH_RESTRICT                                      elindex,
        const block_idx_t *const SMESH_RESTRICT                                        block_number,
        count_t                                                                      **out_rowptr,
        idx_t                                                                        **out_colidx) {
    count_t *rowptr = (count_t *)SMESH_ALLOC((nnodes + 1) * sizeof(count_t));
    idx_t   *colidx = 0;
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
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];
            int lidx = -1;
            for (int edof_i = 0; edof_i < 5; ++edof_i) {
                if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                    lidx = edof_i;
                    break;
                }
            }
            SMESH_ASSERT(lidx != -1);
            const int nn = pyramid5_n_edge_neigh[lidx];
            for (int d = 0; d < nn; d++) {
                SMESH_ASSERT(nneighs < 2048);
                n2nbuff[nneighs++] = block_elems[pyramid5_edge_connectivity[lidx][d]][eidx];
            }
        }
        nneighs          = sort_and_unique(n2nbuff, nneighs);
        rowptr[node + 1] = nneighs;
    }
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        rowptr[node + 1] += rowptr[node];
    }
    colidx = (idx_t *)SMESH_ALLOC(rowptr[nnodes] * sizeof(idx_t));
#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
        idx_t   n2nbuff[2048];
        count_t ebegin  = n2eptr[node];
        count_t eend    = n2eptr[node + 1];
        idx_t   nneighs = 0;
        for (count_t e = ebegin; e < eend; ++e) {
            const block_idx_t   b    = block_number[e];
            const element_idx_t eidx = elindex[e];
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];
            int lidx = -1;
            for (int edof_i = 0; edof_i < 5; ++edof_i) {
                if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                    lidx = edof_i;
                    break;
                }
            }
            const int nn = pyramid5_n_edge_neigh[lidx];
            for (int d = 0; d < nn; d++) {
                n2nbuff[nneighs++] = block_elems[pyramid5_edge_connectivity[lidx][d]][eidx];
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

static void sspyramid_fill_lattice_coords(const int L, int **coords) {
    const int nxe = sspyramid_nxe(L);
    for (int d = 0; d < 3; d++) {
        coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(int));
    }
    for (int k = 0; k <= L; ++k) {
        for (int j = 0; j <= L - k; ++j) {
            for (int i = 0; i <= L - k; ++i) {
                const int lidx  = sspyramid_lidx(L, i, j, k);
                coords[0][lidx] = i;
                coords[1][lidx] = j;
                coords[2][lidx] = k;
            }
        }
    }
}

static int sspyramid_corner_lidx(const int L, const int c) {
    const int i = sspyramid_corner_ijk[c][0] * L;
    const int j = sspyramid_corner_ijk[c][1] * L;
    const int k = sspyramid_corner_ijk[c][2] * L;
    if (c == 4) {
        return sspyramid_lidx(L, 0, 0, L);
    }
    return sspyramid_lidx(L, i, j, k);
}

template <typename idx_t>
static void sspyramid_index_tri_face(const int                                               L,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                     const int *const                                        local_side_table,
                                     const int *const                                        lagr_to_proteus_corners,
                                     int *const *const                                       coords,
                                     const idx_t                                             global_face_offset,
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
    const int lidx_o = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_o]];
    const int lidx_u = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_u]];
    const int lidx_v = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_v]];
    int       Po[3], Pu[3], Pv[3];
    for (int d = 0; d < 3; d++) {
        Po[d] = coords[d][lidx_o];
        Pu[d] = coords[d][lidx_u];
        Pv[d] = coords[d][lidx_v];
    }
    int local_offset = 0;
    for (int t = 1; t <= L - 2; ++t) {
        for (int s = 1; s <= L - 1 - t; ++s) {
            const int w  = L - s - t;
            const int xi = (Po[0] * w + Pu[0] * s + Pv[0] * t) / L;
            const int yi = (Po[1] * w + Pu[1] * s + Pv[1] * t) / L;
            const int zi = (Po[2] * w + Pu[2] * s + Pv[2] * t) / L;
            const int pidx = sspyramid_lidx(L, xi, yi, zi);
            elements[pidx][e] = global_face_offset + local_offset++;
        }
    }
}

template <typename idx_t>
static void sspyramid_index_quad_face(const int                                               L,
                                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                      const int *const                                        local_side_table,
                                      const int *const                                        lagr_to_proteus_corners,
                                      int *const *const                                       coords,
                                      const idx_t                                             global_face_offset,
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
    const int lidx_o = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_o]];
    const int lidx_u = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_u]];
    const int lidx_v = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_v]];
    int       Po[3], Pu[3], Pv[3];
    for (int d = 0; d < 3; d++) {
        Po[d] = coords[d][lidx_o];
        Pu[d] = coords[d][lidx_u];
        Pv[d] = coords[d][lidx_v];
    }
    int local_offset = 0;
    for (int t = 1; t < L; ++t) {
        for (int s = 1; s < L; ++s) {
            const int xi = Po[0] + ((Pu[0] - Po[0]) * s + (Pv[0] - Po[0]) * t) / L;
            const int yi = Po[1] + ((Pu[1] - Po[1]) * s + (Pv[1] - Po[1]) * t) / L;
            const int zi = Po[2] + ((Pu[2] - Po[2]) * s + (Pv[2] - Po[2]) * t) / L;
            const int pidx = sspyramid_lidx(L, xi, yi, zi);
            elements[pidx][e] = global_face_offset + local_offset++;
        }
    }
}

template <typename idx_t>
static void sspyramid_index_owned_edges(const int                                               L,
                                        const ptrdiff_t                                         e,
                                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                        const idx_t *const                                      rowptr,
                                        const idx_t *const                                      colidx,
                                        const idx_t *const                                      edge_idx,
                                        const int *const                                        lagr_to_proteus_corners,
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
        const int          nn        = pyramid5_n_edge_neigh[d1];
        idx_t              g_neigh[4];
        for (int k = 0; k < nn; k++) {
            g_neigh[k] = nodes[pyramid5_edge_connectivity[d1][k]];
        }
        idx_t offsets[4];
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
            const int lid1 = lagr_to_proteus_corners[d1];
            const int lid2 = lagr_to_proteus_corners[pyramid5_edge_connectivity[d1][d2]];
            int       P1[3], P2[3];
            for (int d = 0; d < 3; d++) {
                P1[d] = coords[d][lid1];
                P2[d] = coords[d][lid2];
            }
            const idx_t edge_start = index_base + edge_view[offsets[d2]] * static_cast<idx_t>(nxedge);
            for (int t = 1; t < L; ++t) {
                const int xi   = (P1[0] * (L - t) + P2[0] * t) / L;
                const int yi   = (P1[1] * (L - t) + P2[1] * t) / L;
                const int zi   = (P1[2] * (L - t) + P2[2] * t) / L;
                const int lidx = sspyramid_lidx(L, xi, yi, zi);
                elements[lidx][e] = edge_start + (t - 1);
            }
        }
    }
}

template <typename idx_t>
static void sspyramid_index_volume(const int                                         L,
                                   const ptrdiff_t                                   e,
                                   const idx_t                                       index_base,
                                   const ptrdiff_t                                   e_global,
                                   const int                                         nxvol,
                                   idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements) {
    int en = 0;
    for (int k = 1; k <= L - 2; ++k) {
        for (int j = 1; j <= L - k - 1; ++j) {
            for (int i = 1; i <= L - k - 1; ++i) {
                const int lidx    = sspyramid_lidx(L, i, j, k);
                elements[lidx][e] = index_base + static_cast<idx_t>(e_global) * static_cast<idx_t>(nxvol) +
                                    static_cast<idx_t>(en);
                en++;
            }
        }
    }
}

template <typename idx_t>
int sspyramid_generate_elements(const int                                               L,
                                const ptrdiff_t                                         m_nelements,
                                const ptrdiff_t                                         m_nnodes,
                                const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
                                ptrdiff_t                                              *n_unique_nodes_out,
                                ptrdiff_t                                              *interior_start_out) {
    SMESH_ASSERT(L >= 1);
    int lagr_to_proteus_corners[5];
    for (int d = 0; d < 5; d++) {
        lagr_to_proteus_corners[d] = sspyramid_corner_lidx(L, d);
    }
    int *coords[3];
    sspyramid_fill_lattice_coords(L, coords);
    for (int d = 0; d < 5; d++) {
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            elements[lagr_to_proteus_corners[d]][e] = m_elements[d][e];
        }
    }
    idx_t     index_base = static_cast<idx_t>(m_nnodes);
    const int nxedge     = sspyramid_nxedge(L);
    if (nxedge) {
        idx_t *rowptr{nullptr};
        idx_t *colidx{nullptr};
        pyramid5_build_edge_graph<idx_t, idx_t, idx_t>(m_nelements, m_nnodes, m_elements, &rowptr, &colidx);
        const ptrdiff_t nedges     = rowptr[m_nnodes] / 2;
        const ptrdiff_t nnz        = rowptr[m_nnodes];
        idx_t          *edge_idx   = (idx_t *)SMESH_CALLOC(nnz, sizeof(idx_t));
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
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            sspyramid_index_owned_edges(L,
                                        e,
                                        m_elements,
                                        rowptr,
                                        colidx,
                                        edge_idx,
                                        lagr_to_proteus_corners,
                                        coords,
                                        index_base,
                                        nxedge,
                                        elements);
        }
        SMESH_FREE(rowptr);
        SMESH_FREE(colidx);
        SMESH_FREE(edge_idx);
        index_base += static_cast<idx_t>(nedges * nxedge);
    }

    LocalSideTable lst;
    lst.fill(PYRAMID5);
    element_idx_t *adj_table = 0;
    create_element_adj_table(m_nelements, m_nnodes, PYRAMID5, m_elements, &adj_table);
    const int nxf_tri  = sspyramid_nx_tri_face(L);
    const int nxf_quad = sspyramid_nx_quad_face(L);
    idx_t     n_tri    = 0;
    idx_t     n_quad   = 0;
    if (nxf_tri) {
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            for (int f = 0; f < 4; f++) {
                const element_idx_t neigh = adj_table[e * 5 + f];
                if (neigh != invalid_idx<element_idx_t>() && neigh < e) {
                    continue;
                }
                const idx_t off = index_base + n_tri * nxf_tri;
                sspyramid_index_tri_face(L, m_elements, lst.table, lagr_to_proteus_corners, coords, off, e, f, elements);
                if (neigh != invalid_idx<element_idx_t>()) {
                    int neigh_f = 5;
                    for (int f2 = 0; f2 < 5; f2++) {
                        if (e == adj_table[neigh * 5 + f2]) {
                            neigh_f = f2;
                            break;
                        }
                    }
                    SMESH_ASSERT(neigh_f != 5);
                    sspyramid_index_tri_face(
                            L, m_elements, lst.table, lagr_to_proteus_corners, coords, off, neigh, neigh_f, elements);
                }
                n_tri++;
            }
        }
    }
    index_base += n_tri * nxf_tri;
    if (nxf_quad) {
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            const int               f     = 4;
            const element_idx_t     neigh = adj_table[e * 5 + f];
            if (neigh != invalid_idx<element_idx_t>() && neigh < e) {
                continue;
            }
            const idx_t off = index_base + n_quad * nxf_quad;
            sspyramid_index_quad_face(L, m_elements, lst.table, lagr_to_proteus_corners, coords, off, e, f, elements);
            if (neigh != invalid_idx<element_idx_t>()) {
                int neigh_f = 5;
                for (int f2 = 0; f2 < 5; f2++) {
                    if (e == adj_table[neigh * 5 + f2]) {
                        neigh_f = f2;
                        break;
                    }
                }
                SMESH_ASSERT(neigh_f != 5);
                sspyramid_index_quad_face(
                        L, m_elements, lst.table, lagr_to_proteus_corners, coords, off, neigh, neigh_f, elements);
            }
            n_quad++;
        }
    }
    SMESH_FREE(adj_table);
    index_base += n_quad * nxf_quad;

    const int       nxvol          = sspyramid_nxvol(L);
    const ptrdiff_t interior_start = index_base;
    if (nxvol) {
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            sspyramid_index_volume(L, e, index_base, e, nxvol, elements);
        }
    }
    for (int d = 0; d < 3; d++) {
        SMESH_FREE(coords[d]);
    }
    *n_unique_nodes_out = interior_start + m_nelements * nxvol;
    *interior_start_out = interior_start;
    return SMESH_SUCCESS;
}

template <typename idx_t>
int sspyramid_generate_elements_blocks(
        const int                                                                     L,
        const ptrdiff_t                                                               n_blocks,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               m_nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        const element_idx_t *const *const                                             hft,
        const block_idx_t *const *const                                               hnbb,
        ptrdiff_t                                                                    *n_unique_nodes_out,
        ptrdiff_t                                                                    *interior_start_out) {
    SMESH_ASSERT(L >= 1);
    int lagr_to_proteus_corners[5];
    for (int d = 0; d < 5; d++) {
        lagr_to_proteus_corners[d] = sspyramid_corner_lidx(L, d);
    }
    int *coords[3];
    sspyramid_fill_lattice_coords(L, coords);
    ptrdiff_t n_e_total = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        n_e_total += n_elements[b];
        for (int d = 0; d < 5; d++) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                elements[b][lagr_to_proteus_corners[d]][e] = m_elements[b][d][e];
            }
        }
    }
    idx_t     index_base = static_cast<idx_t>(m_nnodes);
    const int nxedge     = sspyramid_nxedge(L);
    if (nxedge) {
        enum ElemType *types = (enum ElemType *)SMESH_ALLOC(static_cast<size_t>(n_blocks) * sizeof(enum ElemType));
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            types[b] = PYRAMID5;
        }
        idx_t       *n2eptr{nullptr};
        idx_t       *elindex{nullptr};
        block_idx_t *block_number{nullptr};
        create_multiblock_n2e<idx_t, idx_t, idx_t>(static_cast<block_idx_t>(n_blocks),
                                                   types,
                                                   n_elements,
                                                   m_elements,
                                                   m_nnodes,
                                                   &block_number,
                                                   &n2eptr,
                                                   &elindex);
        idx_t *rowptr{nullptr};
        idx_t *colidx{nullptr};
        pyramid5_build_edge_graph_from_multiblock_n2e<idx_t, idx_t, idx_t>(
                n_blocks, n_elements, m_nnodes, m_elements, n2eptr, elindex, block_number, &rowptr, &colidx);
        const ptrdiff_t nedges     = rowptr[m_nnodes] / 2;
        const ptrdiff_t nnz        = rowptr[m_nnodes];
        idx_t          *edge_idx   = (idx_t *)SMESH_CALLOC(nnz, sizeof(idx_t));
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
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                sspyramid_index_owned_edges(L,
                                            e,
                                            m_elements[b],
                                            rowptr,
                                            colidx,
                                            edge_idx,
                                            lagr_to_proteus_corners,
                                            coords,
                                            index_base,
                                            nxedge,
                                            elements[b]);
            }
        }
        SMESH_FREE(rowptr);
        SMESH_FREE(colidx);
        SMESH_FREE(edge_idx);
        SMESH_FREE(n2eptr);
        SMESH_FREE(elindex);
        SMESH_FREE(block_number);
        SMESH_FREE(types);
        index_base += static_cast<idx_t>(nedges * nxedge);
    }

    LocalSideTable lst;
    lst.fill(PYRAMID5);
    const int nxf_tri  = sspyramid_nx_tri_face(L);
    const int nxf_quad = sspyramid_nx_quad_face(L);
    idx_t     n_tri    = 0;
    idx_t     n_quad   = 0;
    if (nxf_tri) {
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                for (int f = 0; f < 4; f++) {
                    const element_idx_t neigh_e = hft[b][e * 5 + f];
                    const block_idx_t   neigh_b = hnbb[b][e * 5 + f];
                    if (neigh_e != invalid_idx<element_idx_t>() &&
                        (neigh_b < static_cast<block_idx_t>(b) ||
                         (neigh_b == static_cast<block_idx_t>(b) && neigh_e < static_cast<element_idx_t>(e)))) {
                        continue;
                    }
                    const idx_t off = index_base + n_tri * nxf_tri;
                    sspyramid_index_tri_face(
                            L, m_elements[b], lst.table, lagr_to_proteus_corners, coords, off, e, f, elements[b]);
                    if (neigh_e != invalid_idx<element_idx_t>()) {
                        int neigh_f = 5;
                        for (int f2 = 0; f2 < 5; f2++) {
                            if (hft[neigh_b][neigh_e * 5 + f2] == static_cast<element_idx_t>(e) &&
                                hnbb[neigh_b][neigh_e * 5 + f2] == static_cast<block_idx_t>(b)) {
                                neigh_f = f2;
                                break;
                            }
                        }
                        SMESH_ASSERT(neigh_f != 5);
                        sspyramid_index_tri_face(L,
                                                 m_elements[neigh_b],
                                                 lst.table,
                                                 lagr_to_proteus_corners,
                                                 coords,
                                                 off,
                                                 neigh_e,
                                                 neigh_f,
                                                 elements[neigh_b]);
                    }
                    n_tri++;
                }
            }
        }
    }
    index_base += n_tri * nxf_tri;
    if (nxf_quad) {
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                const int               f       = 4;
                const element_idx_t     neigh_e = hft[b][e * 5 + f];
                const block_idx_t       neigh_b = hnbb[b][e * 5 + f];
                if (neigh_e != invalid_idx<element_idx_t>() &&
                    (neigh_b < static_cast<block_idx_t>(b) ||
                     (neigh_b == static_cast<block_idx_t>(b) && neigh_e < static_cast<element_idx_t>(e)))) {
                    continue;
                }
                const idx_t off = index_base + n_quad * nxf_quad;
                sspyramid_index_quad_face(
                        L, m_elements[b], lst.table, lagr_to_proteus_corners, coords, off, e, f, elements[b]);
                if (neigh_e != invalid_idx<element_idx_t>()) {
                    int neigh_f = 5;
                    for (int f2 = 0; f2 < 5; f2++) {
                        if (hft[neigh_b][neigh_e * 5 + f2] == static_cast<element_idx_t>(e) &&
                            hnbb[neigh_b][neigh_e * 5 + f2] == static_cast<block_idx_t>(b)) {
                            neigh_f = f2;
                            break;
                        }
                    }
                    SMESH_ASSERT(neigh_f != 5);
                    sspyramid_index_quad_face(L,
                                              m_elements[neigh_b],
                                              lst.table,
                                              lagr_to_proteus_corners,
                                              coords,
                                              off,
                                              neigh_e,
                                              neigh_f,
                                              elements[neigh_b]);
                }
                n_quad++;
            }
        }
    }
    index_base += n_quad * nxf_quad;

    const int       nxvol          = sspyramid_nxvol(L);
    const ptrdiff_t interior_start = index_base;
    if (nxvol) {
        ptrdiff_t e_off = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                sspyramid_index_volume(L, e, index_base, e_off + e, nxvol, elements[b]);
            }
            e_off += n_elements[b];
        }
    }
    for (int d = 0; d < 3; d++) {
        SMESH_FREE(coords[d]);
    }
    *n_unique_nodes_out = interior_start + n_e_total * nxvol;
    *interior_start_out = interior_start;
    return SMESH_SUCCESS;
}

template <typename idx_t>
int sspyramid_hierarchical_renumbering(const int                                         L,
                                       const int                                         nlevels,
                                       int *const                                        levels,
                                       const ptrdiff_t                                   nelements,
                                       const ptrdiff_t                                   nnodes,
                                       idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                       idx_t *const SMESH_RESTRICT                       node_mapping,
                                       const bool                                        preserve_corner_ordering) {
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < nnodes; i++) {
        node_mapping[i] = invalid_idx<idx_t>();
    }
    idx_t next_id = 0;
    for (int c = 0; c < 5; c++) {
        const int v = sspyramid_corner_lidx(L, c);
        for (ptrdiff_t e = 0; e < nelements; e++) {
            const idx_t idx = elements[v][e];
            if (preserve_corner_ordering) {
                node_mapping[idx] = idx;
                next_id           = std::max(next_id, idx);
            } else if (node_mapping[idx] == invalid_idx<idx_t>()) {
                node_mapping[idx] = next_id++;
            }
        }
    }
    if (preserve_corner_ordering) {
        next_id++;
    }
    for (int klev = 1; klev < nlevels; klev++) {
        const int l           = levels[klev];
        const int step_factor = L / l;
        for (ptrdiff_t e = 0; e < nelements; e++) {
            for (int k = 0; k <= l; ++k) {
                for (int j = 0; j <= l - k; ++j) {
                    for (int i = 0; i <= l - k; ++i) {
                        const int   v   = sspyramid_lidx(L, i * step_factor, j * step_factor, k * step_factor);
                        const idx_t idx = elements[v][e];
                        if (node_mapping[idx] == invalid_idx<idx_t>()) {
                            node_mapping[idx] = next_id++;
                        }
                    }
                }
            }
        }
    }
    for (int k = 0; k <= L; ++k) {
        for (int j = 0; j <= L - k; ++j) {
            for (int i = 0; i <= L - k; ++i) {
                const int v = sspyramid_lidx(L, i, j, k);
                for (ptrdiff_t e = 0; e < nelements; e++) {
                    const idx_t idx = elements[v][e];
                    if (node_mapping[idx] == invalid_idx<idx_t>()) {
                        SMESH_ERROR("Uninitialized node mapping\n");
                    }
                    elements[v][e] = node_mapping[idx];
                }
            }
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t>
int sspyramid_hierarchical_renumbering_blocks(
        const int                                                                     L,
        const int                                                                     nlevels,
        int *const                                                                    levels,
        const ptrdiff_t                                                               n_blocks,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               nnodes,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        idx_t *const SMESH_RESTRICT                                                   node_mapping,
        const bool                                                                    preserve_corner_ordering) {
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < nnodes; i++) {
        node_mapping[i] = invalid_idx<idx_t>();
    }
    idx_t next_id = 0;
    for (int c = 0; c < 5; c++) {
        const int v = sspyramid_corner_lidx(L, c);
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                const idx_t idx = elements[b][v][e];
                if (preserve_corner_ordering) {
                    node_mapping[idx] = idx;
                    next_id           = std::max(next_id, idx);
                } else if (node_mapping[idx] == invalid_idx<idx_t>()) {
                    node_mapping[idx] = next_id++;
                }
            }
        }
    }
    if (preserve_corner_ordering) {
        next_id++;
    }
    for (int klev = 1; klev < nlevels; klev++) {
        const int l           = levels[klev];
        const int step_factor = L / l;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                for (int k = 0; k <= l; ++k) {
                    for (int j = 0; j <= l - k; ++j) {
                        for (int i = 0; i <= l - k; ++i) {
                            const int   v   = sspyramid_lidx(L, i * step_factor, j * step_factor, k * step_factor);
                            const idx_t idx = elements[b][v][e];
                            if (node_mapping[idx] == invalid_idx<idx_t>()) {
                                node_mapping[idx] = next_id++;
                            }
                        }
                    }
                }
            }
        }
    }
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        for (int k = 0; k <= L; ++k) {
            for (int j = 0; j <= L - k; ++j) {
                for (int i = 0; i <= L - k; ++i) {
                    const int v = sspyramid_lidx(L, i, j, k);
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        const idx_t idx = elements[b][v][e];
                        if (node_mapping[idx] == invalid_idx<idx_t>()) {
                            SMESH_ERROR("Uninitialized node mapping\n");
                        }
                        elements[b][v][e] = node_mapping[idx];
                    }
                }
            }
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t>
static void sspyramid_fill_tri_side(const int                                               L,
                                    const ptrdiff_t                                         e,
                                    const int                                               s,
                                    const LocalSideTable                                   &lst,
                                    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                    const ptrdiff_t                                         i,
                                    idx_t **const SMESH_RESTRICT                            sides) {
    const int c0    = lst(s, 0);
    const int c1    = lst(s, 1);
    const int c2    = lst(s, 2);
    const int P0[3] = {sspyramid_corner_ijk[c0][0] * L,
                       sspyramid_corner_ijk[c0][1] * L,
                       sspyramid_corner_ijk[c0][2] * L};
    const int P1[3] = {sspyramid_corner_ijk[c1][0] * L,
                       sspyramid_corner_ijk[c1][1] * L,
                       sspyramid_corner_ijk[c1][2] * L};
    const int P2[3] = {sspyramid_corner_ijk[c2][0] * L,
                       sspyramid_corner_ijk[c2][1] * L,
                       sspyramid_corner_ijk[c2][2] * L};
    int       lidx  = 0;
    for (int t = 0; t <= L; ++t) {
        for (int si = 0; si <= L - t; ++si) {
            const int w  = L - si - t;
            const int xi = (P0[0] * w + P1[0] * si + P2[0] * t) / L;
            const int yi = (P0[1] * w + P1[1] * si + P2[1] * t) / L;
            const int zi = (P0[2] * w + P1[2] * si + P2[2] * t) / L;
            sides[lidx++][i] = elems[sspyramid_lidx(L, xi, yi, zi)][e];
        }
    }
}

template <typename idx_t>
static void sspyramid_fill_quad_side(const int                                               L,
                                     const ptrdiff_t                                         e,
                                     const int                                               s,
                                     const LocalSideTable                                   &lst,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                     const ptrdiff_t                                         i,
                                     idx_t **const SMESH_RESTRICT                            sides) {
    const int c0    = lst(s, 0);
    const int c1    = lst(s, 1);
    const int c2    = lst(s, 2);
    const int c3    = lst(s, 3);
    const int P0[3] = {sspyramid_corner_ijk[c0][0] * L,
                       sspyramid_corner_ijk[c0][1] * L,
                       sspyramid_corner_ijk[c0][2] * L};
    const int P1[3] = {sspyramid_corner_ijk[c1][0] * L,
                       sspyramid_corner_ijk[c1][1] * L,
                       sspyramid_corner_ijk[c1][2] * L};
    const int P2[3] = {sspyramid_corner_ijk[c2][0] * L,
                       sspyramid_corner_ijk[c2][1] * L,
                       sspyramid_corner_ijk[c2][2] * L};
    const int P3[3] = {sspyramid_corner_ijk[c3][0] * L,
                       sspyramid_corner_ijk[c3][1] * L,
                       sspyramid_corner_ijk[c3][2] * L};
    int       lidx  = 0;
    for (int t = 0; t <= L; ++t) {
        for (int sxy = 0; sxy <= L; ++sxy) {
            const int xi =
                    (P0[0] * (L - sxy) * (L - t) + P1[0] * sxy * (L - t) + P2[0] * sxy * t + P3[0] * (L - sxy) * t) /
                    (L * L);
            const int yi =
                    (P0[1] * (L - sxy) * (L - t) + P1[1] * sxy * (L - t) + P2[1] * sxy * t + P3[1] * (L - sxy) * t) /
                    (L * L);
            const int zi =
                    (P0[2] * (L - sxy) * (L - t) + P1[2] * sxy * (L - t) + P2[2] * sxy * t + P3[2] * (L - sxy) * t) /
                    (L * L);
            sides[lidx++][i] = elems[sspyramid_lidx(L, xi, yi, zi)][e];
        }
    }
}

template <typename idx_t, typename element_idx_t>
int sspyramid_extract_surface_from_sideset(const int                                               L,
                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                           const ptrdiff_t                                         n_surf_elements,
                                           const element_idx_t *const SMESH_RESTRICT               parent_element,
                                           const i16 *const SMESH_RESTRICT                         side_idx,
                                           idx_t **const SMESH_RESTRICT                            sides) {
    LocalSideTable lst;
    lst.fill(PYRAMID5);
    if (n_surf_elements <= 0) {
        return SMESH_SUCCESS;
    }
    const int  first = side_idx[0];
    const bool quad  = first == 4;
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n_surf_elements; i++) {
        const ptrdiff_t e = parent_element[i];
        const int       s = side_idx[i];
        SMESH_ASSERT(s >= 0 && s < 5);
        SMESH_ASSERT((s == 4) == quad);
        if (s == 4) {
            sspyramid_fill_quad_side(L, e, s, lst, elems, i, sides);
        } else {
            sspyramid_fill_tri_side(L, e, s, lst, elems, i, sides);
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int sspyramid_extract_nodeset_from_sideset(const int                                               L,
                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                           const ptrdiff_t                                         n_surf_elements,
                                           const element_idx_t *const SMESH_RESTRICT               parent_element,
                                           const i16 *const SMESH_RESTRICT                         side_idx,
                                           ptrdiff_t                                              *n_nodes_out,
                                           idx_t **SMESH_RESTRICT                                  nodes_out) {
    if (n_surf_elements <= 0) {
        *n_nodes_out = 0;
        *nodes_out   = nullptr;
        return SMESH_SUCCESS;
    }
    const int       nnxs = (side_idx[0] == 4) ? ((L + 1) * (L + 1)) : sstet4_n_tri(L);
    const ptrdiff_t n    = static_cast<ptrdiff_t>(nnxs) * n_surf_elements;
    auto            tmp  = create_host_buffer<idx_t>(nnxs, n_surf_elements);
    if (sspyramid_extract_surface_from_sideset(L, elems, n_surf_elements, parent_element, side_idx, tmp->data()) !=
        SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    idx_t *nodes = static_cast<idx_t *>(SMESH_ALLOC(static_cast<size_t>(n) * sizeof(idx_t)));
    for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
        for (int v = 0; v < nnxs; ++v) {
            nodes[i * nnxs + v] = tmp->data()[v][i];
        }
    }
    *n_nodes_out = static_cast<ptrdiff_t>(sort_and_unique(nodes, static_cast<size_t>(n)));
    *nodes_out   = static_cast<idx_t *>(SMESH_REALLOC(nodes, static_cast<size_t>(*n_nodes_out) * sizeof(idx_t)));
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
