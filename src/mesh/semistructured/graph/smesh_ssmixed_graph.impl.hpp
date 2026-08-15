#ifndef SMESH_SSMIXED_GRAPH_IMPL_HPP
#define SMESH_SSMIXED_GRAPH_IMPL_HPP

#include "smesh_ssmixed_graph.hpp"

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_multiblock_graph.hpp"
#include "smesh_search.hpp"
#include "smesh_sort.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sstet4.hpp"

#include <algorithm>

namespace smesh {

    static int mixed_hex8_edge_connectivity[8][3] = {{1, 3, 4},
                                                     {0, 2, 5},
                                                     {1, 3, 6},
                                                     {0, 2, 7},
                                                     {0, 5, 7},
                                                     {1, 4, 6},
                                                     {2, 5, 7},
                                                     {3, 4, 6}};

    static int mixed_tet4_edge_connectivity[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

    template <typename idx_t, typename count_t, typename element_idx_t>
    static int mixed_hex_tet_build_edge_graph_from_n2e(
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
                SMESH_ASSERT(b >= 0 && b < n_blocks);
                SMESH_ASSERT(eidx < n_elements[b]);
                const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];
                const enum ElemType                                     et          = block_types[b];

                if (et == HEX8) {
                    int lidx = -1;
                    for (int edof_i = 0; edof_i < 8; ++edof_i) {
                        if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                            lidx = edof_i;
                            break;
                        }
                    }
                    SMESH_ASSERT(lidx != -1);
                    for (int d = 0; d < 3; d++) {
                        SMESH_ASSERT(nneighs < 2048);
                        n2nbuff[nneighs++] = block_elems[mixed_hex8_edge_connectivity[lidx][d]][eidx];
                    }
                } else if (et == TET4) {
                    int lidx = -1;
                    for (int edof_i = 0; edof_i < 4; ++edof_i) {
                        if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                            lidx = edof_i;
                            break;
                        }
                    }
                    SMESH_ASSERT(lidx != -1);
                    for (int d = 0; d < 3; d++) {
                        SMESH_ASSERT(nneighs < 2048);
                        n2nbuff[nneighs++] = block_elems[mixed_tet4_edge_connectivity[lidx][d]][eidx];
                    }
                } else {
                    SMESH_ERROR("ssmixed: unsupported block type %s\n", type_to_string(et));
                }
            }

            nneighs          = sort_and_unique(n2nbuff, nneighs);
            rowptr[node + 1] = nneighs;
        }

        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            rowptr[node + 1] += rowptr[node];
        }

        const ptrdiff_t nnz = rowptr[nnodes];
        colidx              = (idx_t *)SMESH_ALLOC(static_cast<size_t>(nnz) * sizeof(idx_t));

#pragma omp parallel for
        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            idx_t   n2nbuff[2048];
            count_t ebegin  = n2eptr[node];
            count_t eend    = n2eptr[node + 1];
            idx_t   nneighs = 0;

            for (count_t e = ebegin; e < eend; ++e) {
                const block_idx_t                                                   b           = block_number[e];
                const element_idx_t                                                 eidx        = elindex[e];
                const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT             block_elems = elems[b];
                const enum ElemType                                                 et          = block_types[b];

                if (et == HEX8) {
                    int lidx = 0;
                    for (int edof_i = 0; edof_i < 8; ++edof_i) {
                        if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                            lidx = edof_i;
                            break;
                        }
                    }
                    for (int d = 0; d < 3; d++) {
                        n2nbuff[nneighs++] = block_elems[mixed_hex8_edge_connectivity[lidx][d]][eidx];
                    }
                } else if (et == TET4) {
                    int lidx = 0;
                    for (int edof_i = 0; edof_i < 4; ++edof_i) {
                        if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                            lidx = edof_i;
                            break;
                        }
                    }
                    for (int d = 0; d < 3; d++) {
                        n2nbuff[nneighs++] = block_elems[mixed_tet4_edge_connectivity[lidx][d]][eidx];
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
    static void mixed_index_hex_face(const int                                               L,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                     const int *const                                        local_side_table,
                                     int                                                    *lagr_to_proteus_corners,
                                     int                                                   **coords,
                                     const idx_t                                             global_face_offset,
                                     const ptrdiff_t                                         e,
                                     const int                                               f,
                                     idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements) {
        int       argmin = 0;
        const int nnxs   = LocalSideTable::MAX_NUM_NODES_PER_SIDE;
        idx_t     valmin = m_elements[local_side_table[f * nnxs + 0]][e];
        for (int i = 0; i < 4; i++) {
            idx_t temp = m_elements[local_side_table[f * nnxs + i]][e];
            if (temp < valmin) {
                argmin = i;
                valmin = temp;
            }
        }

        int lst_o = argmin;
        int lst_u = ((lst_o + 1) % 4);
        int lst_v = ((lst_o + 3) % 4);
        if (m_elements[local_side_table[f * nnxs + lst_u]][e] > m_elements[local_side_table[f * nnxs + lst_v]][e]) {
            int temp = lst_v;
            lst_v    = lst_u;
            lst_u    = temp;
        }

        int lidx_o = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_o]];
        int lidx_u = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_u]];
        int lidx_v = lagr_to_proteus_corners[local_side_table[f * nnxs + lst_v]];

        int o_start[3];
        int u_len[3], u_dir[3];
        int v_len[3], v_dir[3];

        for (int d = 0; d < 3; d++) {
            o_start[d] = coords[d][lidx_o];
        }

        for (int d = 0; d < 3; d++) {
            int x    = coords[d][lidx_u] - coords[d][lidx_o];
            u_dir[d] = 1;
            u_len[d] = 1;
            if (x > 0) {
                x -= 1;
                u_len[d]   = x;
                o_start[d] = 1;
            } else if (x < 0) {
                x += 1;
                u_len[d]   = x;
                u_dir[d]   = -1;
                o_start[d] = L - 1;
            }
        }

        for (int d = 0; d < 3; d++) {
            int x    = coords[d][lidx_v] - coords[d][lidx_o];
            v_dir[d] = 1;
            v_len[d] = 1;
            if (x > 0) {
                x -= 1;
                v_len[d]   = x;
                o_start[d] = 1;
            } else if (x < 0) {
                x += 1;
                v_len[d]   = x;
                v_dir[d]   = -1;
                o_start[d] = L - 1;
            }
        }

        int local_offset = 0;
        for (int vzi = 0; vzi != v_len[2]; vzi += v_dir[2]) {
            for (int vyi = 0; vyi != v_len[1]; vyi += v_dir[1]) {
                for (int vxi = 0; vxi != v_len[0]; vxi += v_dir[0]) {
                    for (int uzi = 0; uzi != u_len[2]; uzi += u_dir[2]) {
                        for (int uyi = 0; uyi != u_len[1]; uyi += u_dir[1]) {
                            for (int uxi = 0; uxi != u_len[0]; uxi += u_dir[0]) {
                                int   pidx        = sshex8_lidx(L, uxi + vxi + o_start[0], uyi + vyi + o_start[1], uzi + vzi + o_start[2]);
                                idx_t fidx        = global_face_offset + local_offset++;
                                elements[pidx][e] = fidx;
                            }
                        }
                    }
                }
            }
        }
    }

    template <typename idx_t>
    static void mixed_index_tet_face(const int                                               L,
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

        int Po[3], Pu[3], Pv[3];
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
                const int pidx = sstet4_lidx(L, xi, yi, zi);
                elements[pidx][e] = global_face_offset + local_offset++;
            }
        }
    }

    template <typename idx_t>
    static void mixed_write_hex_owned_edges(const int                                               L,
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
        idx_t nodes[8];
        for (int d = 0; d < 8; d++) {
            nodes[d] = m_elements[d][e];
        }

        for (int d1 = 0; d1 < 8; d1++) {
            const idx_t        node1     = nodes[d1];
            const idx_t *const columns   = &colidx[rowptr[node1]];
            const idx_t *const edge_view = &edge_idx[rowptr[node1]];

            idx_t g_neigh[3];
            for (int k = 0; k < 3; k++) {
                g_neigh[k] = nodes[mixed_hex8_edge_connectivity[d1][k]];
            }

            idx_t offsets[3];
            find<3>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);
            idx_t g_edges[3];
            for (int d = 0; d < 3; d++) {
                g_edges[d] = edge_view[offsets[d]];
            }

            for (int d2 = 0; d2 < 3; d2++) {
                const idx_t node2 = g_neigh[d2];
                if (node1 > node2) {
                    continue;
                }

                const int lid1 = lagr_to_proteus_corners[d1];
                const int lid2 = lagr_to_proteus_corners[mixed_hex8_edge_connectivity[d1][d2]];

                int start[3], len[3], dir[3];
                for (int d = 0; d < 3; d++) {
                    start[d] = coords[d][lid1];
                }
                for (int d = 0; d < 3; d++) {
                    int x  = coords[d][lid2] - coords[d][lid1];
                    dir[d] = 1;
                    len[d] = 1;
                    if (x > 0) {
                        x -= 1;
                        len[d]   = x;
                        start[d] = 1;
                    } else if (x < 0) {
                        x += 1;
                        len[d]   = x;
                        dir[d]   = -1;
                        start[d] = L - 1;
                    }
                }

                const idx_t edge_start = index_base + g_edges[d2] * static_cast<idx_t>(nxedge);
                int         en         = 0;
                for (int zi = 0; zi != len[2]; zi += dir[2]) {
                    for (int yi = 0; yi != len[1]; yi += dir[1]) {
                        for (int xi = 0; xi != len[0]; xi += dir[0]) {
                            const int lidx_edge     = sshex8_lidx(L, start[0] + xi, start[1] + yi, start[2] + zi);
                            elements[lidx_edge][e] = edge_start + en;
                            en += 1;
                        }
                    }
                }
            }
        }
    }

    template <typename idx_t>
    static void mixed_write_tet_owned_edges(const int                                               L,
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
        idx_t nodes[4];
        for (int d = 0; d < 4; d++) {
            nodes[d] = m_elements[d][e];
        }

        for (int d1 = 0; d1 < 4; d1++) {
            const idx_t        node1     = nodes[d1];
            const idx_t *const columns   = &colidx[rowptr[node1]];
            const idx_t *const edge_view = &edge_idx[rowptr[node1]];

            idx_t g_neigh[3];
            for (int k = 0; k < 3; k++) {
                g_neigh[k] = nodes[mixed_tet4_edge_connectivity[d1][k]];
            }

            idx_t offsets[3];
            find<3>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);

            for (int d2 = 0; d2 < 3; d2++) {
                const idx_t node2 = g_neigh[d2];
                if (node1 > node2) {
                    continue;
                }

                const int lid1 = lagr_to_proteus_corners[d1];
                const int lid2 = lagr_to_proteus_corners[mixed_tet4_edge_connectivity[d1][d2]];

                int P1[3], P2[3];
                for (int d = 0; d < 3; d++) {
                    P1[d] = coords[d][lid1];
                    P2[d] = coords[d][lid2];
                }

                const idx_t edge_start = index_base + edge_view[offsets[d2]] * static_cast<idx_t>(nxedge);
                for (int t = 1; t < L; ++t) {
                    const int xi   = (P1[0] * (L - t) + P2[0] * t) / L;
                    const int yi   = (P1[1] * (L - t) + P2[1] * t) / L;
                    const int zi   = (P1[2] * (L - t) + P2[2] * t) / L;
                    const int lidx = sstet4_lidx(L, xi, yi, zi);
                    elements[lidx][e] = edge_start + (t - 1);
                }
            }
        }
    }

    template <typename idx_t>
    int ssmixed_hex_tet_generate_elements_blocks(
            const int                                                                     L,
            const ptrdiff_t                                                               n_blocks,
            const enum ElemType *const                                                    block_types,
            const ptrdiff_t *const                                                        n_elements,
            const ptrdiff_t                                                               m_nnodes,
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
            idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
            const element_idx_t *const *const                                             hft,
            const block_idx_t *const *const                                               hnbb,
            ptrdiff_t                                                                    *n_unique_nodes_out,
            ptrdiff_t                                                                    *interior_start_out) {
        SMESH_ASSERT(L >= 1);
        SMESH_ASSERT(n_blocks >= 1);
        SMESH_ASSERT(block_types != nullptr);

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

        const int nxe_hex = sshex8_nxe(L);
        const int nxe_tet = sstet4_nxe(L);

        int *hex_coords[3];
        int *tet_coords[3];
        for (int d = 0; d < 3; d++) {
            hex_coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe_hex) * sizeof(int));
            tet_coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe_tet) * sizeof(int));
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

        ptrdiff_t n_e_total = 0;
        ptrdiff_t n_e_hex   = 0;
        ptrdiff_t n_e_tet   = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            n_e_total += n_elements[b];
            if (block_types[b] == HEX8) {
                n_e_hex += n_elements[b];
                for (int d = 0; d < 8; d++) {
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        elements[b][hex_corners[d]][e] = m_elements[b][d][e];
                    }
                }
            } else if (block_types[b] == TET4) {
                n_e_tet += n_elements[b];
                for (int d = 0; d < 4; d++) {
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        elements[b][tet_corners[d]][e] = m_elements[b][d][e];
                    }
                }
            } else {
                SMESH_ERROR("ssmixed: unsupported block type %s\n", type_to_string(block_types[b]));
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
            mixed_hex_tet_build_edge_graph_from_n2e<idx_t, idx_t, idx_t>(n_blocks,
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

            ptrdiff_t edge_count = 0;
            idx_t     next_id    = 0;
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
                if (block_types[b] != HEX8) {
                    continue;
                }
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    mixed_write_hex_owned_edges(L,
                                                e,
                                                m_elements[b],
                                                rowptr,
                                                colidx,
                                                edge_idx,
                                                hex_corners,
                                                hex_coords,
                                                index_base,
                                                nxedge,
                                                elements[b]);
                }
            }
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                if (block_types[b] != TET4) {
                    continue;
                }
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    mixed_write_tet_owned_edges(L,
                                                e,
                                                m_elements[b],
                                                rowptr,
                                                colidx,
                                                edge_idx,
                                                tet_corners,
                                                tet_coords,
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

            index_base += static_cast<idx_t>(nedges * nxedge);
        }

        const int nxf_hex = (L - 1) * (L - 1);
        if (nxf_hex) {
            LocalSideTable lst;
            lst.fill(HEX8);

            idx_t n_unique_faces = 0;
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                if (block_types[b] != HEX8) {
                    continue;
                }
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    for (int f = 0; f < 6; f++) {
                        const element_idx_t neigh_element = hft[b][e * 6 + f];
                        const block_idx_t   neigh_block   = hnbb[b][e * 6 + f];
                        const bool          same_family =
                                neigh_element != invalid_idx<element_idx_t>() && block_types[neigh_block] == HEX8;
                        if (same_family &&
                            (neigh_block < static_cast<block_idx_t>(b) ||
                             (neigh_block == static_cast<block_idx_t>(b) && neigh_element < static_cast<element_idx_t>(e)))) {
                            continue;
                        }

                        const idx_t global_face_offset = index_base + n_unique_faces * nxf_hex;
                        mixed_index_hex_face(L, m_elements[b], lst.table, hex_corners, hex_coords, global_face_offset, e, f, elements[b]);

                        if (same_family) {
                            int neigh_f = 6;
                            for (int f2 = 0; f2 < 6; f2++) {
                                if (hft[neigh_block][neigh_element * 6 + f2] == static_cast<element_idx_t>(e) &&
                                    hnbb[neigh_block][neigh_element * 6 + f2] == static_cast<block_idx_t>(b)) {
                                    neigh_f = f2;
                                    break;
                                }
                            }
                            SMESH_ASSERT(neigh_f != 6);
                            mixed_index_hex_face(L,
                                                 m_elements[neigh_block],
                                                 lst.table,
                                                 hex_corners,
                                                 hex_coords,
                                                 global_face_offset,
                                                 neigh_element,
                                                 neigh_f,
                                                 elements[neigh_block]);
                        }

                        n_unique_faces++;
                    }
                }
            }

            index_base += n_unique_faces * nxf_hex;
        }

        const int nxf_tet = sstet4_nxface(L);
        if (nxf_tet) {
            LocalSideTable lst;
            lst.fill(TET4);

            idx_t n_unique_faces = 0;
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                if (block_types[b] != TET4) {
                    continue;
                }
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    for (int f = 0; f < 4; f++) {
                        const element_idx_t neigh_element = hft[b][e * 4 + f];
                        const block_idx_t   neigh_block   = hnbb[b][e * 4 + f];
                        const bool          same_family =
                                neigh_element != invalid_idx<element_idx_t>() && block_types[neigh_block] == TET4;
                        if (same_family &&
                            (neigh_block < static_cast<block_idx_t>(b) ||
                             (neigh_block == static_cast<block_idx_t>(b) && neigh_element < static_cast<element_idx_t>(e)))) {
                            continue;
                        }

                        const idx_t global_face_offset = index_base + n_unique_faces * nxf_tet;
                        mixed_index_tet_face(L,
                                             m_elements[b],
                                             lst.table,
                                             tet_corners,
                                             tet_coords,
                                             global_face_offset,
                                             e,
                                             f,
                                             elements[b]);

                        if (same_family) {
                            int neigh_f = 4;
                            for (int f2 = 0; f2 < 4; f2++) {
                                if (hft[neigh_block][neigh_element * 4 + f2] == static_cast<element_idx_t>(e) &&
                                    hnbb[neigh_block][neigh_element * 4 + f2] == static_cast<block_idx_t>(b)) {
                                    neigh_f = f2;
                                    break;
                                }
                            }
                            SMESH_ASSERT(neigh_f != 4);
                            mixed_index_tet_face(L,
                                                 m_elements[neigh_block],
                                                 lst.table,
                                                 tet_corners,
                                                 tet_coords,
                                                 global_face_offset,
                                                 neigh_element,
                                                 neigh_f,
                                                 elements[neigh_block]);
                        }

                        n_unique_faces++;
                    }
                }
            }

            index_base += n_unique_faces * nxf_tet;
        }

        const int       nxvol_hex      = (L - 1) * (L - 1) * (L - 1);
        const int       nxvol_tet      = sstet4_nxvol(L);
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

        for (int d = 0; d < 3; d++) {
            SMESH_FREE(hex_coords[d]);
            SMESH_FREE(tet_coords[d]);
        }

        *n_unique_nodes_out = index_base;
        *interior_start_out = interior_start;

        const double tock = time_seconds();
        printf("Create idx (HEX8+TET4 mixed) took\t%g [s]\n", tock - tick);
        printf("#macroelements %ld, #macronodes %ld\n", n_e_total, m_nnodes);
        printf("#micronodes %ld\n", *n_unique_nodes_out);

        return SMESH_SUCCESS;
    }

    template <typename idx_t>
    int ssmixed_hex_tet_hierarchical_renumbering_blocks(
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

        auto visit_hex_corners = [&](const bool assign_sequential) {
            for (int zi = 0; zi <= 1; zi++) {
                for (int yi = 0; yi <= 1; yi++) {
                    for (int xi = 0; xi <= 1; xi++) {
                        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                            if (block_types[b] != HEX8) {
                                continue;
                            }
                            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                                const int   v   = sshex8_lidx(L, xi * L, yi * L, zi * L);
                                const idx_t idx = elements[b][v][e];
                                SMESH_ASSERT(idx < nnodes);
                                if (assign_sequential) {
                                    if (node_mapping[idx] == invalid_idx<idx_t>()) {
                                        node_mapping[idx] = next_id++;
                                    }
                                } else {
                                    node_mapping[idx] = idx;
                                    next_id           = std::max(next_id, idx);
                                }
                            }
                        }
                    }
                }
            }
        };

        auto visit_tet_corners = [&](const bool assign_sequential) {
            for (int c = 0; c < 4; c++) {
                const int v = sstet4_lidx(L, tet_corner_xyz[c][0] * L, tet_corner_xyz[c][1] * L, tet_corner_xyz[c][2] * L);
                for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                    if (block_types[b] != TET4) {
                        continue;
                    }
                    for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                        const idx_t idx = elements[b][v][e];
                        SMESH_ASSERT(idx < nnodes);
                        if (assign_sequential) {
                            if (node_mapping[idx] == invalid_idx<idx_t>()) {
                                node_mapping[idx] = next_id++;
                            }
                        } else {
                            node_mapping[idx] = idx;
                            next_id           = std::max(next_id, idx);
                        }
                    }
                }
            }
        };

        if (preserve_corner_ordering) {
            visit_hex_corners(false);
            visit_tet_corners(false);
            next_id++;
        } else {
            visit_hex_corners(true);
            visit_tet_corners(true);
        }

        for (int k = 1; k < nlevels; k++) {
            const int l           = levels[k];
            const int step_factor = L / l;

            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                if (block_types[b] != HEX8) {
                    continue;
                }
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    for (int zi = 0; zi <= l; zi++) {
                        for (int yi = 0; yi <= l; yi++) {
                            for (int xi = 0; xi <= l; xi++) {
                                const int   v   = sshex8_lidx(L, xi * step_factor, yi * step_factor, zi * step_factor);
                                const idx_t idx = elements[b][v][e];
                                SMESH_ASSERT(idx < nnodes);
                                if (node_mapping[idx] == invalid_idx<idx_t>()) {
                                    node_mapping[idx] = next_id++;
                                }
                            }
                        }
                    }
                }
            }

            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                if (block_types[b] != TET4) {
                    continue;
                }
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    for (int z = 0; z <= l; ++z) {
                        for (int y = 0; y <= l - z; ++y) {
                            for (int x = 0; x <= l - z - y; ++x) {
                                const int   v   = sstet4_lidx(L, x * step_factor, y * step_factor, z * step_factor);
                                const idx_t idx = elements[b][v][e];
                                SMESH_ASSERT(idx < nnodes);
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
            if (block_types[b] == HEX8) {
                for (int zi = 0; zi <= L; zi++) {
                    for (int yi = 0; yi <= L; yi++) {
                        for (int xi = 0; xi <= L; xi++) {
                            const int v = sshex8_lidx(L, xi, yi, zi);
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
            } else {
                for (int z = 0; z <= L; ++z) {
                    for (int y = 0; y <= L - z; ++y) {
                        for (int x = 0; x <= L - z - y; ++x) {
                            const int v = sstet4_lidx(L, x, y, z);
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
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif
