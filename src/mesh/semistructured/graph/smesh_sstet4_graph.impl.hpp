#ifndef SMESH_SSTET4_GRAPH_IMPL_HPP
#define SMESH_SSTET4_GRAPH_IMPL_HPP

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_graph.hpp"
#include "smesh_multiblock_graph.hpp"
#include "smesh_search.hpp"
#include "smesh_sort.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sstet4_graph.hpp"

#include <algorithm>

namespace smesh {

    static int tet4_edge_connectivity[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

    template <typename idx_t, typename count_t, typename element_idx_t>
    static int tet4_build_edge_graph_from_n2e(const ptrdiff_t                                         nelements,
                                              const ptrdiff_t                                         nnodes,
                                              const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                              const count_t *const SMESH_RESTRICT                     n2eptr,
                                              const element_idx_t *const SMESH_RESTRICT               elindex,
                                              count_t                                               **out_rowptr,
                                              idx_t                                                 **out_colidx) {
        SMESH_UNUSED(nelements);
        count_t *rowptr = (count_t *)SMESH_ALLOC((nnodes + 1) * sizeof(count_t));
        idx_t   *colidx = 0;
        static const int nnodesxelem = 4;

        rowptr[0] = 0;

#pragma omp parallel for
        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            idx_t n2nbuff[2048];
            count_t ebegin = n2eptr[node];
            count_t eend   = n2eptr[node + 1];
            idx_t   nneighs = 0;

            for (count_t e = ebegin; e < eend; ++e) {
                element_idx_t eidx = elindex[e];
                SMESH_ASSERT(eidx < nelements);

                int lidx = -1;
                for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
                    if (elems[edof_i][eidx] == node) {
                        lidx = edof_i;
                        break;
                    }
                }
                SMESH_ASSERT(lidx != -1);

                for (int d = 0; d < 3; d++) {
                    idx_t neighnode = elems[tet4_edge_connectivity[lidx][d]][eidx];
                    SMESH_ASSERT(nneighs < 2048);
                    n2nbuff[nneighs++] = neighnode;
                }
            }

            nneighs          = sort_and_unique(n2nbuff, nneighs);
            rowptr[node + 1] = nneighs;
        }

        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            rowptr[node + 1] += rowptr[node];
        }

        const ptrdiff_t nnz = rowptr[nnodes];
        colidx              = (idx_t *)SMESH_ALLOC(nnz * sizeof(idx_t));

#pragma omp parallel for
        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            idx_t n2nbuff[2048];
            count_t ebegin  = n2eptr[node];
            count_t eend    = n2eptr[node + 1];
            idx_t   nneighs = 0;

            for (count_t e = ebegin; e < eend; ++e) {
                element_idx_t eidx = elindex[e];
                int           lidx = 0;
                for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
                    if (elems[edof_i][eidx] == node) {
                        lidx = edof_i;
                        break;
                    }
                }
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = elems[tet4_edge_connectivity[lidx][d]][eidx];
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
    static int tet4_build_edge_graph_from_multiblock_n2e(
            const ptrdiff_t                                                                n_blocks,
            const ptrdiff_t *const                                                         n_elements,
            const ptrdiff_t                                                                nnodes,
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
            const count_t *const SMESH_RESTRICT                                            n2eptr,
            const element_idx_t *const SMESH_RESTRICT                                      elindex,
            const block_idx_t *const SMESH_RESTRICT                                        block_number,
            count_t                                                                      **out_rowptr,
            idx_t                                                                        **out_colidx) {
        SMESH_UNUSED(n_elements);
        count_t *rowptr = (count_t *)SMESH_ALLOC((nnodes + 1) * sizeof(count_t));
        idx_t   *colidx = 0;
        static const int nnodesxelem = 4;

        rowptr[0] = 0;

#pragma omp parallel for
        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            idx_t n2nbuff[2048];
            count_t ebegin  = n2eptr[node];
            count_t eend    = n2eptr[node + 1];
            idx_t   nneighs = 0;

            for (count_t e = ebegin; e < eend; ++e) {
                const block_idx_t   b    = block_number[e];
                const element_idx_t eidx = elindex[e];
                SMESH_ASSERT(b >= 0 && b < n_blocks);
                SMESH_ASSERT(eidx < n_elements[b]);
                const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];

                int lidx = -1;
                for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
                    if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                        lidx = edof_i;
                        break;
                    }
                }
                SMESH_ASSERT(lidx != -1);

                for (int d = 0; d < 3; d++) {
                    SMESH_ASSERT(nneighs < 2048);
                    n2nbuff[nneighs++] = block_elems[tet4_edge_connectivity[lidx][d]][eidx];
                }
            }

            nneighs          = sort_and_unique(n2nbuff, nneighs);
            rowptr[node + 1] = nneighs;
        }

        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            rowptr[node + 1] += rowptr[node];
        }

        const ptrdiff_t nnz = rowptr[nnodes];
        colidx              = (idx_t *)SMESH_ALLOC(nnz * sizeof(idx_t));

#pragma omp parallel for
        for (ptrdiff_t node = 0; node < nnodes; ++node) {
            idx_t n2nbuff[2048];
            count_t ebegin  = n2eptr[node];
            count_t eend    = n2eptr[node + 1];
            idx_t   nneighs = 0;

            for (count_t e = ebegin; e < eend; ++e) {
                const block_idx_t   b    = block_number[e];
                const element_idx_t eidx = elindex[e];
                const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems = elems[b];

                int lidx = 0;
                for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
                    if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                        lidx = edof_i;
                        break;
                    }
                }
                for (int d = 0; d < 3; d++) {
                    n2nbuff[nneighs++] = block_elems[tet4_edge_connectivity[lidx][d]][eidx];
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

    template <typename element_idx_t, typename idx_t, typename count_t>
    static int tet4_build_edge_graph(const ptrdiff_t                                         nelements,
                                     const ptrdiff_t                                         nnodes,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                     count_t                                               **out_rowptr,
                                     idx_t                                                 **out_colidx) {
        count_t       *n2eptr;
        element_idx_t *elindex;
        create_n2e<idx_t, count_t, element_idx_t>(nelements, nnodes, 4, elems, &n2eptr, &elindex);
        int err = tet4_build_edge_graph_from_n2e(nelements, nnodes, elems, n2eptr, elindex, out_rowptr, out_colidx);
        SMESH_FREE(n2eptr);
        SMESH_FREE(elindex);
        return err;
    }

    static void sstet4_fill_lattice_coords(const int L, int **coords) {
        const int nxe = sstet4_nxe(L);
        for (int d = 0; d < 3; d++) {
            coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(int));
        }
        for (int z = 0; z <= L; ++z) {
            for (int y = 0; y <= L - z; ++y) {
                for (int x = 0; x <= L - z - y; ++x) {
                    const int lidx = sstet4_lidx(L, x, y, z);
                    coords[0][lidx] = x;
                    coords[1][lidx] = y;
                    coords[2][lidx] = z;
                }
            }
        }
    }

    template <typename idx_t>
    static void sstet4_index_face(const int                                               L,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                  const int *const                                        local_side_table,
                                  const int *const                                        lagr_to_proteus_corners,
                                  int *const *const                                       coords,
                                  const idx_t                                             global_face_offset,
                                  const ptrdiff_t                                         e,
                                  const int                                               f,
                                  idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements) {
        static const int nnxs = LocalSideTable::MAX_NUM_NODES_PER_SIDE;
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
                const int w  = L - s - t;
                const int xi = (Po[0] * w + Pu[0] * s + Pv[0] * t) / L;
                const int yi = (Po[1] * w + Pu[1] * s + Pv[1] * t) / L;
                const int zi = (Po[2] * w + Pu[2] * s + Pv[2] * t) / L;
                const int pidx = sstet4_lidx(L, xi, yi, zi);
                elements[pidx][e] = global_face_offset + local_offset++;
            }
        }
    }

    template <typename idx_t>
    static void sstet4_index_owned_edges(const int                                               L,
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
                g_neigh[k] = nodes[tet4_edge_connectivity[d1][k]];
            }

            idx_t offsets[3];
            find<3>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);

            for (int d2 = 0; d2 < 3; d2++) {
                const idx_t node2 = g_neigh[d2];
                if (node1 > node2) {
                    continue;
                }

                const int lid1 = lagr_to_proteus_corners[d1];
                const int lid2 = lagr_to_proteus_corners[tet4_edge_connectivity[d1][d2]];

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
    int sstet4_generate_elements(const int                                               L,
                                 const ptrdiff_t                                         m_nelements,
                                 const ptrdiff_t                                         m_nnodes,
                                 const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                                 idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
                                 ptrdiff_t                                              *n_unique_nodes_out,
                                 ptrdiff_t                                              *interior_start_out) {
        SMESH_ASSERT(L >= 1);

        double tick = time_seconds();
        const int nxe = sstet4_nxe(L);

        int lagr_to_proteus_corners[4] = {sstet4_lidx(L, 0, 0, 0),
                                          sstet4_lidx(L, L, 0, 0),
                                          sstet4_lidx(L, 0, L, 0),
                                          sstet4_lidx(L, 0, 0, L)};

        int *coords[3];
        sstet4_fill_lattice_coords(L, coords);

        for (int d = 0; d < 4; d++) {
            for (ptrdiff_t e = 0; e < m_nelements; e++) {
                elements[lagr_to_proteus_corners[d]][e] = m_elements[d][e];
            }
        }

        idx_t index_base = static_cast<idx_t>(m_nnodes);
        const int nxedge = sstet4_nxedge(L);

        if (nxedge) {
            idx_t *rowptr{nullptr};
            idx_t *colidx{nullptr};
            tet4_build_edge_graph<idx_t, idx_t, idx_t>(m_nelements, m_nnodes, m_elements, &rowptr, &colidx);

            const ptrdiff_t nedges   = rowptr[m_nnodes] / 2;
            const ptrdiff_t nnz      = rowptr[m_nnodes];
            idx_t          *edge_idx = (idx_t *)SMESH_CALLOC(nnz, sizeof(idx_t));

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

            for (ptrdiff_t e = 0; e < m_nelements; e++) {
                sstet4_index_owned_edges(L,
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

        const int nxf = sstet4_nxface(L);
        if (nxf) {
            LocalSideTable lst;
            lst.fill(TET4);

            element_idx_t *adj_table = 0;
            create_element_adj_table(m_nelements, m_nnodes, TET4, m_elements, &adj_table);

            idx_t n_unique_faces = 0;
            for (ptrdiff_t e = 0; e < m_nelements; e++) {
                for (int f = 0; f < 4; f++) {
                    const element_idx_t neigh_element = adj_table[e * 4 + f];
                    if (neigh_element != invalid_idx<element_idx_t>() && neigh_element < e) {
                        continue;
                    }

                    const idx_t global_face_offset = index_base + n_unique_faces * nxf;
                    sstet4_index_face(L, m_elements, lst.table, lagr_to_proteus_corners, coords, global_face_offset, e, f, elements);

                    if (neigh_element != invalid_idx<element_idx_t>()) {
                        int neigh_f = 4;
                        for (int f2 = 0; f2 < 4; f2++) {
                            if (e == adj_table[neigh_element * 4 + f2]) {
                                neigh_f = f2;
                                break;
                            }
                        }
                        SMESH_ASSERT(neigh_f != 4);
                        sstet4_index_face(L,
                                          m_elements,
                                          lst.table,
                                          lagr_to_proteus_corners,
                                          coords,
                                          global_face_offset,
                                          neigh_element,
                                          neigh_f,
                                          elements);
                    }
                    n_unique_faces++;
                }
            }

            index_base += n_unique_faces * nxf;
            SMESH_FREE(adj_table);
        }

        const int       nxvol          = sstet4_nxvol(L);
        const ptrdiff_t interior_start = index_base;
        if (nxvol) {
            for (int z = 1; z <= L - 3; ++z) {
                for (int y = 1; y <= L - 2 - z; ++y) {
                    for (int x = 1; x <= L - 1 - z - y; ++x) {
                        const int lidx_vol = sstet4_lidx(L, x, y, z);
                        const int en = sstet4_lidx(L - 4, x - 1, y - 1, z - 1);
                        for (ptrdiff_t e = 0; e < m_nelements; e++) {
                            elements[lidx_vol][e] = index_base + e * nxvol + en;
                        }
                    }
                }
            }
        }

        for (int d = 0; d < 3; d++) {
            SMESH_FREE(coords[d]);
        }

        *n_unique_nodes_out = interior_start + m_nelements * nxvol;
        *interior_start_out = interior_start;

        const double tock = time_seconds();
        printf("Create idx (%s) took\t%g [s]\n", type_to_string(TET4), tock - tick);
        printf("#macroelements %ld, #macronodes %ld\n", m_nelements, m_nnodes);
        printf("#microelements %ld, #micronodes %ld\n", m_nelements * sstet4_txe(L), *n_unique_nodes_out);

        SMESH_UNUSED(nxe);
        return SMESH_SUCCESS;
    }

    template <typename idx_t>
    int sstet4_generate_elements_blocks(
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
        SMESH_ASSERT(n_blocks >= 1);

        double tick = time_seconds();

        int lagr_to_proteus_corners[4] = {sstet4_lidx(L, 0, 0, 0),
                                          sstet4_lidx(L, L, 0, 0),
                                          sstet4_lidx(L, 0, L, 0),
                                          sstet4_lidx(L, 0, 0, L)};

        int *coords[3];
        sstet4_fill_lattice_coords(L, coords);

        ptrdiff_t n_e_total = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            n_e_total += n_elements[b];
            for (int d = 0; d < 4; d++) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    elements[b][lagr_to_proteus_corners[d]][e] = m_elements[b][d][e];
                }
            }
        }

        idx_t     index_base = static_cast<idx_t>(m_nnodes);
        const int nxedge     = sstet4_nxedge(L);

        if (nxedge) {
            enum ElemType *types = (enum ElemType *)SMESH_ALLOC(static_cast<size_t>(n_blocks) * sizeof(enum ElemType));
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                types[b] = TET4;
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
            tet4_build_edge_graph_from_multiblock_n2e<idx_t, idx_t, idx_t>(
                    n_blocks, n_elements, m_nnodes, m_elements, n2eptr, elindex, block_number, &rowptr, &colidx);

            const ptrdiff_t nedges   = rowptr[m_nnodes] / 2;
            const ptrdiff_t nnz      = rowptr[m_nnodes];
            idx_t          *edge_idx = (idx_t *)SMESH_CALLOC(nnz, sizeof(idx_t));

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
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    sstet4_index_owned_edges(L,
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

        const int nxf = sstet4_nxface(L);
        if (nxf) {
            LocalSideTable lst;
            lst.fill(TET4);

            idx_t n_unique_faces = 0;
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    for (int f = 0; f < 4; f++) {
                        const element_idx_t neigh_element = hft[b][e * 4 + f];
                        const block_idx_t   neigh_block   = hnbb[b][e * 4 + f];
                        if (neigh_element != invalid_idx<element_idx_t>() &&
                            (neigh_block < static_cast<block_idx_t>(b) ||
                             (neigh_block == static_cast<block_idx_t>(b) && neigh_element < static_cast<element_idx_t>(e)))) {
                            continue;
                        }

                        const idx_t global_face_offset = index_base + n_unique_faces * nxf;
                        sstet4_index_face(L,
                                          m_elements[b],
                                          lst.table,
                                          lagr_to_proteus_corners,
                                          coords,
                                          global_face_offset,
                                          e,
                                          f,
                                          elements[b]);

                        if (neigh_element != invalid_idx<element_idx_t>()) {
                            int neigh_f = 4;
                            for (int f2 = 0; f2 < 4; f2++) {
                                if (hft[neigh_block][neigh_element * 4 + f2] == static_cast<element_idx_t>(e) &&
                                    hnbb[neigh_block][neigh_element * 4 + f2] == static_cast<block_idx_t>(b)) {
                                    neigh_f = f2;
                                    break;
                                }
                            }
                            SMESH_ASSERT(neigh_f != 4);
                            sstet4_index_face(L,
                                              m_elements[neigh_block],
                                              lst.table,
                                              lagr_to_proteus_corners,
                                              coords,
                                              global_face_offset,
                                              neigh_element,
                                              neigh_f,
                                              elements[neigh_block]);
                        }

                        n_unique_faces++;
                    }
                }
            }

            index_base += n_unique_faces * nxf;
        }

        const int       nxvol          = sstet4_nxvol(L);
        const ptrdiff_t interior_start = index_base;
        if (nxvol) {
            ptrdiff_t e_off = 0;
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                for (int z = 1; z <= L - 3; ++z) {
                    for (int y = 1; y <= L - 2 - z; ++y) {
                        for (int x = 1; x <= L - 1 - z - y; ++x) {
                            const int lidx_vol = sstet4_lidx(L, x, y, z);
                            const int en       = sstet4_lidx(L - 4, x - 1, y - 1, z - 1);
                            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                                elements[b][lidx_vol][e] = index_base + (e_off + e) * nxvol + en;
                            }
                        }
                    }
                }
                e_off += n_elements[b];
            }
        }

        for (int d = 0; d < 3; d++) {
            SMESH_FREE(coords[d]);
        }

        *n_unique_nodes_out = interior_start + n_e_total * nxvol;
        *interior_start_out = interior_start;

        const double tock = time_seconds();
        printf("Create idx (TET4 blocks) took\t%g [s]\n", tock - tick);
        printf("#macroelements %ld, #macronodes %ld\n", n_e_total, m_nnodes);
        printf("#microelements %ld, #micronodes %ld\n", n_e_total * sstet4_txe(L), *n_unique_nodes_out);

        return SMESH_SUCCESS;
    }

    static const int sstet4_corner_xyz[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    template <typename idx_t>
    int sstet4_hierarchical_renumbering(const int                                         L,
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

        for (int c = 0; c < 4; c++) {
            const int v = sstet4_lidx(L, sstet4_corner_xyz[c][0] * L, sstet4_corner_xyz[c][1] * L, sstet4_corner_xyz[c][2] * L);
            for (ptrdiff_t e = 0; e < nelements; e++) {
                const idx_t idx = elements[v][e];
                SMESH_ASSERT(idx < nnodes);
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

        for (int k = 1; k < nlevels; k++) {
            const int l           = levels[k];
            const int step_factor = L / l;
            for (ptrdiff_t e = 0; e < nelements; e++) {
                for (int z = 0; z <= l; ++z) {
                    for (int y = 0; y <= l - z; ++y) {
                        for (int x = 0; x <= l - z - y; ++x) {
                            const int   v   = sstet4_lidx(L, x * step_factor, y * step_factor, z * step_factor);
                            const idx_t idx = elements[v][e];
                            SMESH_ASSERT(idx < nnodes);
                            if (node_mapping[idx] == invalid_idx<idx_t>()) {
                                node_mapping[idx] = next_id++;
                            }
                        }
                    }
                }
            }
        }

        for (int z = 0; z <= L; ++z) {
            for (int y = 0; y <= L - z; ++y) {
                for (int x = 0; x <= L - z - y; ++x) {
                    const int v = sstet4_lidx(L, x, y, z);
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
    int sstet4_hierarchical_renumbering_blocks(
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

        for (int c = 0; c < 4; c++) {
            const int v = sstet4_lidx(L, sstet4_corner_xyz[c][0] * L, sstet4_corner_xyz[c][1] * L, sstet4_corner_xyz[c][2] * L);
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    const idx_t idx = elements[b][v][e];
                    SMESH_ASSERT(idx < nnodes);
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

        for (int k = 1; k < nlevels; k++) {
            const int l           = levels[k];
            const int step_factor = L / l;
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
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

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_SSTET4_GRAPH_IMPL_HPP
