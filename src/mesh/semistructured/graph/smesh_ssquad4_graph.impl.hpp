#ifndef SMESH_SSQUAD4_GRAPH_IMPL_HPP
#define SMESH_SSQUAD4_GRAPH_IMPL_HPP

#include "smesh_alloc.hpp"
#include "smesh_graph.hpp"
#include "smesh_multiblock_graph.hpp"
#include "smesh_search.hpp"
#include "smesh_sort.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_graph.hpp"
#include "smesh_types.hpp"

#include <algorithm>

namespace smesh {

static int quad4_lagr_conn[4][2] = {{1, 3}, {0, 2}, {1, 3}, {0, 2}};

template <typename idx_t, typename count_t, typename element_idx_t>
int quad4_build_edge_graph_from_n2e(
    const ptrdiff_t nelements, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex, count_t **out_rowptr,
    idx_t **out_colidx) {
  SMESH_UNUSED(nelements);
  count_t *rowptr = (count_t *)SMESH_ALLOC((nnodes + 1) * sizeof(count_t));
  idx_t *colidx = 0;

  static const int nnodesxelem = 4;

  {
    rowptr[0] = 0;

#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
      idx_t n2nbuff[2048];

      count_t ebegin = n2eptr[node];
      count_t eend = n2eptr[node + 1];

      idx_t nneighs = 0;

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
        SMESH_ASSERT(lidx < 4);

        for (int d = 0; d < 2; d++) {
          idx_t neighnode = elems[quad4_lagr_conn[lidx][d]][eidx];
          SMESH_ASSERT(nneighs < 2048);

          n2nbuff[nneighs++] = neighnode;
        }
      }

      nneighs = sort_and_unique(n2nbuff, nneighs);
      rowptr[node + 1] = nneighs;
    }

    for (ptrdiff_t node = 0; node < nnodes; ++node) {
      rowptr[node + 1] += rowptr[node];
    }

    const ptrdiff_t nnz = rowptr[nnodes];
    colidx = (idx_t *)SMESH_ALLOC(nnz * sizeof(idx_t));

#pragma omp parallel for
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
      idx_t n2nbuff[2048];

      count_t ebegin = n2eptr[node];
      count_t eend = n2eptr[node + 1];

      idx_t nneighs = 0;

      for (count_t e = ebegin; e < eend; ++e) {
        element_idx_t eidx = elindex[e];
        SMESH_ASSERT(eidx < nelements);

        int lidx = 0;
        for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
          if (elems[edof_i][eidx] == node) {
            lidx = edof_i;
            break;
          }
        }

        for (int d = 0; d < 2; d++) {
          idx_t neighnode = elems[quad4_lagr_conn[lidx][d]][eidx];
          SMESH_ASSERT(nneighs < 2048);

          n2nbuff[nneighs++] = neighnode;
        }
      }

      nneighs = sort_and_unique(n2nbuff, nneighs);

      for (idx_t i = 0; i < nneighs; ++i) {
        colidx[rowptr[node] + i] = n2nbuff[i];
      }
    }
  }

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return 0;
}

template <typename idx_t, typename count_t, typename element_idx_t>
static int quad4_build_edge_graph_from_multiblock_n2e(
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

            int lidx = -1;
            for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
                if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
                    lidx = edof_i;
                    break;
                }
            }
            SMESH_ASSERT(lidx != -1);

            for (int d = 0; d < 2; d++) {
                SMESH_ASSERT(nneighs < 2048);
                const idx_t neighnode = block_elems[quad4_lagr_conn[lidx][d]][eidx];
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
        idx_t   n2nbuff[2048];
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
            for (int d = 0; d < 2; d++) {
                const idx_t neighnode = block_elems[quad4_lagr_conn[lidx][d]][eidx];
                n2nbuff[nneighs++] = neighnode;
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
static int quad4_build_edge_graph(const ptrdiff_t                                         nelements,
                                  const ptrdiff_t                                         nnodes,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                  count_t                                               **out_rowptr,
                                  idx_t                                                 **out_colidx) {
    count_t       *n2eptr;
    element_idx_t *elindex;
    create_n2e<idx_t, count_t, element_idx_t>(nelements, nnodes, 4, elems, &n2eptr, &elindex);
    int err = quad4_build_edge_graph_from_n2e(nelements, nnodes, elems, n2eptr, elindex, out_rowptr, out_colidx);
    SMESH_FREE(n2eptr);
    SMESH_FREE(elindex);
    return err;
}

static void ssquad4_fill_lattice_coords(const int L, int *coords[2]) {
    const int nxe = ssquad4_nxe(L);
    for (int d = 0; d < 2; ++d) {
        coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(int));
    }
    for (int yi = 0; yi <= L; ++yi) {
        for (int xi = 0; xi <= L; ++xi) {
            const int lidx    = ssquad4_lidx(L, xi, yi);
            coords[0][lidx] = xi;
            coords[1][lidx] = yi;
        }
    }
}

template <typename idx_t>
static void ssquad4_index_owned_edges(const int                                               L,
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

        idx_t g_neigh[2];
        for (int k = 0; k < 2; k++) {
            g_neigh[k] = nodes[quad4_lagr_conn[d1][k]];
        }

        idx_t offsets[2];
        find<2>(g_neigh, columns, static_cast<int>(rowptr[node1 + 1] - rowptr[node1]), offsets);

        for (int d2 = 0; d2 < 2; d2++) {
            const idx_t node2 = g_neigh[d2];
            if (node1 > node2) {
                continue;
            }

            const int lid1 = lagr_to_proteus_corners[d1];
            const int lid2 = lagr_to_proteus_corners[quad4_lagr_conn[d1][d2]];

            int P1[2], P2[2];
            for (int d = 0; d < 2; d++) {
                P1[d] = coords[d][lid1];
                P2[d] = coords[d][lid2];
            }

            const idx_t edge_start = index_base + edge_view[offsets[d2]] * static_cast<idx_t>(nxedge);
            for (int t = 1; t < L; ++t) {
                const int xi   = (P1[0] * (L - t) + P2[0] * t) / L;
                const int yi   = (P1[1] * (L - t) + P2[1] * t) / L;
                const int lidx = ssquad4_lidx(L, xi, yi);
                elements[lidx][e] = edge_start + (t - 1);
            }
        }
    }
}

template <typename idx_t>
static void ssquad4_index_owned_interiors(const int                                         L,
                                          const ptrdiff_t                                   e,
                                          const idx_t                                       index_base,
                                          const ptrdiff_t                                   e_global,
                                          const int                                         nxf,
                                          idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements) {
    const int Lm1 = L - 1;
    for (int yi = 1; yi < L; ++yi) {
        for (int xi = 1; xi < L; ++xi) {
            const int lidx            = ssquad4_lidx(L, xi, yi);
            const int en              = (yi - 1) * Lm1 + (xi - 1);
            elements[lidx][e]         = index_base + static_cast<idx_t>(e_global) * static_cast<idx_t>(nxf) +
                                static_cast<idx_t>(en);
        }
    }
}

template <typename idx_t>
int ssquad4_generate_elements(const int                                               L,
                              const ptrdiff_t                                         m_nelements,
                              const ptrdiff_t                                         m_nnodes,
                              const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                              idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
                              ptrdiff_t                                              *n_unique_nodes_out,
                              ptrdiff_t                                              *interior_start_out) {
    SMESH_ASSERT(L >= 1);

    int lagr_to_proteus_corners[4] = {ssquad4_lidx(L, 0, 0),
                                      ssquad4_lidx(L, L, 0),
                                      ssquad4_lidx(L, L, L),
                                      ssquad4_lidx(L, 0, L)};

    int *coords[2];
    ssquad4_fill_lattice_coords(L, coords);

    for (int d = 0; d < 4; d++) {
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            elements[lagr_to_proteus_corners[d]][e] = m_elements[d][e];
        }
    }

    idx_t     index_base = static_cast<idx_t>(m_nnodes);
    const int nxedge     = ssquad4_nxedge(L);

    if (nxedge) {
        idx_t *rowptr{nullptr};
        idx_t *colidx{nullptr};
        quad4_build_edge_graph<idx_t, idx_t, idx_t>(m_nelements, m_nnodes, m_elements, &rowptr, &colidx);

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
            ssquad4_index_owned_edges(L,
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

    const int       nxf            = ssquad4_nxface(L);
    const ptrdiff_t interior_start = index_base;
    if (nxf) {
        for (ptrdiff_t e = 0; e < m_nelements; e++) {
            ssquad4_index_owned_interiors(L, e, index_base, e, nxf, elements);
        }
        index_base += static_cast<idx_t>(m_nelements) * static_cast<idx_t>(nxf);
    }

    for (int d = 0; d < 2; d++) {
        SMESH_FREE(coords[d]);
    }

    *n_unique_nodes_out = index_base;
    *interior_start_out = interior_start;
    return SMESH_SUCCESS;
}

template <typename idx_t>
int ssquad4_generate_elements_blocks(
        const int                                                                     L,
        const ptrdiff_t                                                               n_blocks,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               m_nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        ptrdiff_t                                                                    *n_unique_nodes_out,
        ptrdiff_t                                                                    *interior_start_out) {
    SMESH_ASSERT(L >= 1);
    SMESH_ASSERT(n_blocks >= 1);

    int lagr_to_proteus_corners[4] = {ssquad4_lidx(L, 0, 0),
                                      ssquad4_lidx(L, L, 0),
                                      ssquad4_lidx(L, L, L),
                                      ssquad4_lidx(L, 0, L)};

    int *coords[2];
    ssquad4_fill_lattice_coords(L, coords);

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
    const int nxedge     = ssquad4_nxedge(L);

    if (nxedge) {
        enum ElemType *types = (enum ElemType *)SMESH_ALLOC(static_cast<size_t>(n_blocks) * sizeof(enum ElemType));
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            types[b] = QUAD4;
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
        quad4_build_edge_graph_from_multiblock_n2e<idx_t, idx_t, idx_t>(
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
                ssquad4_index_owned_edges(L,
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

    const int       nxf            = ssquad4_nxface(L);
    const ptrdiff_t interior_start = index_base;
    if (nxf) {
        ptrdiff_t e_off = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                ssquad4_index_owned_interiors(L, e, index_base, e_off + e, nxf, elements[b]);
            }
            e_off += n_elements[b];
        }
        index_base += static_cast<idx_t>(n_e_total) * static_cast<idx_t>(nxf);
    }

    for (int d = 0; d < 2; d++) {
        SMESH_FREE(coords[d]);
    }

    *n_unique_nodes_out = index_base;
    *interior_start_out = interior_start;
    return SMESH_SUCCESS;
}

template <typename idx_t>
int ssquad4_hierarchical_renumbering(const int                                         L,
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

    for (int yi = 0; yi <= 1; yi++) {
        for (int xi = 0; xi <= 1; xi++) {
            for (ptrdiff_t e = 0; e < nelements; e++) {
                const int   v   = ssquad4_lidx(L, xi * L, yi * L);
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
    }

    if (preserve_corner_ordering) {
        next_id++;
    }

    int stride = 1;
    for (int k = 1; k < nlevels; k++) {
        const int l           = levels[k];
        const int step_factor = L / l;
        for (ptrdiff_t e = 0; e < nelements; e++) {
            for (int yi = 0; yi <= l; yi += stride) {
                for (int xi = 0; xi <= l; xi += stride) {
                    const int   v   = ssquad4_lidx(L, xi * step_factor, yi * step_factor);
                    const idx_t idx = elements[v][e];
                    SMESH_ASSERT(idx < nnodes);
                    if (node_mapping[idx] == invalid_idx<idx_t>()) {
                        node_mapping[idx] = next_id++;
                    }
                }
            }
        }
    }

    for (int yi = 0; yi <= L; yi++) {
        for (int xi = 0; xi <= L; xi++) {
            for (ptrdiff_t e = 0; e < nelements; e++) {
                const int   v   = ssquad4_lidx(L, xi, yi);
                const idx_t idx = elements[v][e];
                if (node_mapping[idx] == invalid_idx<idx_t>()) {
                    SMESH_ERROR("Uninitialized node mapping\n");
                }
                elements[v][e] = node_mapping[idx];
            }
        }
    }

    return SMESH_SUCCESS;
}

template <typename idx_t>
int ssquad4_hierarchical_renumbering_blocks(
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

    for (int yi = 0; yi <= 1; yi++) {
        for (int xi = 0; xi <= 1; xi++) {
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    const int   v   = ssquad4_lidx(L, xi * L, yi * L);
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
    }

    if (preserve_corner_ordering) {
        next_id++;
    }

    int stride = 1;
    for (int k = 1; k < nlevels; k++) {
        const int l           = levels[k];
        const int step_factor = L / l;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                for (int yi = 0; yi <= l; yi += stride) {
                    for (int xi = 0; xi <= l; xi += stride) {
                        const int   v   = ssquad4_lidx(L, xi * step_factor, yi * step_factor);
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

    for (int yi = 0; yi <= L; yi++) {
        for (int xi = 0; xi <= L; xi++) {
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                for (ptrdiff_t e = 0; e < n_elements[b]; e++) {
                    const int   v   = ssquad4_lidx(L, xi, yi);
                    const idx_t idx = elements[b][v][e];
                    if (node_mapping[idx] == invalid_idx<idx_t>()) {
                        SMESH_ERROR("Uninitialized node mapping\n");
                    }
                    elements[b][v][e] = node_mapping[idx];
                }
            }
        }
    }

    return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int ssquad4_extract_nodeset_from_sideset(const int                                               L,
                                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                         const ptrdiff_t                                         n_surf_elements,
                                         const element_idx_t *const SMESH_RESTRICT               parent_element,
                                         const i16 *const SMESH_RESTRICT                         side_idx,
                                         ptrdiff_t                                              *n_nodes_out,
                                         idx_t **SMESH_RESTRICT                                  nodes_out) {
    const int       nnxs = L + 1;
    const ptrdiff_t n    = static_cast<ptrdiff_t>(nnxs) * n_surf_elements;
    if (n == 0) {
        *n_nodes_out = 0;
        *nodes_out   = nullptr;
        return SMESH_SUCCESS;
    }

    idx_t *nodes = static_cast<idx_t *>(SMESH_ALLOC(static_cast<size_t>(n) * sizeof(idx_t)));

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n_surf_elements; i++) {
        const ptrdiff_t e = parent_element[i];
        const int       s = side_idx[i];
        SMESH_ASSERT(s >= 0 && s < 4);
        idx_t *const out = &nodes[i * nnxs];
        for (int t = 0; t <= L; ++t) {
            int x = 0;
            int y = 0;
            switch (s) {
                case 0:
                    x = t;
                    y = 0;
                    break;
                case 1:
                    x = L;
                    y = t;
                    break;
                case 2:
                    x = L - t;
                    y = L;
                    break;
                case 3:
                    x = 0;
                    y = L - t;
                    break;
                default:
                    break;
            }
            out[t] = elems[ssquad4_lidx(L, x, y)][e];
        }
    }

    *n_nodes_out = static_cast<ptrdiff_t>(sort_and_unique(nodes, static_cast<size_t>(n)));
    *nodes_out   = static_cast<idx_t *>(SMESH_REALLOC(nodes, static_cast<size_t>(*n_nodes_out) * sizeof(idx_t)));
    return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int ssquad4_extract_surface_from_sideset(const int                                               L,
                                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                         const ptrdiff_t                                         n_surf_elements,
                                         const element_idx_t *const SMESH_RESTRICT               parent_element,
                                         const i16 *const SMESH_RESTRICT                         side_idx,
                                         idx_t **const SMESH_RESTRICT                            sides) {
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n_surf_elements; i++) {
        const ptrdiff_t e = parent_element[i];
        const int       s = side_idx[i];
        SMESH_ASSERT(s >= 0 && s < 4);
        for (int t = 0; t <= L; ++t) {
            int x = 0;
            int y = 0;
            switch (s) {
                case 0:
                    x = t;
                    y = 0;
                    break;
                case 1:
                    x = L;
                    y = t;
                    break;
                case 2:
                    x = L - t;
                    y = L;
                    break;
                case 3:
                    x = 0;
                    y = L - t;
                    break;
                default:
                    break;
            }
            sides[ssedge_lidx(L, t)][i] = elems[ssquad4_lidx(L, x, y)][e];
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t>
int ssedge_hierarchical_remapping(const int                                         L,
                                  const int                                         nlevels,
                                  int *const                                        levels,
                                  const ptrdiff_t                                   nelements,
                                  const ptrdiff_t                                   nnodes,
                                  idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                  idx_t **SMESH_RESTRICT                            node_mapping_out,
                                  ptrdiff_t                                        *count_out) {
    idx_t *node_mapping = (idx_t *)SMESH_ALLOC(nnodes * sizeof(idx_t));

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < nnodes; i++) {
        node_mapping[i] = invalid_idx<idx_t>();
    }

    idx_t next_id = 0;

    const int corners[2] = {0, L};
    for (int c = 0; c < 2; ++c) {
        for (ptrdiff_t e = 0; e < nelements; e++) {
            const idx_t idx = elements[ssedge_lidx(L, corners[c])][e];
            SMESH_ASSERT(idx < nnodes);
            if (node_mapping[idx] == invalid_idx<idx_t>()) {
                node_mapping[idx] = next_id++;
            }
        }
    }

    for (int k = 1; k < nlevels; k++) {
        const int l           = levels[k];
        const int step_factor = L / l;
        for (ptrdiff_t e = 0; e < nelements; e++) {
            for (int xi = 0; xi <= l; xi++) {
                const idx_t idx = elements[ssedge_lidx(L, xi * step_factor)][e];
                SMESH_ASSERT(idx < nnodes);
                if (node_mapping[idx] == invalid_idx<idx_t>()) {
                    node_mapping[idx] = next_id++;
                }
            }
        }
    }

    for (int xi = 0; xi <= L; xi++) {
        for (ptrdiff_t e = 0; e < nelements; e++) {
            const idx_t idx = elements[ssedge_lidx(L, xi)][e];
            if (node_mapping[idx] == invalid_idx<idx_t>()) {
                SMESH_ERROR("Uninitialized node mapping\n");
            }
            elements[ssedge_lidx(L, xi)][e] = node_mapping[idx];
        }
    }

    *count_out        = next_id;
    *node_mapping_out = (idx_t *)SMESH_ALLOC(*count_out * sizeof(idx_t));

#ifndef NDEBUG
    for (ptrdiff_t i = 0; i < *count_out; i++) {
        (*node_mapping_out)[i] = invalid_idx<idx_t>();
    }
#endif

    for (ptrdiff_t i = 0; i < nnodes; i++) {
        if (node_mapping[i] != invalid_idx<idx_t>()) {
            (*node_mapping_out)[node_mapping[i]] = i;
        }
    }

#ifndef NDEBUG
    for (ptrdiff_t i = 0; i < *count_out; i++) {
        SMESH_ASSERT((*node_mapping_out)[i] != invalid_idx<idx_t>());
    }
#endif

    SMESH_FREE(node_mapping);
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif

