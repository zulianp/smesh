#ifndef SMESH_SSQUAD4_GRAPH_IMPL_HPP
#define SMESH_SSQUAD4_GRAPH_IMPL_HPP

#include "smesh_alloc.hpp"
#include "smesh_sort.hpp"
#include "smesh_ssquad4_graph.hpp"

namespace smesh {

static int quad4_edge_connectivity[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

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
          idx_t neighnode = elems[quad4_edge_connectivity[lidx][d]][eidx];
          SMESH_ASSERT(nneighs < 2048);

          if (node < neighnode) {
            n2nbuff[nneighs++] = neighnode;
          }
        }
      }

      nneighs = sort_and_unique(n2nbuff, nneighs);
      rowptr[node + 1] = nneighs;
    }

    // Cumulative sum
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
          idx_t neighnode = elems[quad4_edge_connectivity[lidx][d]][eidx];
          SMESH_ASSERT(nneighs < 2048);

          if (node < neighnode) {
            n2nbuff[nneighs++] = neighnode;
          }
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

} // namespace smesh

#endif // SMESH_SSQUAD4_GRAPH_IMPL_HPP
