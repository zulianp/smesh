#ifndef SMESH_MULTIBLOCK_GRAPH_IMPL_HPP
#define SMESH_MULTIBLOCK_GRAPH_IMPL_HPP

#include "smesh_multiblock_graph.hpp"
#include "smesh_sort.hpp"

#include <cstring>

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_n2e(const block_idx_t n_blocks,
                         const enum ElemType element_types[],
                         const ptrdiff_t n_elements[],
                         idx_t **const SMESH_RESTRICT elements[],
                         const ptrdiff_t n_nodes,
                         block_idx_t **out_block_number, count_t **out_n2eptr,
                         element_idx_t **out_elindex) {
  count_t *n2eptr = (count_t *)malloc((n_nodes + 1) * sizeof(count_t));
  std::memset(n2eptr, 0, (n_nodes + 1) * sizeof(count_t));

  int *book_keeping = (int *)malloc((n_nodes) * sizeof(int));
  std::memset(book_keeping, 0, (n_nodes) * sizeof(int));

  for (block_idx_t i = 0; i < n_blocks; i++) {
    enum ElemType element_type = element_types[i];
    int nnodesxelem = elem_num_nodes(element_type);

    for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
      for (ptrdiff_t j = 0; j < n_elements[i]; ++j) {
        const ptrdiff_t node = static_cast<ptrdiff_t>(elements[i][edof_i][j]);
        n2eptr[node + 1]++;
      }
    }
  }

  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    n2eptr[i + 1] += n2eptr[i];
  }

  element_idx_t *elindex =
      (element_idx_t *)malloc(n2eptr[n_nodes] * sizeof(element_idx_t));
  block_idx_t *block_number =
      (block_idx_t *)malloc(n2eptr[n_nodes] * sizeof(block_idx_t));

  for (block_idx_t i = 0; i < n_blocks; i++) {
    enum ElemType element_type = element_types[i];
    int nnodesxelem = elem_num_nodes(element_type);

    for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
      for (ptrdiff_t j = 0; j < n_elements[i]; ++j) {
        const ptrdiff_t node = static_cast<ptrdiff_t>(elements[i][edof_i][j]);

        SMESH_ASSERT(n2eptr[node] + book_keeping[node] < n2eptr[node + 1]);

        elindex[n2eptr[node] + book_keeping[node]] =
            static_cast<element_idx_t>(j);
        block_number[n2eptr[node] + book_keeping[node]++] = i;
      }
    }
  }

  free(book_keeping);

  *out_n2eptr = n2eptr;
  *out_elindex = elindex;
  *out_block_number = block_number;

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_crs_graph_from_n2e(
    const block_idx_t n_blocks, const enum ElemType element_types[],
    const ptrdiff_t n_elements[], const ptrdiff_t n_nodes,
    idx_t **const SMESH_RESTRICT elems[],
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex,
    const block_idx_t *const SMESH_RESTRICT block_number, count_t **out_rowptr,
    idx_t **out_colidx) {
  SMESH_UNUSED(n_blocks);
  SMESH_UNUSED(n_elements);
  count_t *rowptr = (count_t *)malloc((n_nodes + 1) * sizeof(count_t));
  idx_t *colidx = 0;

  {
    rowptr[0] = 0;

#pragma omp parallel
    {
      idx_t n2nbuff[4096];
#pragma omp for
      for (ptrdiff_t node = 0; node < n_nodes; ++node) {
        count_t ebegin = n2eptr[node];
        count_t eend = n2eptr[node + 1];

        count_t nneighs = 0;

        for (count_t e = ebegin; e < eend; ++e) {
          element_idx_t eidx = elindex[e];
          block_idx_t b = block_number[e];

          SMESH_ASSERT(b < n_blocks);
          SMESH_ASSERT(eidx < n_elements[b]);

          int nnodesxelem = elem_num_nodes(element_types[b]);

          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            idx_t neighnode = elems[b][edof_i][eidx];
            SMESH_ASSERT(nneighs < 4096);
            n2nbuff[nneighs++] = neighnode;
          }
        }

        nneighs = static_cast<count_t>(
            sort_and_unique(n2nbuff, static_cast<size_t>(nneighs)));
        rowptr[node + 1] = nneighs;
      }
    }

    // Cumulative sum
    for (ptrdiff_t node = 0; node < n_nodes; ++node) {
      rowptr[node + 1] += rowptr[node];
    }

    const ptrdiff_t nnz = rowptr[n_nodes];
    colidx = (idx_t *)malloc(nnz * sizeof(idx_t));

#pragma omp parallel
    {
      idx_t n2nbuff[4096];
#pragma omp for
      for (ptrdiff_t node = 0; node < n_nodes; ++node) {
        count_t ebegin = n2eptr[node];
        count_t eend = n2eptr[node + 1];

        count_t nneighs = 0;

        for (count_t e = ebegin; e < eend; ++e) {
          element_idx_t eidx = elindex[e];
          block_idx_t b = block_number[e];
          SMESH_ASSERT(eidx < n_elements[b]);

          int nnodesxelem = elem_num_nodes(element_types[b]);

          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            idx_t neighnode = elems[b][edof_i][eidx];
            SMESH_ASSERT(nneighs < 4096);
            n2nbuff[nneighs++] = neighnode;
          }
        }

        nneighs = static_cast<count_t>(
            sort_and_unique(n2nbuff, static_cast<size_t>(nneighs)));

        for (count_t i = 0; i < nneighs; ++i) {
          colidx[rowptr[node] + i] = n2nbuff[i];
        }
      }
    }
  }

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return 0;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_crs_graph(const block_idx_t n_blocks,
                               const enum ElemType element_types[],
                               const ptrdiff_t n_elements[],
                               idx_t **const SMESH_RESTRICT elems[],
                               const ptrdiff_t n_nodes,
                               count_t **out_rowptr, idx_t **out_colidx) {
  block_idx_t *block_number = 0;
  count_t *n2eptr = 0;
  element_idx_t *elindex = 0;

  create_multiblock_n2e(n_blocks, element_types, n_elements, elems, n_nodes,
                       &block_number, &n2eptr, &elindex);
  create_multiblock_crs_graph_from_n2e(n_blocks, element_types, n_elements,
                                      n_nodes, elems, n2eptr, elindex,
                                      block_number, out_rowptr, out_colidx);

  free(block_number);
  free(n2eptr);
  free(elindex);

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_crs_graph_upper_triangular_from_n2e(
    const block_idx_t n_blocks, const enum ElemType element_types[],
    const ptrdiff_t n_elements[], const ptrdiff_t n_nodes,
    idx_t **const SMESH_RESTRICT elems[],
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex,
    const block_idx_t *const SMESH_RESTRICT block_number, count_t **out_rowptr,
    idx_t **out_colidx) {
  SMESH_UNUSED(n_blocks);
  SMESH_UNUSED(n_elements);

  count_t *rowptr = (count_t *)malloc((n_nodes + 1) * sizeof(count_t));
  idx_t *colidx = 0;

  {
    rowptr[0] = 0;

    {
#pragma omp parallel for
      for (ptrdiff_t node = 0; node < n_nodes; ++node) {
        idx_t n2nbuff[4096];

        count_t ebegin = n2eptr[node];
        count_t eend = n2eptr[node + 1];

        count_t nneighs = 0;

        for (count_t e = ebegin; e < eend; ++e) {
          element_idx_t eidx = elindex[e];
          block_idx_t b = block_number[e];
          SMESH_ASSERT(b < n_blocks);
          SMESH_ASSERT(eidx < n_elements[b]);

          int nnodesxelem = elem_num_nodes(element_types[b]);

          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            idx_t neighnode = elems[b][edof_i][eidx];
            if (neighnode > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = neighnode;
            }
          }
        }

        nneighs = static_cast<count_t>(
            sort_and_unique(n2nbuff, static_cast<size_t>(nneighs)));
        rowptr[node + 1] = nneighs;
      }

      // Cumulative sum
      for (ptrdiff_t node = 0; node < n_nodes; ++node) {
        rowptr[node + 1] += rowptr[node];
      }

      const ptrdiff_t nnz = rowptr[n_nodes];
      colidx = (idx_t *)malloc(nnz * sizeof(idx_t));

      {
#pragma omp parallel for
        for (ptrdiff_t node = 0; node < n_nodes; ++node) {
          idx_t n2nbuff[4096];
          count_t ebegin = n2eptr[node];
          count_t eend = n2eptr[node + 1];

          count_t nneighs = 0;

          for (count_t e = ebegin; e < eend; ++e) {
            element_idx_t eidx = elindex[e];
            block_idx_t b = block_number[e];
            SMESH_ASSERT(eidx < n_elements[b]);

            int nnodesxelem = elem_num_nodes(element_types[b]);

            for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
              idx_t neighnode = elems[b][edof_i][eidx];
              if (neighnode > node) {
                SMESH_ASSERT(nneighs < 4096);
                n2nbuff[nneighs++] = neighnode;
              }
            }
          }

          nneighs = static_cast<count_t>(
              sort_and_unique(n2nbuff, static_cast<size_t>(nneighs)));

          for (count_t i = 0; i < nneighs; ++i) {
            colidx[rowptr[node] + i] = n2nbuff[i];
          }
        }
      }
    }
  }

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return 0;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_crs_graph_upper_triangular(
    const block_idx_t n_blocks, const enum ElemType element_types[],
    const ptrdiff_t n_elements[], idx_t **const SMESH_RESTRICT elems[],
    const ptrdiff_t n_nodes, count_t **out_rowptr, idx_t **out_colidx) {
  block_idx_t *block_number = 0;
  count_t *n2eptr = 0;
  element_idx_t *elindex = 0;

  create_multiblock_n2e(n_blocks, element_types, n_elements, elems, n_nodes,
                       &block_number, &n2eptr, &elindex);
  create_multiblock_crs_graph_upper_triangular_from_n2e(
      n_blocks, element_types, n_elements, n_nodes, elems, n2eptr, elindex,
      block_number, out_rowptr, out_colidx);

  free(block_number);
  free(n2eptr);
  free(elindex);

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_MULTIBLOCK_GRAPH_IMPL_HPP
