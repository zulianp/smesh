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
                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements[],
                         const ptrdiff_t n_nodes,
                         block_idx_t **out_block_number, count_t **out_n2eptr,
                         element_idx_t **out_elindex) {
  count_t *n2eptr = (count_t *)malloc((n_nodes + 1) * sizeof(count_t));
  std::memset(n2eptr, 0, (n_nodes + 1) * sizeof(count_t));

  int *book_keeping = (int *)malloc((n_nodes) * sizeof(int));
  std::memset(book_keeping, 0, (n_nodes) * sizeof(int));

  const bool write_block_number = (out_block_number != nullptr);

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
  block_idx_t *block_number = nullptr;
  if (write_block_number) {
    block_number = (block_idx_t *)malloc(n2eptr[n_nodes] * sizeof(block_idx_t));
  }

  element_idx_t global_element_base = 0;
  for (block_idx_t i = 0; i < n_blocks; i++) {
    enum ElemType element_type = element_types[i];
    int nnodesxelem = elem_num_nodes(element_type);

    for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
      for (ptrdiff_t j = 0; j < n_elements[i]; ++j) {
        const ptrdiff_t node = static_cast<ptrdiff_t>(elements[i][edof_i][j]);

        SMESH_ASSERT(n2eptr[node] + book_keeping[node] < n2eptr[node + 1]);

        const count_t pos =
            static_cast<count_t>(n2eptr[node] + book_keeping[node]);
        if (write_block_number) {
          elindex[pos] = static_cast<element_idx_t>(j);
          block_number[pos] = i;
        } else {
          elindex[pos] = global_element_base + static_cast<element_idx_t>(j);
        }

        book_keeping[node]++;
      }
    }

    global_element_base += static_cast<element_idx_t>(n_elements[i]);
  }

  free(book_keeping);

  *out_n2eptr = n2eptr;
  *out_elindex = elindex;
  if (write_block_number) {
    *out_block_number = block_number;
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_crs_graph_from_n2e(
    const block_idx_t n_blocks, const enum ElemType element_types[],
    const ptrdiff_t n_elements[], const ptrdiff_t n_nodes,
    const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT elems[],
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
                               const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT elems[],
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
    const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT elems[],
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
    const ptrdiff_t n_elements[], const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT elems[],
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

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_dual_graph_from_n2e(
    const block_idx_t n_blocks, const enum ElemType element_types[],
    const ptrdiff_t n_elements[], const ptrdiff_t n_nodes,
    const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT elems[],
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex,
    const block_idx_t *const SMESH_RESTRICT block_number,
    count_t **out_dual_eptr, element_idx_t **out_dual_eidx,
    block_idx_t **out_dual_eblock) {
  SMESH_UNUSED(n_nodes);

  element_idx_t *block_base =
      (element_idx_t *)malloc((n_blocks + 1) * sizeof(element_idx_t));

  ptrdiff_t total_elements = 0;
  element_idx_t global_base = 0;
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    block_base[b] = global_base;
    global_base += static_cast<element_idx_t>(n_elements[b]);
    total_elements += n_elements[b];
  }
  block_base[n_blocks] = global_base;

  int *nnodesxelem = (int *)malloc(n_blocks * sizeof(int));
  int *nnodesxside = (int *)malloc(n_blocks * sizeof(int));
  int *nsides = (int *)malloc(n_blocks * sizeof(int));

  ptrdiff_t n_overestimated_connections = 0;
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    enum ElemType element_type_for_algo = element_types[b];
    if (element_type_for_algo == TET10) {
      element_type_for_algo = TET4;
    } else if (element_type_for_algo == TRI6) {
      element_type_for_algo = TRI3;
    }

    nsides[b] = elem_num_sides(element_type_for_algo);
    nnodesxelem[b] = elem_num_nodes(element_type_for_algo);
    nnodesxside[b] = elem_num_nodes(side_type(element_type_for_algo));

    n_overestimated_connections += n_elements[b] * nsides[b];
  }

  const ptrdiff_t extra_buffer_space = 1000;

  int *connection_counter = (int *)calloc(
      static_cast<size_t>(total_elements), sizeof(int));

  count_t *dual_e_ptr =
      (count_t *)malloc(static_cast<size_t>(total_elements + 1) * sizeof(count_t));
  dual_e_ptr[0] = 0;

  element_idx_t *dual_eidx = (element_idx_t *)malloc(
      static_cast<size_t>(n_overestimated_connections + extra_buffer_space) *
      sizeof(element_idx_t));
  block_idx_t *dual_eblock = (block_idx_t *)malloc(
      static_cast<size_t>(n_overestimated_connections + extra_buffer_space) *
      sizeof(block_idx_t));

  element_idx_t g = 0;
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    for (ptrdiff_t e = 0; e < n_elements[b]; ++e, ++g) {
      const count_t offset = dual_e_ptr[g];
      element_idx_t *elist = &dual_eidx[offset];
      block_idx_t *blist = &dual_eblock[offset];

      int count_common = 0;

      for (int en = 0; en < nnodesxelem[b]; ++en) {
        const idx_t node = elems[b][en][e];

        for (count_t eii = n2eptr[node]; eii < n2eptr[node + 1]; ++eii) {
          const block_idx_t b_adj = block_number[eii];
          const element_idx_t e_adj = elindex[eii];

          SMESH_ASSERT(b_adj < n_blocks);
          SMESH_ASSERT(e_adj < n_elements[b_adj]);

          const element_idx_t g_adj = block_base[b_adj] + e_adj;
          const ptrdiff_t g_adj_p = static_cast<ptrdiff_t>(g_adj);
          SMESH_ASSERT(g_adj_p < total_elements);

          if (connection_counter[g_adj_p] == 0) {
            const ptrdiff_t write_pos =
                static_cast<ptrdiff_t>(offset) + count_common;
            SMESH_ASSERT(write_pos <
                         n_overestimated_connections + extra_buffer_space);

            elist[count_common] = g_adj;
            blist[count_common] = b_adj;
            ++count_common;
          }

          connection_counter[g_adj_p]++;
        }
      }

      connection_counter[static_cast<ptrdiff_t>(g)] = 0;

      const int required_src = nnodesxside[b];
      int actual_count = 0;

      for (int ec = 0; ec < count_common; ++ec) {
        const element_idx_t g_adj = elist[ec];
        const ptrdiff_t g_adj_p = static_cast<ptrdiff_t>(g_adj);
        const int overlap = connection_counter[g_adj_p];
        const block_idx_t b_adj = blist[ec];

        if (overlap == required_src && overlap == nnodesxside[b_adj]) {
          elist[actual_count] = g_adj - block_base[b_adj];
          blist[actual_count] = b_adj;
          ++actual_count;
        }

        connection_counter[g_adj_p] = 0;
      }

      dual_e_ptr[g + 1] = offset + actual_count;
    }
  }

  free(block_base);
  free(nnodesxelem);
  free(nnodesxside);
  free(nsides);
  free(connection_counter);

  *out_dual_eptr = dual_e_ptr;
  *out_dual_eidx = dual_eidx;
  *out_dual_eblock = dual_eblock;

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_MULTIBLOCK_GRAPH_IMPL_HPP
