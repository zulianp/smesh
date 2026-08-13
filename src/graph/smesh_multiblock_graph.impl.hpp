#ifndef SMESH_MULTIBLOCK_GRAPH_IMPL_HPP
#define SMESH_MULTIBLOCK_GRAPH_IMPL_HPP

#include "smesh_multiblock_graph.hpp"
#include "smesh_alloc.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_graph.hpp"
#include "smesh_sort.hpp"
#include "smesh_ssquad4_graph.hpp"

#include <algorithm>
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
  count_t *n2eptr = (count_t *)SMESH_ALLOC((n_nodes + 1) * sizeof(count_t));
  std::memset(n2eptr, 0, (n_nodes + 1) * sizeof(count_t));

  int *book_keeping = (int *)SMESH_ALLOC((n_nodes) * sizeof(int));
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
      (element_idx_t *)SMESH_ALLOC(n2eptr[n_nodes] * sizeof(element_idx_t));
  block_idx_t *block_number = nullptr;
  if (write_block_number) {
    block_number = (block_idx_t *)SMESH_ALLOC(n2eptr[n_nodes] * sizeof(block_idx_t));
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

  SMESH_FREE(book_keeping);

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
  count_t *rowptr = (count_t *)SMESH_ALLOC((n_nodes + 1) * sizeof(count_t));
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
    colidx = (idx_t *)SMESH_ALLOC(nnz * sizeof(idx_t));

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

  SMESH_FREE(block_number);
  SMESH_FREE(n2eptr);
  SMESH_FREE(elindex);

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

  count_t *rowptr = (count_t *)SMESH_ALLOC((n_nodes + 1) * sizeof(count_t));
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
      colidx = (idx_t *)SMESH_ALLOC(nnz * sizeof(idx_t));

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

  SMESH_FREE(block_number);
  SMESH_FREE(n2eptr);
  SMESH_FREE(elindex);

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
      (element_idx_t *)SMESH_ALLOC((n_blocks + 1) * sizeof(element_idx_t));

  ptrdiff_t total_elements = 0;
  element_idx_t global_base = 0;
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    block_base[b] = global_base;
    global_base += static_cast<element_idx_t>(n_elements[b]);
    total_elements += n_elements[b];
  }
  block_base[n_blocks] = global_base;

  int *nnodesxelem = (int *)SMESH_ALLOC(n_blocks * sizeof(int));
  int *nnodesxside = (int *)SMESH_ALLOC(n_blocks * sizeof(int));
  int *nsides = (int *)SMESH_ALLOC(n_blocks * sizeof(int));

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

  int *connection_counter = (int *)SMESH_CALLOC(
      static_cast<size_t>(total_elements), sizeof(int));

  count_t *dual_e_ptr =
      (count_t *)SMESH_ALLOC(static_cast<size_t>(total_elements + 1) * sizeof(count_t));
  dual_e_ptr[0] = 0;

  element_idx_t *dual_eidx = (element_idx_t *)SMESH_ALLOC(
      static_cast<size_t>(n_overestimated_connections + extra_buffer_space) *
      sizeof(element_idx_t));
  block_idx_t *dual_eblock = (block_idx_t *)SMESH_ALLOC(
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
            SMESH_UNUSED(write_pos);

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

  SMESH_FREE(block_base);
  SMESH_FREE(nnodesxelem);
  SMESH_FREE(nnodesxside);
  SMESH_FREE(nsides);
  SMESH_FREE(connection_counter);

  *out_dual_eptr = dual_e_ptr;
  *out_dual_eidx = dual_eidx;
  *out_dual_eblock = dual_eblock;

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
static enum ElemType element_type_for_adjacency(enum ElemType element_type) {
  if (element_type == TET10) {
    return TET4;
  }
  if (element_type == TRI6) {
    return TRI3;
  }
  return element_type;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_half_face_table_for_block(
    const block_idx_t target_block, const block_idx_t n_blocks,
    const enum ElemType element_types[], const ptrdiff_t n_elements[],
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems[],
    const count_t *const SMESH_RESTRICT global_dual_ptr,
    const element_idx_t *const SMESH_RESTRICT global_dual_idx,
    const block_idx_t *const SMESH_RESTRICT global_dual_block,
    const element_idx_t *const SMESH_RESTRICT block_base,
    element_idx_t **out_table, block_idx_t **out_neighbor_block) {
  SMESH_ASSERT(target_block < n_blocks);

  enum ElemType element_type =
      element_type_for_adjacency<idx_t, count_t, element_idx_t>(
          element_types[target_block]);
  if (is_semistructured_type(element_type)) {
    SMESH_ERROR(
        "create_multiblock_half_face_table_for_block: semistructured type "
        "not supported: %s\n",
        type_to_string(element_types[target_block]));
    return SMESH_FAILURE;
  }

  LocalSideTable lst_src;
  lst_src.fill(element_type);

  const enum ElemType st = side_type(element_type);
  const int nn_src = elem_num_nodes(st);
  const int ns_src = elem_num_sides(element_type);
  const ptrdiff_t n_block_elements = n_elements[target_block];
  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems_src =
      elems[target_block];

  element_idx_t *table = (element_idx_t *)SMESH_ALLOC(
      static_cast<size_t>(n_block_elements * ns_src) * sizeof(element_idx_t));
  block_idx_t *neighbor_block = (block_idx_t *)SMESH_ALLOC(
      static_cast<size_t>(n_block_elements * ns_src) * sizeof(block_idx_t));

  std::vector<LocalSideTable> lst_blocks(n_blocks);
  std::vector<int> ns_blocks(n_blocks);
  std::vector<int> nn_blocks(n_blocks);
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    enum ElemType et =
        element_type_for_adjacency<idx_t, count_t, element_idx_t>(
            element_types[b]);
    lst_blocks[b].fill(et);
    ns_blocks[b] = elem_num_sides(et);
    nn_blocks[b] = elem_num_nodes(side_type(et));
  }

#pragma omp parallel
  {
    idx_t *nodes1 = (idx_t *)SMESH_ALLOC(
        LocalSideTable::MAX_NUM_NODES_PER_SIDE * sizeof(idx_t));
    idx_t *nodes2 = (idx_t *)SMESH_ALLOC(
        LocalSideTable::MAX_NUM_NODES_PER_SIDE * sizeof(idx_t));
    int *assigned =
        (int *)SMESH_ALLOC(LocalSideTable::MAX_NUM_SIDES * sizeof(int));

#pragma omp for
    for (ptrdiff_t e = 0; e < n_block_elements; ++e) {
      const element_idx_t g = block_base[target_block] + static_cast<element_idx_t>(e);
      const count_t begin = global_dual_ptr[g];
      const count_t end = global_dual_ptr[g + 1];
      const count_t range = end - begin;

      memset(assigned, 0, static_cast<size_t>(range) * sizeof(int));

      for (int s1 = 0; s1 < ns_src; ++s1) {
        table[static_cast<size_t>(e) * static_cast<size_t>(ns_src) +
              static_cast<size_t>(s1)] = invalid_idx<element_idx_t>();
        neighbor_block[static_cast<size_t>(e) * static_cast<size_t>(ns_src) +
                       static_cast<size_t>(s1)] = target_block;

        for (int j = 0; j < nn_src; ++j) {
          nodes1[j] = elems_src[lst_src(s1, j)][e];
        }
        std::sort(nodes1, nodes1 + nn_src);

        for (count_t k = 0; k < range; ++k) {
          if (assigned[k]) {
            continue;
          }

          const element_idx_t e_adj = global_dual_idx[begin + k];
          const block_idx_t b_adj = global_dual_block[begin + k];
          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems_adj =
              elems[b_adj];
          const LocalSideTable &lst_adj = lst_blocks[b_adj];
          const int ns_adj = ns_blocks[b_adj];
          const int nn_adj = nn_blocks[b_adj];

          for (int s2 = 0; s2 < ns_adj; ++s2) {
            for (int j = 0; j < nn_adj; ++j) {
              nodes2[j] = elems_adj[lst_adj(s2, j)][e_adj];
            }
            std::sort(nodes2, nodes2 + nn_adj);

            if (nn_adj != nn_src) {
              continue;
            }

            int diffs = 0;
            for (int j = 0; j < nn_src; ++j) {
              diffs += nodes1[j] != nodes2[j];
            }

            if (!diffs) {
              table[static_cast<size_t>(e) * static_cast<size_t>(ns_src) +
                    static_cast<size_t>(s1)] = e_adj;
              neighbor_block[static_cast<size_t>(e) * static_cast<size_t>(ns_src) +
                             static_cast<size_t>(s1)] = b_adj;
              assigned[k] = 1;
              break;
            }
          }
        }
      }
    }

    SMESH_FREE(nodes1);
    SMESH_FREE(nodes2);
    SMESH_FREE(assigned);
  }

  *out_table = table;
  *out_neighbor_block = neighbor_block;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_multiblock_edge_graph_from_n2e(
    const block_idx_t n_blocks, const enum ElemType element_types[],
    const ptrdiff_t n_elements[],
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems[],
    const ptrdiff_t n_nodes, const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex,
    const block_idx_t *const SMESH_RESTRICT block_number, count_t **out_rowptr,
    idx_t **out_colidx) {
  count_t *rowptr = (count_t *)SMESH_ALLOC((n_nodes + 1) * sizeof(count_t));
  rowptr[0] = 0;

#pragma omp parallel
  {
    idx_t n2nbuff[4096];
#pragma omp for
    for (ptrdiff_t node = 0; node < n_nodes; ++node) {
      count_t nneighs = 0;
      const count_t ebegin = n2eptr[node];
      const count_t eend = n2eptr[node + 1];

      for (count_t ei = ebegin; ei < eend; ++ei) {
        const block_idx_t b = block_number[ei];
        const element_idx_t eidx = elindex[ei];
        const enum ElemType element_type = element_types[b];
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems =
            elems[b];

        if (element_type == TET4 || element_type == TRI3 ||
            element_type == TRISHELL3) {
          const int nnodesxelem = elem_num_nodes(element_type);
          int lidx = -1;
          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
              lidx = edof_i;
              break;
            }
          }
          SMESH_ASSERT(lidx != -1);
          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            if (edof_i == lidx) {
              continue;
            }
            const idx_t neighnode = block_elems[edof_i][eidx];
            if (static_cast<ptrdiff_t>(neighnode) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = neighnode;
            }
          }
        } else if (element_type == QUAD4 || element_type == QUADSHELL4) {
          static const int nnodesxelem = 4;
          int lidx = -1;
          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
              lidx = edof_i;
              break;
            }
          }
          SMESH_ASSERT(lidx != -1);
          static int quad4_edge_connectivity[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
          for (int d = 0; d < 2; ++d) {
            const idx_t neighnode =
                block_elems[quad4_edge_connectivity[lidx][d]][eidx];
            if (static_cast<ptrdiff_t>(neighnode) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = neighnode;
            }
          }
        } else if (element_type == HEX8) {
          static int hex8_edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                          {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
          for (int e = 0; e < 12; ++e) {
            const idx_t n0 = block_elems[hex8_edges[e][0]][eidx];
            const idx_t n1 = block_elems[hex8_edges[e][1]][eidx];
            if (n0 == static_cast<idx_t>(node) && static_cast<ptrdiff_t>(n1) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = n1;
            }
            if (n1 == static_cast<idx_t>(node) && static_cast<ptrdiff_t>(n0) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = n0;
            }
          }
        } else {
          SMESH_ERROR(
              "create_multiblock_edge_graph_from_n2e: unsupported element "
              "type %s\n",
              type_to_string(element_type));
        }
      }

      nneighs = static_cast<count_t>(
          sort_and_unique(n2nbuff, static_cast<size_t>(nneighs)));
      rowptr[node + 1] = nneighs;
    }
  }

  for (ptrdiff_t node = 0; node < n_nodes; ++node) {
    rowptr[node + 1] += rowptr[node];
  }

  const ptrdiff_t nnz = rowptr[n_nodes];
  idx_t *colidx = (idx_t *)SMESH_ALLOC(static_cast<size_t>(nnz) * sizeof(idx_t));

#pragma omp parallel
  {
    idx_t n2nbuff[4096];
#pragma omp for
    for (ptrdiff_t node = 0; node < n_nodes; ++node) {
      count_t nneighs = 0;
      const count_t ebegin = n2eptr[node];
      const count_t eend = n2eptr[node + 1];

      for (count_t ei = ebegin; ei < eend; ++ei) {
        const block_idx_t b = block_number[ei];
        const element_idx_t eidx = elindex[ei];
        const enum ElemType element_type = element_types[b];
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT block_elems =
            elems[b];

        if (element_type == TET4 || element_type == TRI3 ||
            element_type == TRISHELL3) {
          const int nnodesxelem = elem_num_nodes(element_type);
          int lidx = -1;
          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
              lidx = edof_i;
              break;
            }
          }
          SMESH_ASSERT(lidx != -1);
          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            if (edof_i == lidx) {
              continue;
            }
            const idx_t neighnode = block_elems[edof_i][eidx];
            if (static_cast<ptrdiff_t>(neighnode) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = neighnode;
            }
          }
        } else if (element_type == QUAD4 || element_type == QUADSHELL4) {
          static const int nnodesxelem = 4;
          int lidx = -1;
          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            if (block_elems[edof_i][eidx] == static_cast<idx_t>(node)) {
              lidx = edof_i;
              break;
            }
          }
          SMESH_ASSERT(lidx != -1);
          static int quad4_edge_connectivity[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
          for (int d = 0; d < 2; ++d) {
            const idx_t neighnode =
                block_elems[quad4_edge_connectivity[lidx][d]][eidx];
            if (static_cast<ptrdiff_t>(neighnode) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = neighnode;
            }
          }
        } else if (element_type == HEX8) {
          static int hex8_edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                          {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
          for (int e = 0; e < 12; ++e) {
            const idx_t n0 = block_elems[hex8_edges[e][0]][eidx];
            const idx_t n1 = block_elems[hex8_edges[e][1]][eidx];
            if (n0 == static_cast<idx_t>(node) && static_cast<ptrdiff_t>(n1) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = n1;
            }
            if (n1 == static_cast<idx_t>(node) && static_cast<ptrdiff_t>(n0) > node) {
              SMESH_ASSERT(nneighs < 4096);
              n2nbuff[nneighs++] = n0;
            }
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

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_MULTIBLOCK_GRAPH_IMPL_HPP
