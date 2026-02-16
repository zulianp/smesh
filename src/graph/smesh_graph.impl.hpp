#ifndef SMESH_GRAPH_IMPL_HPP
#define SMESH_GRAPH_IMPL_HPP

#include "smesh_graph.hpp"
#include "smesh_tracer.hpp"

#include "smesh_base.hpp"
#include "smesh_sort.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include "matrixio_array.h"
// #include "matrixio_crs.h"
// #include "utils.h"

// #include "sort_and_unique.h"
// #include "bitonic.h"

// https://dirtyhandscoding.github.io/posts/vectorizing-small-fixed-size-sort.html
// https://xhad1234.github.io/Parallel-Sort-Merge-Join-in-Peloton/
// https://github.com/sid1607/avx2-merge-sort/blob/master/merge_sort.h
// https://onlinelibrary.wiley.com/doi/full/10.1002/spe.2922

namespace smesh {

template <typename T> SMESH_INLINE bool ispow2(T n) {
  return n && (!(n & (n - 1)));
}

template <typename idx_t>
idx_t find_idx(const idx_t target, const idx_t *x, idx_t n) {
  for (idx_t i = 0; i < n; ++i) {
    if (target == x[i]) {
      return i;
    }
  }

  return n;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(const ptrdiff_t nelements, const ptrdiff_t nnodes,
               const int nnodesxelem,
               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
               count_t **out_n2eptr, element_idx_t **out_elindex) {
  SMESH_TRACE_SCOPE("create_n2e");

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("build_n2e: allocating %g GB\n",
         (nnodes + 1) * sizeof(count_t) * 1e-9);
#endif

  count_t *n2eptr = (count_t *)malloc((nnodes + 1) * sizeof(count_t));
  memset(n2eptr, 0, (nnodes + 1) * sizeof(count_t));

  int *book_keeping = (int *)malloc((nnodes) * sizeof(int));
  memset(book_keeping, 0, (nnodes) * sizeof(int));

  for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
    for (ptrdiff_t i = 0; i < nelements; ++i) {
      assert(elems[edof_i][i] < nnodes);
      assert(elems[edof_i][i] >= 0);

      ++n2eptr[elems[edof_i][i] + 1];
    }
  }

  for (ptrdiff_t i = 0; i < nnodes; ++i) {
    n2eptr[i + 1] += n2eptr[i];
  }

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("build_n2e: allocating %g GB\n",
         n2eptr[nnodes] * sizeof(element_idx_t) * 1e-9);
#endif
  element_idx_t *elindex = (element_idx_t *)malloc(
      static_cast<size_t>(n2eptr[nnodes]) * sizeof(element_idx_t));

  for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
    for (ptrdiff_t i = 0; i < nelements; ++i) {
      const idx_t node = elems[edof_i][i];

      assert(n2eptr[node] + book_keeping[node] < n2eptr[node + 1]);

      elindex[n2eptr[node] + book_keeping[node]++] =
          static_cast<element_idx_t>(i);
    }
  }

  free(book_keeping);

  *out_n2eptr = n2eptr;
  *out_elindex = elindex;

  return SMESH_SUCCESS;
}

template <typename count_t, typename element_idx_t>
int sort_n2e(const ptrdiff_t nnodes, const count_t *const SMESH_RESTRICT n2eptr,
             element_idx_t *const SMESH_RESTRICT elindex) {
  SMESH_TRACE_SCOPE("sort_n2e");

#pragma omp parallel for
  for (ptrdiff_t node = 0; node < nnodes; ++node) {
    const count_t ebegin = n2eptr[node];
    const count_t eend = n2eptr[node + 1];

    std::sort(elindex + ebegin, elindex + eend);
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
static int create_n2e_for_elem_type(
    const enum ElemType element_type, const ptrdiff_t nelements,
    const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_n2eptr, element_idx_t **out_elindex) {
  SMESH_TRACE_SCOPE("create_n2e_for_elem_type");
  // TODO (maybe)
  if (element_type != MACRO_TET4 /*&& element_type != MACRO_TRI3*/) {
    return create_n2e<idx_t, count_t, element_idx_t>(
        nelements, nnodes, elem_num_nodes(element_type), elems, out_n2eptr,
        out_elindex);
  }

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("build_n2e_for_elem_type: allocating %g GB\n",
         (nnodes + 1) * sizeof(count_t) * 1e-9);
#endif

  count_t *n2eptr = (count_t *)malloc((nnodes + 1) * sizeof(count_t));
  memset(n2eptr, 0, (nnodes + 1) * sizeof(count_t));

  int *book_keeping = (int *)malloc((nnodes) * sizeof(int));
  memset(book_keeping, 0, (nnodes) * sizeof(int));

  if (element_type == MACRO_TET4) {
    static const int tet4_refine_pattern[8][4] = {// Corner tests
                                                  {0, 4, 6, 7},
                                                  {4, 1, 5, 8},
                                                  {6, 5, 2, 9},
                                                  {7, 8, 9, 3},
                                                  // Octahedron tets
                                                  {4, 5, 6, 8},
                                                  {7, 4, 6, 8},
                                                  {6, 5, 9, 8},
                                                  {7, 6, 9, 8}};

    for (int sub_elem = 0; sub_elem < 8; sub_elem++) {
      for (int sub_elem_node = 0; sub_elem_node < 4; ++sub_elem_node) {
        int node_number = tet4_refine_pattern[sub_elem][sub_elem_node];

        for (ptrdiff_t i = 0; i < nelements; ++i) {
          assert(elems[node_number][i] < nnodes);
          assert(elems[node_number][i] >= 0);

          ++n2eptr[elems[node_number][i] + 1];
        }
      }
    }
  }

  for (ptrdiff_t i = 0; i < nnodes; ++i) {
    n2eptr[i + 1] += n2eptr[i];
  }

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("build_n2e: allocating %g GB\n",
         n2eptr[nnodes] * sizeof(element_idx_t) * 1e-9);
#endif
  element_idx_t *elindex = (element_idx_t *)malloc(
      static_cast<size_t>(n2eptr[nnodes]) * sizeof(element_idx_t));

  const int nnodesxelem = elem_num_nodes(element_type);
  for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
    for (ptrdiff_t i = 0; i < nelements; ++i) {
      const idx_t node = elems[edof_i][i];

      assert(n2eptr[node] + book_keeping[node] < n2eptr[node + 1]);

      elindex[n2eptr[node] + book_keeping[node]++] =
          static_cast<element_idx_t>(i);
    }
  }

  free(book_keeping);

  *out_n2eptr = n2eptr;
  *out_elindex = elindex;

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int n2n_from_n2e(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                 const int nnodesxelem,
                 const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                 const count_t *const SMESH_RESTRICT n2eptr,
                 const element_idx_t *const SMESH_RESTRICT elindex,
                 count_t **out_rowptr, idx_t **out_colidx) {
  SMESH_UNUSED(nelements);
  count_t *rowptr = (count_t *)malloc((nnodes + 1) * sizeof(count_t));
  idx_t *colidx = 0;

  {
    rowptr[0] = 0;

#pragma omp parallel
    {
      idx_t n2nbuff[4096];
#pragma omp for
      for (ptrdiff_t node = 0; node < nnodes; ++node) {
        count_t ebegin = n2eptr[node];
        count_t eend = n2eptr[node + 1];

        idx_t nneighs = 0;

        for (count_t e = ebegin; e < eend; ++e) {
          element_idx_t eidx = elindex[e];
          assert(eidx < nelements);

          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            idx_t neighnode = elems[edof_i][eidx];
            assert(nneighs < 4096);
            n2nbuff[nneighs++] = neighnode;
          }
        }

        nneighs = sort_and_unique(n2nbuff, nneighs);
        rowptr[node + 1] = nneighs;
      }
    }

    // Cumulative sum
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
      rowptr[node + 1] += rowptr[node];
    }

    const ptrdiff_t nnz = rowptr[nnodes];
    colidx = (idx_t *)malloc(nnz * sizeof(idx_t));

#pragma omp parallel
    {
      idx_t n2nbuff[4096];
#pragma omp for
      for (ptrdiff_t node = 0; node < nnodes; ++node) {
        count_t ebegin = n2eptr[node];
        count_t eend = n2eptr[node + 1];

        idx_t nneighs = 0;

        for (count_t e = ebegin; e < eend; ++e) {
          element_idx_t eidx = elindex[e];
          assert(eidx < nelements);

          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            idx_t neighnode = elems[edof_i][eidx];
            assert(nneighs < 4096);
            n2nbuff[nneighs++] = neighnode;
          }
        }

        nneighs = sort_and_unique(n2nbuff, nneighs);

        for (idx_t i = 0; i < nneighs; ++i) {
          colidx[rowptr[node] + i] = n2nbuff[i];
        }
      }
    }
  }

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
static int create_crs_graph_mem_conservative(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx) {

  SMESH_TRACE_SCOPE("create_crs_graph_mem_conservative");
  using element_idx_t_local = ptrdiff_t;

  count_t *n2eptr;
  element_idx_t_local *elindex;
  create_n2e<idx_t, count_t, element_idx_t_local>(
      nelements, nnodes, nnodesxelem, elems, &n2eptr, &elindex);

  int err = n2n_from_n2e<idx_t, count_t, element_idx_t_local>(
      nelements, nnodes, nnodesxelem, elems, n2eptr, elindex, out_rowptr,
      out_colidx);

  free(n2eptr);
  free(elindex);

  return err;
}

template <typename idx_t, typename count_t>
static int create_crs_graph_faster(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx) {
  SMESH_TRACE_SCOPE("create_crs_graph_faster");
  using element_idx_t_local = ptrdiff_t;

  ptrdiff_t nnz = 0;
  count_t *rowptr = (count_t *)malloc((nnodes + 1) * sizeof(count_t));
  idx_t *colidx = 0;

  {
    count_t *n2eptr;
    element_idx_t_local *elindex;
    create_n2e<idx_t, count_t, element_idx_t_local>(
        nelements, nnodes, nnodesxelem, elems, &n2eptr, &elindex);

    rowptr[0] = 0;

    ptrdiff_t overestimated_nnz = 0;
#pragma omp parallel for reduction(+ : overestimated_nnz)
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
      const count_t ebegin = n2eptr[node];
      const count_t eend = n2eptr[node + 1];
      idx_t nneighs = (eend - ebegin) * nnodesxelem;
      overestimated_nnz += nneighs;
    }

    colidx = (idx_t *)malloc(overestimated_nnz * sizeof(idx_t));

    ptrdiff_t coloffset = 0;
    idx_t n2nbuff[2048];
    for (ptrdiff_t node = 0; node < nnodes; ++node) {
      const count_t ebegin = n2eptr[node];
      const count_t eend = n2eptr[node + 1];

      idx_t nneighs = 0;
      for (count_t e = ebegin; e < eend; ++e) {
        const element_idx_t_local eidx = elindex[e];
        assert(eidx < nelements);

        for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
          const idx_t neighnode = elems[edof_i][eidx];
          assert(nneighs < 2048);
          n2nbuff[nneighs++] = neighnode;
        }
      }

      nneighs = sort_and_unique(n2nbuff, nneighs);

      nnz += nneighs;
      rowptr[node + 1] = nnz;

      for (idx_t i = 0; i < nneighs; ++i) {
        colidx[coloffset + i] = n2nbuff[i];
      }

      coloffset += nneighs;
    }

    free(n2eptr);
    free(elindex);
  }

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
int create_crs_graph_for_elem_type(
    const enum ElemType element_type, const ptrdiff_t nelements,
    const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx) {
  int SMESH_CRS_FAST_SERIAL = 0;
  SMESH_READ_ENV(SMESH_CRS_FAST_SERIAL, atoi);

  if (SMESH_CRS_FAST_SERIAL) {
    return create_crs_graph_faster<idx_t, count_t>(
        nelements, nnodes, elem_num_nodes(element_type), elems, out_rowptr,
        out_colidx);
  }

  return create_crs_graph_mem_conservative<idx_t, count_t>(
      nelements, nnodes, elem_num_nodes(element_type), elems, out_rowptr,
      out_colidx);
}

template <typename idx_t, typename count_t>
int create_crs_graph_from_element(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, int nxe,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx) {
  int SMESH_CRS_FAST_SERIAL = 0;
  SMESH_READ_ENV(SMESH_CRS_FAST_SERIAL, atoi);

  if (SMESH_CRS_FAST_SERIAL) {
    return create_crs_graph_faster<idx_t, count_t>(
        nelements, nnodes, nxe, elems, out_rowptr, out_colidx);
  }

  return create_crs_graph_mem_conservative<idx_t, count_t>(
      nelements, nnodes, nxe, elems, out_rowptr, out_colidx);
}

template <typename idx_t, typename count_t>
int create_crs_graph(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                         elems,
                     count_t **out_rowptr, idx_t **out_colidx) {
  return create_crs_graph_for_elem_type<idx_t, count_t>(
      TET4, nelements, nnodes, elems, out_rowptr, out_colidx);
}

template <typename idx_t, typename count_t>
int create_crs_graph_3(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                           elems,
                       count_t **out_rowptr, idx_t **out_colidx) {
  return create_crs_graph_for_elem_type<idx_t, count_t>(
      TRI3, nelements, nnodes, elems, out_rowptr, out_colidx);
}

template <typename idx_t, typename count_t, typename real_t>
int block_crs_to_crs(const ptrdiff_t nnodes, const int block_size,
                     const count_t *const SMESH_RESTRICT block_rowptr,
                     const idx_t *const SMESH_RESTRICT block_colidx,
                     const real_t *const SMESH_RESTRICT block_values,
                     count_t *const SMESH_RESTRICT rowptr,
                     idx_t *const SMESH_RESTRICT colidx,
                     real_t *const SMESH_RESTRICT values) {
  for (ptrdiff_t i = 0; i < nnodes; ++i) {
    count_t k = block_rowptr[i] * (block_size * block_size);
    count_t ncols = block_rowptr[i + 1] - block_rowptr[i];

    for (int b = 0; b < block_size; ++b) {
      rowptr[i * block_size + b] = k + ncols * (b * block_size);
    }
  }

  rowptr[nnodes * block_size] =
      2 * rowptr[nnodes * block_size - 1] - rowptr[nnodes * block_size - 2];

  for (ptrdiff_t i = 0; i < nnodes; ++i) {
    // Block row
    const count_t bstart = block_rowptr[i];
    const count_t bend = block_rowptr[i + 1];

    for (int brow = 0; brow < block_size; ++brow) {
      const idx_t row = i * block_size + brow;
      // Scalar row
      const count_t start = rowptr[row];
#ifndef NDEBUG
      const count_t end = rowptr[row + 1];
#endif

      for (count_t bk = bstart, k = start; bk < bend; ++bk) {
        // Block column
        const idx_t bcolidx = block_colidx[bk];
        // Data block
        const real_t *block = &block_values[bk * block_size * block_size];

        for (int bcol = 0; bcol < block_size; ++bcol, ++k) {
          assert(k < end);

          colidx[k] = bcolidx * block_size + bcol;
          values[k] = block[bcol * block_size + brow];
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
int crs_graph_block_to_scalar(const ptrdiff_t nnodes, const int block_size,
                              const count_t *const SMESH_RESTRICT block_rowptr,
                              const idx_t *const SMESH_RESTRICT block_colidx,
                              count_t *const SMESH_RESTRICT rowptr,
                              idx_t *const SMESH_RESTRICT colidx) {
  for (ptrdiff_t i = 0; i < nnodes; ++i) {
    count_t k = block_rowptr[i] * (block_size * block_size);
    count_t ncols = block_rowptr[i + 1] - block_rowptr[i];

    for (int b = 0; b < block_size; ++b) {
      rowptr[i * block_size + b] = k + ncols * (b * block_size);
    }
  }

  rowptr[nnodes * block_size] =
      2 * rowptr[nnodes * block_size - 1] - rowptr[nnodes * block_size - 2];

  for (ptrdiff_t i = 0; i < nnodes; ++i) {
    // Block row
    const count_t bstart = block_rowptr[i];
    const count_t bend = block_rowptr[i + 1];

    for (int brow = 0; brow < block_size; ++brow) {
      const idx_t row = i * block_size + brow;
      // Scalar row
      const count_t start = rowptr[row];
#ifndef NDEBUG
      const count_t end = rowptr[row + 1];
#endif

      for (count_t bk = bstart, k = start; bk < bend; ++bk) {
        // Block column
        const idx_t bcolidx = block_colidx[bk];
        // Data block
        for (int bcol = 0; bcol < block_size; ++bcol, ++k) {
          assert(k < end);

          colidx[k] = bcolidx * block_size + bcol;
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
static int create_dual_graph_mem_conservative(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_dual_eptr, element_idx_t **out_dual_eidx) {
  count_t *n2eptr = 0;
  element_idx_t *elindex = 0;

  if (element_type == TET10) {
    create_n2e<idx_t, count_t, element_idx_t>(
        n_elements, n_nodes, elem_num_nodes(TET4), elems, &n2eptr, &elindex);
  } else {
    create_n2e<idx_t, count_t, element_idx_t>(n_elements, n_nodes,
                                              elem_num_nodes(element_type),
                                              elems, &n2eptr, &elindex);
  }

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("create_dual_graph_mem_conservative: allocating %g GB\n",
         n_elements * sizeof(int) * 1e-9);
#endif

  int *connection_counter = (int *)malloc(n_elements * sizeof(int));
  memset(connection_counter, 0, n_elements * sizeof(int));

  const int n_sides = elem_num_sides(element_type);
  int n_nodes_per_elem = elem_num_nodes(element_type);

  // Optimize for Tet10
  if (element_type == TET10) {
    n_nodes_per_elem = 4;
  }

  enum ElemType st = side_type(element_type);
  int n_nodes_per_side = elem_num_nodes(st);

  if (element_type == TET10) {
    n_nodes_per_side = 3;
  }

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("create_dual_graph_mem_conservative: allocating %g GB\n",
         (n_elements + 1) * sizeof(count_t) * 1e-9);
#endif
  count_t *dual_e_ptr = (count_t *)calloc((n_elements + 1), sizeof(count_t));

  const ptrdiff_t n_overestimated_connections = n_elements * n_sides;
  // +1 more to avoid illegal access when counting self
  const ptrdiff_t extra_buffer_space = 1000;

#ifdef SMESH_ENABLE_MEM_DIAGNOSTICS
  printf("create_dual_graph_mem_conservative: allocating %g GB\n",
         (n_overestimated_connections + extra_buffer_space) *
             sizeof(element_idx_t) * 1e-9);
#endif

  element_idx_t *dual_eidx = (element_idx_t *)calloc(
      static_cast<size_t>(n_overestimated_connections + extra_buffer_space),
      sizeof(element_idx_t));

  for (ptrdiff_t e = 0; e < n_elements; e++) {
    count_t offset = dual_e_ptr[e];
    element_idx_t *elist = &dual_eidx[offset];

    int count_common = 0;
    for (int en = 0; en < n_nodes_per_elem; en++) {
      const idx_t node = elems[en][e];

      for (count_t eii = n2eptr[node]; eii < n2eptr[node + 1]; eii++) {
        const element_idx_t e_adj = elindex[eii];
        assert(e_adj < n_elements);

        if (connection_counter[e_adj] == 0) {
          const ptrdiff_t write_pos = static_cast<ptrdiff_t>(offset) +
                                      static_cast<ptrdiff_t>(count_common);
          assert(write_pos < n_overestimated_connections + extra_buffer_space);
          SMESH_UNUSED(write_pos);

          elist[count_common++] = e_adj;
        }

        connection_counter[e_adj]++;
      }
    }

    connection_counter[e] = 0;

    int actual_count = 0;
    for (int ec = 0; ec < count_common; ec++) {
      element_idx_t l = elist[ec];
      int overlap = connection_counter[l];
      assert(overlap <= n_nodes_per_elem);

      if (overlap == n_nodes_per_side) {
        elist[actual_count++] = l;
      }

      connection_counter[l] = 0;
    }

    dual_e_ptr[e + 1] = actual_count + offset;
  }

  free(n2eptr);
  free(elindex);
  free(connection_counter);

  *out_dual_eptr = dual_e_ptr;
  *out_dual_eidx = dual_eidx;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_dual_graph(const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
                      const enum ElemType element_type,
                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                          elems,
                      count_t **out_rowptr, element_idx_t **out_colidx) {
  SMESH_TRACE_SCOPE("create_dual_graph");
  const int ret =
      create_dual_graph_mem_conservative<idx_t, count_t, element_idx_t>(
          n_elements, n_nodes, element_type, elems, out_rowptr, out_colidx);

  return ret;
}

template <typename idx_t, typename count_t, typename element_idx_t>
static int create_crs_graph_upper_triangular_from_n2e(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex, count_t **out_rowptr,
    idx_t **out_colidx) {
  SMESH_UNUSED(nelements);
  count_t *rowptr = (count_t *)malloc((nnodes + 1) * sizeof(count_t));
  idx_t *colidx = 0;

  {
    rowptr[0] = 0;

    {
#pragma omp parallel for
      for (ptrdiff_t node = 0; node < nnodes; ++node) {
        idx_t n2nbuff[4096];

        count_t ebegin = n2eptr[node];
        count_t eend = n2eptr[node + 1];

        idx_t nneighs = 0;

        for (count_t e = ebegin; e < eend; ++e) {
          element_idx_t eidx = elindex[e];
          assert(eidx < nelements);

          for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            idx_t neighnode = elems[edof_i][eidx];
            if (neighnode > node) {
              assert(nneighs < 4096);
              n2nbuff[nneighs++] = neighnode;
            }
          }

          nneighs = sort_and_unique(n2nbuff, nneighs);
          rowptr[node + 1] = nneighs;
        }
      }

      // Cumulative sum
      for (ptrdiff_t node = 0; node < nnodes; ++node) {
        rowptr[node + 1] += rowptr[node];
      }

      const ptrdiff_t nnz = rowptr[nnodes];
      colidx = (idx_t *)malloc(nnz * sizeof(idx_t));

      {
#pragma omp parallel for
        for (ptrdiff_t node = 0; node < nnodes; ++node) {
          idx_t n2nbuff[4096];
          count_t ebegin = n2eptr[node];
          count_t eend = n2eptr[node + 1];

          idx_t nneighs = 0;

          for (count_t e = ebegin; e < eend; ++e) {
            element_idx_t eidx = elindex[e];
            assert(eidx < nelements);

            for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
              idx_t neighnode = elems[edof_i][eidx];
              if (neighnode > node) {
                assert(nneighs < 4096);
                n2nbuff[nneighs++] = neighnode;
              }
            }

            nneighs = sort_and_unique(n2nbuff, nneighs);

            for (idx_t i = 0; i < nneighs; ++i) {
              colidx[rowptr[node] + i] = n2nbuff[i];
            }
          }
        }
      }
    }
  }

  *out_rowptr = rowptr;
  *out_colidx = colidx;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
int create_crs_graph_upper_triangular_from_element(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, int nxe,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx) {
  SMESH_TRACE_SCOPE("create_crs_graph_upper_triangular_from_element");
  using element_idx_t_local = ptrdiff_t;

  count_t *n2eptr;
  element_idx_t_local *elindex;
  create_n2e<idx_t, count_t, element_idx_t_local>(nelements, nnodes, nxe, elems,
                                                  &n2eptr, &elindex);

  int err = create_crs_graph_upper_triangular_from_n2e<idx_t, count_t,
                                                       element_idx_t_local>(
      nelements, nnodes, nxe, elems, n2eptr, elindex, out_rowptr, out_colidx);

  free(n2eptr);
  free(elindex);

  return err;
}

template <typename idx_t, typename count_t>
int crs_to_coo(const ptrdiff_t n, const count_t *const SMESH_RESTRICT rowptr,
               idx_t *const SMESH_RESTRICT row_idx) {
#pragma omp parallel for
  for (ptrdiff_t row = 0; row < n; row++) {
    for (count_t k = rowptr[row]; k < rowptr[row + 1]; k++) {
      row_idx[k] = row;
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
int sorted_coo_to_crs(const count_t nnz,
                      const idx_t *const SMESH_RESTRICT row_idx,
                      const ptrdiff_t n, count_t *const SMESH_RESTRICT rowptr) {
  memset(rowptr, 0, (n + 1) * sizeof(count_t));

#pragma omp parallel for
  for (count_t i = 0; i < nnz; i++) {
    // Check sorted
    assert(i == 0 || row_idx[i - 1] <= row_idx[i]);

#pragma omp atomic update
    rowptr[row_idx[i] + 1]++;
  }

  // Cumulative sum
  for (ptrdiff_t i = 0; i < n; i++) {
    rowptr[i + 1] += rowptr[i];
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_GRAPH_IMPL_HPP
