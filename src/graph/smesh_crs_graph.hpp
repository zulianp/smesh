#ifndef SMESH_CRS_GRAPH_HPP
#define SMESH_CRS_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(const ptrdiff_t nelements, const ptrdiff_t nnodes,
               const int nnodesxelem,
               idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
               count_t **out_n2eptr, element_idx_t **out_elindex);

template <typename idx_t, typename count_t>
int create_crs_graph_for_elem_type(
    const enum ElemType element_type, const ptrdiff_t nelements,
    const ptrdiff_t nnodes,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx);

template <typename idx_t, typename count_t>
int create_crs_graph(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                     idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                     count_t **out_rowptr, idx_t **out_colidx);

template <typename idx_t, typename count_t>
int create_crs_graph_3(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                       idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                       count_t **out_rowptr, idx_t **out_colidx);

template <typename idx_t, typename count_t, typename real_t>
int block_crs_to_crs(const ptrdiff_t nnodes, const int block_size,
                     const count_t *const SMESH_RESTRICT block_rowptr,
                     const idx_t *const SMESH_RESTRICT block_colidx,
                     const real_t *const SMESH_RESTRICT block_values,
                     count_t *const SMESH_RESTRICT rowptr,
                     idx_t *const SMESH_RESTRICT colidx,
                     real_t *const SMESH_RESTRICT values);

template <typename idx_t, typename count_t>
int crs_to_coo(const ptrdiff_t n, const count_t *const rowptr,
               idx_t *const SMESH_RESTRICT row_idx);

template <typename idx_t, typename count_t>
int sorted_coo_to_crs(const count_t nnz,
                      const idx_t *const SMESH_RESTRICT row_idx,
                      const ptrdiff_t n, count_t *const SMESH_RESTRICT rowptr);

// for crs insertion
template <typename idx_t>
idx_t find_idx(const idx_t key, const idx_t *const SMESH_RESTRICT arr,
               idx_t size);
// idx_t find_idx_binary_search(const idx_t key, const idx_t *arr, idx_t size);

template <typename idx_t, typename count_t, typename element_idx_t>
int create_dual_graph(const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
                      const enum ElemType element_type,
                      idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                      count_t **out_rowptr, element_idx_t **out_colidx);

template <typename idx_t, typename count_t>
int crs_graph_block_to_scalar(const ptrdiff_t nnodes, const int block_size,
                              const count_t *const SMESH_RESTRICT block_rowptr,
                              const idx_t *const SMESH_RESTRICT block_colidx,
                              count_t *const SMESH_RESTRICT rowptr,
                              idx_t *const SMESH_RESTRICT colidx);

template <typename idx_t, typename count_t>
int create_crs_graph_from_element(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, int nxe,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx);

template <typename idx_t, typename count_t>
int create_crs_graph_upper_triangular_from_element(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, int nxe,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx);

} // namespace smesh

#endif // SMESH_CRS_GRAPH_HPP
