#ifndef SMESH_GRAPH_HPP
#define SMESH_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {

/**
 * @brief Build a node-to-element adjacency in CSR form.
 *
 * Produces a compressed sparse row (CSR) structure mapping each node to the
 * incident element indices:
 * - `out_n2eptr` is an array of length `nnodes + 1`
 * - `out_elindex` is an array of length `(*out_n2eptr)[nnodes]`
 *
 * @tparam idx_t          Node index type used in `elems`.
 * @tparam count_t        CSR pointer type (counts/offsets).
 * @tparam element_idx_t  Element index type stored in `out_elindex`.
 *
 * @param nelements     Number of elements.
 * @param nnodes        Number of nodes.
 * @param nnodesxelem   Number of nodes per element (connectivity arity).
 * @param elems         Connectivity in SoA layout: `elems[dof][e]` is the node
 *                     index of local dof `dof` of element `e`.
 * @param out_n2eptr    Output CSR row pointer (allocated with `malloc`).
 * @param out_elindex   Output element indices (allocated with `malloc`).
 *
 * @return `SMESH_SUCCESS` on success.
 *
 * @note The caller owns `*out_n2eptr` and `*out_elindex` and must `free()`
 * them.
 * @note Preconditions: all `elems[*][e]` satisfy `0 <= node < nnodes`.
 */
template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(const ptrdiff_t nelements, const ptrdiff_t nnodes,
               const int nnodesxelem,
               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
               count_t **out_n2eptr, element_idx_t **out_elindex);

template <typename count_t, typename element_idx_t>
int sort_n2e(const ptrdiff_t nnodes, const count_t *const SMESH_RESTRICT n2eptr,
             element_idx_t *const SMESH_RESTRICT elindex);

template <typename idx_t, typename count_t, typename element_idx_t>
 int
n2n_from_n2e(const ptrdiff_t nelements, const ptrdiff_t nnodes,
             const int nnodesxelem,
             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
             const count_t *const SMESH_RESTRICT n2eptr,
             const element_idx_t *const SMESH_RESTRICT elindex,
             count_t **out_rowptr, idx_t **out_colidx);

    /**
     * @brief Build a node adjacency graph (CSR) induced by mesh elements.
     *
     * For each node, neighbors are the unique set of nodes appearing in any
     * element incident to the node (including the node itself).
     *
     * @tparam idx_t    Node index type.
     * @tparam count_t  CSR pointer type (counts/offsets).
     *
     * @param element_type Element type used to interpret `elems` (see
     * `ElemType`).
     * @param nelements    Number of elements.
     * @param nnodes       Number of nodes.
     * @param elems        Connectivity in SoA layout: `elems[local_node][e]`.
     * @param out_rowptr   Output CSR row pointer, length `nnodes + 1`
     * (malloc'ed).
     * @param out_colidx   Output CSR column indices, length
     * `(*out_rowptr)[nnodes]` (malloc'ed). Each row is sorted and unique.
     *
     * @return `SMESH_SUCCESS` on success.
     *
     * @note The caller owns `*out_rowptr` and `*out_colidx` and must `free()`
     * them.
     * @note Implementation may use OpenMP.
     * @note Setting environment variable `SMESH_CRS_FAST_SERIAL=1` selects an
     *       alternative construction path optimized for serial execution.
     */
    template <typename idx_t, typename count_t>
    int create_crs_graph_for_elem_type(
        const enum ElemType element_type, const ptrdiff_t nelements,
        const ptrdiff_t nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
        count_t **out_rowptr, idx_t **out_colidx);

/**
 * @brief Convenience wrapper for `create_crs_graph_for_elem_type(..., TET4,
 * ...)`.
 *
 * @tparam idx_t    Node index type.
 * @tparam count_t  CSR pointer type (counts/offsets).
 *
 * @param nelements  Number of elements.
 * @param nnodes     Number of nodes.
 * @param elems      Connectivity in SoA layout.
 * @param out_rowptr Output CSR row pointer (malloc'ed).
 * @param out_colidx Output CSR column indices (malloc'ed).
 *
 * @return `SMESH_SUCCESS` on success.
 */
template <typename idx_t, typename count_t>
int create_crs_graph(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                         elems,
                     count_t **out_rowptr, idx_t **out_colidx);

/**
 * @brief Convenience wrapper for `create_crs_graph_for_elem_type(..., TRI3,
 * ...)`.
 *
 * @tparam idx_t    Node index type.
 * @tparam count_t  CSR pointer type (counts/offsets).
 */
template <typename idx_t, typename count_t>
int create_crs_graph_3(const ptrdiff_t nelements, const ptrdiff_t nnodes,
                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                           elems,
                       count_t **out_rowptr, idx_t **out_colidx);

/**
 * @brief Expand a block-CSR matrix into scalar-CSR (including values).
 *
 * Input is a block-sparse matrix with `nnodes` block rows/cols. Each nonzero
 * block is a dense `block_size x block_size` block stored in row-major order in
 * `block_values`.
 *
 * Output is the scalar CSR for the expanded matrix of dimension
 * `(nnodes * block_size) x (nnodes * block_size)`.
 *
 * @tparam idx_t    Column index type.
 * @tparam count_t  CSR pointer type (counts/offsets).
 * @tparam real_t   Scalar value type.
 *
 * @param nnodes        Number of block rows/cols.
 * @param block_size    Block dimension.
 * @param block_rowptr  Block CSR row pointer, length `nnodes + 1`.
 * @param block_colidx  Block CSR column indices, length `block_rowptr[nnodes]`.
 * @param block_values  Dense blocks, length
 *                      `block_rowptr[nnodes] * block_size * block_size`.
 * @param rowptr        Output scalar CSR row pointer, length `nnodes*block_size
 * + 1` (must be preallocated).
 * @param colidx        Output scalar CSR colidx, length
 *                      `block_rowptr[nnodes] * block_size * block_size`
 *                      (must be preallocated).
 * @param values        Output scalar CSR values, same length as `colidx`
 *                      (must be preallocated).
 *
 * @return `SMESH_SUCCESS` on success.
 */
template <typename idx_t, typename count_t, typename real_t>
int block_crs_to_crs(const ptrdiff_t nnodes, const int block_size,
                     const count_t *const SMESH_RESTRICT block_rowptr,
                     const idx_t *const SMESH_RESTRICT block_colidx,
                     const real_t *const SMESH_RESTRICT block_values,
                     count_t *const SMESH_RESTRICT rowptr,
                     idx_t *const SMESH_RESTRICT colidx,
                     real_t *const SMESH_RESTRICT values);

/**
 * @brief Convert CRS row pointers to COO row indices.
 *
 * Computes `row_idx[k] = row` for all `k` in `[rowptr[row], rowptr[row+1])`.
 *
 * @tparam idx_t    Row index type for COO.
 * @tparam count_t  CRS pointer type (counts/offsets).
 *
 * @param n       Number of rows.
 * @param rowptr  CRS row pointer, length `n + 1`.
 * @param row_idx Output COO row index array, length `rowptr[n]` (preallocated).
 *
 * @return `SMESH_SUCCESS` on success.
 *
 * @note Implementation may use OpenMP.
 */
template <typename idx_t, typename count_t>
int crs_to_coo(const ptrdiff_t n, const count_t *const rowptr,
               idx_t *const SMESH_RESTRICT row_idx);

/**
 * @brief Build CRS row pointers from sorted COO row indices.
 *
 * Given a COO `row_idx` array sorted nondecreasing by row, constructs `rowptr`.
 *
 * @tparam idx_t    Row index type in COO.
 * @tparam count_t  CRS pointer type (counts/offsets).
 *
 * @param nnz      Number of nonzeros (length of `row_idx`).
 * @param row_idx  Sorted COO row indices, length `nnz`.
 * @param n        Number of rows.
 * @param rowptr   Output CRS row pointer, length `n + 1` (preallocated).
 *
 * @return `SMESH_SUCCESS` on success.
 *
 * @note Preconditions: `row_idx` is sorted; `0 <= row_idx[i] < n`.
 * @note Implementation may use OpenMP atomics.
 */
template <typename idx_t, typename count_t>
int sorted_coo_to_crs(const count_t nnz,
                      const idx_t *const SMESH_RESTRICT row_idx,
                      const ptrdiff_t n, count_t *const SMESH_RESTRICT rowptr);

/**
 * @brief Linear search for `key` in `arr`.
 *
 * @tparam idx_t Element type.
 *
 * @param key  Target value.
 * @param arr  Array to search.
 * @param size Number of elements in `arr`.
 *
 * @return Index `i` with `arr[i] == key`, or `size` if not found.
 */
template <typename idx_t>
idx_t find_idx(const idx_t key, const idx_t *const SMESH_RESTRICT arr,
               idx_t size);
// idx_t find_idx_binary_search(const idx_t key, const idx_t *arr, idx_t size);

/**
 * @brief Build the element dual graph (element adjacency) in CSR form.
 *
 * Two elements are adjacent if they share a full side (facet/edge depending on
 * `element_type`). The output stores, for each element `e`, the list of
 * adjacent element indices.
 *
 * @tparam idx_t          Node index type used in `elems`.
 * @tparam count_t        CSR pointer type (counts/offsets).
 * @tparam element_idx_t  Element index type stored in the adjacency list.
 *
 * @param n_elements   Number of elements.
 * @param n_nodes      Number of nodes.
 * @param element_type Element type (controls the definition of a "side").
 * @param elems        Connectivity in SoA layout.
 * @param out_rowptr   Output CSR row pointer, length `n_elements + 1`
 * (malloc'ed).
 * @param out_colidx   Output adjacency list (malloc'ed). Valid entries are in
 *                     `[0, (*out_rowptr)[n_elements])`.
 *
 * @return `SMESH_SUCCESS` on success.
 *
 * @note The caller owns `*out_rowptr` and `*out_colidx` and must `free()` them.
 */
template <typename idx_t, typename count_t, typename element_idx_t>
int create_dual_graph(const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
                      const enum ElemType element_type,
                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                          elems,
                      count_t **out_rowptr, element_idx_t **out_colidx);

/**
 * @brief Expand a block-CSR sparsity pattern into scalar-CSR (structure only).
 *
 * Same as `block_crs_to_crs` but without values.
 *
 * @tparam idx_t    Column index type.
 * @tparam count_t  CSR pointer type (counts/offsets).
 */
template <typename idx_t, typename count_t>
int crs_graph_block_to_scalar(const ptrdiff_t nnodes, const int block_size,
                              const count_t *const SMESH_RESTRICT block_rowptr,
                              const idx_t *const SMESH_RESTRICT block_colidx,
                              count_t *const SMESH_RESTRICT rowptr,
                              idx_t *const SMESH_RESTRICT colidx);

/**
 * @brief Build a node adjacency graph (CSR) from generic element connectivity.
 *
 * This is the element-type-agnostic entry point: `nxe` controls the number of
 * nodes per element.
 *
 * @tparam idx_t    Node index type.
 * @tparam count_t  CSR pointer type (counts/offsets).
 *
 * @param nelements  Number of elements.
 * @param nnodes     Number of nodes.
 * @param nxe        Nodes per element.
 * @param elems      Connectivity in SoA layout.
 * @param out_rowptr Output CSR row pointer (malloc'ed).
 * @param out_colidx Output CSR column indices (malloc'ed, sorted and unique per
 * row).
 *
 * @return `SMESH_SUCCESS` on success.
 */
template <typename idx_t, typename count_t>
int create_crs_graph_from_element(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, int nxe,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx);

/**
 * @brief Build an upper-triangular node adjacency graph (CSR) from elements.
 *
 * For each node `i`, only neighbors `j` with `j > i` are stored. Rows are
 * sorted and unique.
 *
 * @tparam idx_t    Node index type.
 * @tparam count_t  CSR pointer type (counts/offsets).
 */
template <typename idx_t, typename count_t>
int create_crs_graph_upper_triangular_from_element(
    const ptrdiff_t nelements, const ptrdiff_t nnodes, int nxe,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx);

} // namespace smesh

#endif // SMESH_GRAPH_HPP
