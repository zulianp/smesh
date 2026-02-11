#ifndef SMESH_MULTIBLOCK_GRAPH_HPP
#define SMESH_MULTIBLOCK_GRAPH_HPP

#include "smesh_types.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {

    template <typename idx_t, typename count_t, typename element_idx_t = idx_t>
    int create_multiblock_n2e(const block_idx_t n_blocks,
                             const enum ElemType element_types[],
                             const ptrdiff_t n_elements[], idx_t **const SMESH_RESTRICT elements[],
                             const ptrdiff_t n_nodes,
                             block_idx_t **out_block_number, count_t **out_n2eptr,
                             element_idx_t **out_elindex);

    template <typename idx_t, typename count_t, typename element_idx_t = idx_t>
    int create_multiblock_crs_graph_from_n2e(
        const block_idx_t n_blocks, const enum ElemType element_types[],
        const ptrdiff_t n_elements[], const ptrdiff_t n_nodes,
        idx_t **const SMESH_RESTRICT elems[],
        const count_t *const SMESH_RESTRICT n2eptr,
        const element_idx_t *const SMESH_RESTRICT elindex,
        const block_idx_t *const SMESH_RESTRICT block_number, count_t **out_rowptr,
        idx_t **out_colidx);

    template <typename idx_t, typename count_t, typename element_idx_t = idx_t>
    int create_multiblock_crs_graph(
        const block_idx_t n_blocks, const enum ElemType element_types[],
        const ptrdiff_t n_elements[], idx_t **const SMESH_RESTRICT elems[],
        const ptrdiff_t n_nodes, count_t **out_rowptr, idx_t **out_colidx);

    template <typename idx_t, typename count_t, typename element_idx_t = idx_t>
    int create_multiblock_crs_graph_upper_triangular_from_n2e(
        const block_idx_t n_blocks, const enum ElemType element_types[],
        const ptrdiff_t n_elements[], const ptrdiff_t n_nodes,
        idx_t **const SMESH_RESTRICT elems[],
        const count_t *const SMESH_RESTRICT n2eptr,
        const element_idx_t *const SMESH_RESTRICT elindex,
        const block_idx_t *const SMESH_RESTRICT block_number, count_t **out_rowptr,
        idx_t **out_colidx);

    template <typename idx_t, typename count_t, typename element_idx_t = idx_t>
    int create_multiblock_crs_graph_upper_triangular(
        const block_idx_t n_blocks, const enum ElemType element_types[],
        const ptrdiff_t n_elements[], idx_t **const SMESH_RESTRICT elems[],
        const ptrdiff_t n_nodes, count_t **out_rowptr, idx_t **out_colidx);

} // namespace smesh

#endif // SMESH_MULTIBLOCK_GRAPH_HPP
