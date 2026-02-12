#ifndef SMESH_REORDER_HPP
#define SMESH_REORDER_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename element_idx_t>
int mesh_block_reorder(
    const int nxe, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT in_elements,
    const element_idx_t *const SMESH_RESTRICT e2e_gather,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT out_elements);

template <typename idx_t>
int mesh_block_renumber_element_nodes(
    const int nxe, const ptrdiff_t n_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *const SMESH_RESTRICT next_node_idx,
    idx_t *const SMESH_RESTRICT n2n_scatter);

template <typename idx_t, typename T>
int reorder_scatter(const ptrdiff_t n,
                    const idx_t *const SMESH_RESTRICT scatter,
                    const T *const SMESH_RESTRICT in_array,
                    T *const SMESH_RESTRICT out_array);

template <typename idx_t, typename T>
int reorder_gather(const ptrdiff_t n,
                   const idx_t *const SMESH_RESTRICT gather,
                   const T *const SMESH_RESTRICT in_array,
                   T *const SMESH_RESTRICT out_array);
} // namespace smesh

#endif // SMESH_REORDER_HPP
