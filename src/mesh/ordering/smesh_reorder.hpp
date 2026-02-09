#ifndef SMESH_REORDER_HPP
#define SMESH_REORDER_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename geom_t, typename element_idx_t>
int mesh_block_reorder(
    const int nxe, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT in_elements,
    const element_idx_t *const SMESH_RESTRICT e2e_gather,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT out_elements);
} // namespace smesh

#endif // SMESH_REORDER_HPP
