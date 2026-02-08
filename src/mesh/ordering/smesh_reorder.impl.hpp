#ifndef SMESH_REORDER_IMPL_HPP
#define SMESH_REORDER_IMPL_HPP

#include "smesh_reorder.hpp"

namespace smesh {

template <typename geom_t>
int mesh_block_reorder(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT in_elements,
    const element_idx_t *const SMESH_RESTRICT e2e_in2out,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT out_elements) {
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_REORDER_IMPL_HPP
