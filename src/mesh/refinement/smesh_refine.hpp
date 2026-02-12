#ifndef SMESH_REFINE_HPP
#define SMESH_REFINE_HPP

#include "smesh_base.hpp"

namespace smesh {

template <typename idx_t, typename count_t,
          typename geom_t>
int mesh_refine(
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT coarse_elements,
    const int spatial_dim, const ptrdiff_t n_coarse_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT coarse_points,
    const count_t *const SMESH_RESTRICT n2n_ptr,
    const idx_t *const SMESH_RESTRICT n2n_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT refined_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT refined_points);

} // namespace smesh

#endif // SMESH_REFINE_HPP
