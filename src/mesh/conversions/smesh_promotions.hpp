#ifndef SMESH_PROMOTIONS_HPP
#define SMESH_PROMOTIONS_HPP

#include "smesh_base.hpp"

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
void mesh_tet4_to_tet15(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements,
    const count_t *const SMESH_RESTRICT n2n_upper_triangular_ptr,
    const idx_t *const SMESH_RESTRICT n2n_upper_triangular_idx,
    const element_idx_t *const SMESH_RESTRICT e2e_table,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_elements,
    ptrdiff_t *const SMESH_RESTRICT tet15_n_nodes);

template <typename idx_t, typename count_t, typename geom_t>
void mesh_tet4_to_tet15_points(
    const ptrdiff_t n_elements, const ptrdiff_t tet4_n_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_pts,
    const count_t *const SMESH_RESTRICT n2n_upper_triangular_ptr,
    const idx_t *const SMESH_RESTRICT n2n_upper_triangular_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_pts);


} // namespace smesh

#endif // SMESH_PROMOTIONS_HPP