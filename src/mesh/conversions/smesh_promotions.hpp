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

template <typename count_t, typename idx_t, typename geom_t>
int p1_to_p2(const enum ElemType element_type, const ptrdiff_t n_elements,
             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                 p1_elements,
             const int spatial_dim, const ptrdiff_t p1_n_nodes,
             const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p1_points,
             const count_t *const SMESH_RESTRICT n2n_ptr,
             const idx_t *const SMESH_RESTRICT n2n_idx,
             idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT p2_elements,
             geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p2_points);

template <typename idx_t, typename geom_t>
void quad4_to_hex8_extrude(
    const ptrdiff_t nsides, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t nlayers, const geom_t height,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT extruded_points);

template <typename idx_t, typename geom_t>
int tri3_to_wedge6_extrude(
    const ptrdiff_t nsides, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tri3_elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t nlayers, const geom_t height,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT wedge6_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT extruded_points);
} // namespace smesh

#endif // SMESH_PROMOTIONS_HPP