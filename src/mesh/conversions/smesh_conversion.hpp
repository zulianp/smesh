#ifndef SMESH_CONVERSION_HPP
#define SMESH_CONVERSION_HPP

#include "smesh_base.hpp"

#include <stddef.h>

namespace smesh {
template <typename idx_t>
void mesh_hex8_to_6x_tet4(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements);

/// LocalSideTable HEX8 faces. Each quad is wound with an outward normal.
inline constexpr int hex8_outward_face_nodes[6][4] = {
    {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6},
    {3, 0, 4, 7}, {3, 2, 1, 0}, {4, 5, 6, 7}};

/// One hex → 24 tets: 4 per face, meeting at the face center and the cell center.
/// Face edges are reversed from `hex8_outward_face_nodes` so the cell center is
/// on the positive side of each face triangle.
template <typename idx_t>
void mesh_hex8_emit_24_tet4(
    const idx_t *corner, const idx_t *face_center, const idx_t cell_center,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements,
    const ptrdiff_t tet_begin);

template <typename idx_t>
void mesh_tet15_to_4x_hex8(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements);

template <typename idx_t>
void mesh_wedge6_to_3x_tet4(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT wedge6_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements);

template <typename idx_t>
void mesh_pyramid5_to_2x_tet4(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT pyramid5_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements);

template <typename idx_t>
void mesh_quad4_to_2x_tri3(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tri3_elements);
} // namespace smesh

#endif // SMESH_CONVERSION_HPP