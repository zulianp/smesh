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
} // namespace smesh

#endif // SMESH_CONVERSION_HPP