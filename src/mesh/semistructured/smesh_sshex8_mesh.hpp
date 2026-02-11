#ifndef SMESH_SSHEX8_MESH_HPP
#define SMESH_SSHEX8_MESH_HPP

#include "smesh_base.hpp"

#include <stddef.h>

namespace smesh {

template <typename idx_t>
int sshex8_to_standard_hex8_mesh(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements);

template <typename idx_t, typename geom_t>
int sshex8_fill_points(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
    geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points);

template <typename idx_t, typename geom_t, typename ref_t>
int sshex8_fill_points_1D_map(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
    const ref_t *const SMESH_RESTRICT ref_points,
    geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points);

} // namespace smesh

#endif // SMESH_SSHEX8_MESH_HPP