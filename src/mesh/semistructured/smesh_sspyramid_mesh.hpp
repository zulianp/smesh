#ifndef SMESH_SSPYRAMID_MESH_HPP
#define SMESH_SSPYRAMID_MESH_HPP

#include "smesh_base.hpp"

#include <stddef.h>

namespace smesh {

template <typename idx_t, typename geom_t>
int sspyramid_fill_points(const int                                                level,
                          const ptrdiff_t                                          nelements,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                          const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                          geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points);

/// Explode a PROTEUS_PYRAMID* SS block into linear PYRAMID5 and TET4 connectivity.
/// The two output SoA buffers must be pre-allocated to (sspyramid_n_pyr(L)*nelements)
/// and (sspyramid_n_tet(L)*nelements) columns respectively.
/// No new nodes are introduced; global node ids come from \p elements.
template <typename idx_t>
int sspyramid_to_pyramid5_and_tet4(const int                                               level,
                                   const ptrdiff_t                                         nelements,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                   idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             pyr5_elements,
                                   idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             tet4_elements);

}  // namespace smesh

#endif
