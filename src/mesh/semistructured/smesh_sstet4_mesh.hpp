#ifndef SMESH_SSTET4_MESH_HPP
#define SMESH_SSTET4_MESH_HPP

#include "smesh_base.hpp"

#include <stddef.h>

namespace smesh {

template <typename idx_t>
int sstet4_to_standard_tet4_mesh(const int                                                level,
                                 const ptrdiff_t                                          nelements,
                                 const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                                 idx_t *SMESH_RESTRICT *const SMESH_RESTRICT              tet4_elements);

template <typename idx_t, typename geom_t>
int sstet4_fill_points(const int                                                 level,
                       const ptrdiff_t                                           nelements,
                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT   elements,
                       const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT  macro_mesh_points,
                       geom_t *SMESH_RESTRICT *const SMESH_RESTRICT              points);

template <typename idx_t>
int sstri_to_standard_tri3_mesh(const int                                               level,
                                const ptrdiff_t                                         nelements,
                                const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             tri3_elements);

}  // namespace smesh

#endif  // SMESH_SSTET4_MESH_HPP
