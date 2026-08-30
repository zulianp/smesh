#ifndef SMESH_SSQUAD4_MESH_HPP
#define SMESH_SSQUAD4_MESH_HPP

#include "smesh_base.hpp"

#include <stddef.h>

namespace smesh {

template <typename idx_t>
int ssquad4_to_standard_quad4_mesh(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements);

template <typename idx_t, typename geom_t>
int ssquad4_fill_points(const int                                                level,
                        const ptrdiff_t                                          nelements,
                        const int                                                n_dims,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                        geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points);

template <typename idx_t>
int ssedge_to_standard_edge2_mesh(const int                                               level,
                                  const ptrdiff_t                                         nelements,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                  idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             edge2_elements);

} // namespace smesh

#endif // SMESH_SSQUAD4_MESH_HPP
