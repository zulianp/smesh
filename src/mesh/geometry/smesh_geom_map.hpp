#ifndef SMESH_GEOM_MAP_HPP
#define SMESH_GEOM_MAP_HPP

#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

/// Relative tolerance used when the caller does not pass one.
inline geom_t geom_map_default_rel_tol() {
    return sizeof(geom_t) <= 4 ? geom_t(1e-5) : geom_t(1e-12);
}

/// Classify every element in a block from geometry. Empty input is
/// `ISOPARAMETRIC`. Linear simplices are `AFFINE`. Pyramids are
/// `ISOPARAMETRIC`. `AXIS_ALIGNED` is HEX/QUAD only.
enum GeomMap detect_geom_map(const enum ElemType                                     type,
                             const ptrdiff_t                                         nelements,
                             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                             const int                                               sdim,
                             const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                             const geom_t                                            rel_tol);

inline enum GeomMap detect_geom_map(const enum ElemType                                     type,
                                    const ptrdiff_t                                         nelements,
                                    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                    const int                                               sdim,
                                    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
    return detect_geom_map(type, nelements, elements, sdim, points, geom_map_default_rel_tol());
}

}  // namespace smesh

#endif  // SMESH_GEOM_MAP_HPP
