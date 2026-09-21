#ifndef SMESH_IMPROVE_HPP
#define SMESH_IMPROVE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

inline bool improve_type_supported(const enum ElemType t) {
    return t == TET4 || t == TRI3 || t == TRISHELL3 || t == QUAD4 || t == QUADSHELL4;
}

/// Quality-driven local remesh. `lock` is 0/1/2 over nodes; `surface` (nullable)
/// marks band-constrained nodes; `x0` is the snapshot SoA. Output arrays are
/// SMESH_ALLOC'd; free with `mesh_improve_free`.
template <typename idx_t, typename count_t, typename geom_t>
int mesh_improve(const enum ElemType                                      element_type,
                 const ptrdiff_t                                          n_elements,
                 const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                 const int                                                sdim,
                 const ptrdiff_t                                          n_nodes,
                 const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                 const uint8_t *const SMESH_RESTRICT                      lock,
                 const uint8_t *const SMESH_RESTRICT                      surface,
                 const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT x0,
                 const ptrdiff_t                                          n_sharp,
                 const idx_t *const SMESH_RESTRICT                        se0,
                 const idx_t *const SMESH_RESTRICT                        se1,
                 const geom_t                                             q_min,
                 const geom_t                                             max_abs_dev,
                 const geom_t                                             max_normal_dev,
                 const int                                                max_passes,
                 const int                                                allow_split,
                 const int                                                allow_collapse,
                 const int                                                allow_swap,
                 ptrdiff_t                                               *n_elements_out,
                 idx_t                                                 ***elements_out,
                 ptrdiff_t                                               *n_nodes_out,
                 geom_t                                                ***points_out,
                 uint8_t                                                **lock_out,
                 uint8_t                                                **surface_out,
                 geom_t                                                ***x0_out,
                 ptrdiff_t                                               *n_ops_out);

void mesh_improve_free(const int nxe, const int sdim, idx_t **elements, geom_t **points, uint8_t *lock,
                       uint8_t *surface, geom_t **x0);

}  // namespace smesh

#endif
