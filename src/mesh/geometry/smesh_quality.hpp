#ifndef SMESH_QUALITY_HPP
#define SMESH_QUALITY_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

/// Isotropic mean-ratio in (0, 1]. Inverted or degenerate → 0.
/// TRI: 4√3 Area / Σ e². TET: 12 (3|V|)^{2/3} / Σ e². QUAD: min of the
/// better diagonal's two corner triangles.
template <typename idx_t, typename real_t>
real_t mesh_elem_mean_ratio(const enum ElemType                                     element_type,
                            const int                                               sdim,
                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                            const real_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                            const ptrdiff_t                                         e);

template <typename idx_t, typename real_t>
int mesh_element_quality(const enum ElemType                                     element_type,
                         const ptrdiff_t                                         n_elements,
                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                         const int                                               sdim,
                         const real_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                         real_t *const SMESH_RESTRICT                            q);

template <typename real_t>
real_t mesh_quality_min(const ptrdiff_t n_elements, const real_t *const SMESH_RESTRICT q);

/// Clip `dx = x - x0` so `|n·dx| <= max_normal_dev` (if > 0) and `|dx| <= max_abs_dev`
/// (if > 0). Negative limits skip that clip. `nz` may be null when `sdim < 3`.
template <typename real_t>
void mesh_clamp_to_band(const int    sdim,
                        real_t      *x,
                        real_t      *y,
                        real_t      *z,
                        const real_t x0,
                        const real_t y0,
                        const real_t z0,
                        const real_t nx,
                        const real_t ny,
                        const real_t nz,
                        const real_t max_abs_dev,
                        const real_t max_normal_dev);

/// Area-weighted vertex normals from TRI/QUAD faces. Interior / unused nodes stay 0.
/// Output arrays are length `n_nodes`. `nz` may be null when `sdim < 3`.
template <typename idx_t, typename geom_t>
int mesh_vertex_normals_from_faces(const enum ElemType                                     element_type,
                                   const ptrdiff_t                                         n_elements,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                   const int                                               sdim,
                                   const ptrdiff_t                                         n_nodes,
                                   const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                   geom_t *const SMESH_RESTRICT                            nx,
                                   geom_t *const SMESH_RESTRICT                            ny,
                                   geom_t *const SMESH_RESTRICT                            nz);

/// Move masked query nodes to the closest point on the TRI/QUAD surface
/// (`surf_points`, `elements` index `n_surf_nodes`). `query_mask` null → all.
/// `max_move` > 0 skips a query whose closest point is farther than that.
template <typename idx_t, typename geom_t>
int mesh_project_to_surface(const enum ElemType                                     element_type,
                            const ptrdiff_t                                         n_elements,
                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                            const int                                               sdim,
                            const ptrdiff_t                                         n_surf_nodes,
                            const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT surf_points,
                            const ptrdiff_t                                         n_query,
                            const uint8_t *const SMESH_RESTRICT                     query_mask,
                            geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT      query_points,
                            const geom_t                                            max_move = 0);

/// Move masked query nodes to the closest point on segments `(e0[k], e1[k])`
/// whose endpoints live in `seg_points`.
template <typename idx_t, typename geom_t>
int mesh_project_to_segments(const int                                               sdim,
                             const ptrdiff_t                                         n_seg,
                             const idx_t *const SMESH_RESTRICT                       e0,
                             const idx_t *const SMESH_RESTRICT                       e1,
                             const ptrdiff_t                                         n_seg_nodes,
                             const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT seg_points,
                             const ptrdiff_t                                         n_query,
                             const uint8_t *const SMESH_RESTRICT                     query_mask,
                             geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT      query_points,
                             const geom_t                                            max_move = 0);

}  // namespace smesh

#endif
