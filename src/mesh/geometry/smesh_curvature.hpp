#ifndef SMESH_CURVATURE_HPP
#define SMESH_CURVATURE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

/// Discrete mean curvature |κ| at each node (Meyer cotangent, |T|/3 mixed area).
/// Surface types: TRI3/TRISHELL3, QUAD4/QUADSHELL4. `sharp_node` (nullable) zeros
/// those nodes. sdim==2 uses boundary turning angle; interior nodes are 0.
template <typename idx_t, typename geom_t>
int mesh_node_curvature(const enum ElemType                                     element_type,
                        const ptrdiff_t                                         n_elements,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                        const int                                               sdim,
                        const ptrdiff_t                                         n_nodes,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                        const uint8_t *const SMESH_RESTRICT                     sharp_node,
                        geom_t *const SMESH_RESTRICT                            kappa);

/// `h = clip(2π / (|κ| N_cpr), h_min, h_max)`. κ≈0 → h_max. h_min/h_max ≤ 0
/// are replaced by `bbox_diag * 1e-4` and `bbox_diag`.
template <typename geom_t>
int mesh_size_from_curvature(const ptrdiff_t              n_nodes,
                             const geom_t *const SMESH_RESTRICT kappa,
                             const geom_t                     cells_per_radius,
                             const geom_t                     h_min,
                             const geom_t                     h_max,
                             const geom_t                     bbox_diag,
                             geom_t *const SMESH_RESTRICT     h);

/// Sagitta estimator `η = |κ| ℓ² / 8`. Sets `h_i = sqrt(8 ε / κ_i)` so the
/// kernel mark `ℓ > h` is `η > ε`. κ≈0 stays unconstrained (`h = 1e6 diag`)
/// so flats/interior are not length-refined; 2:1 skips those neighbors.
/// `geom_error ≤ 0` uses `ε = π² bbox_diag / (2 N_cpr²)`. `h_max` clips only
/// curved nodes.
template <typename geom_t>
int mesh_size_from_curvature_error(const ptrdiff_t              n_nodes,
                                   const geom_t *const SMESH_RESTRICT kappa,
                                   const geom_t                     cells_per_radius,
                                   const geom_t                     geom_error,
                                   const geom_t                     h_min,
                                   const geom_t                     h_max,
                                   const geom_t                     bbox_diag,
                                   geom_t *const SMESH_RESTRICT     h);

/// Jacobi 2:1 on h: neighboring nodes differ by at most a factor two.
/// Unconstrained nodes (h ≫ min curved h) take `min(h, 2 h_j)` from
/// constrained neighbors (fine-to-coarse halo). Constrained neighbors
/// also apply `max(h, h_j/2)` so the fine side cannot jump more than 2×.
template <typename idx_t, typename count_t, typename geom_t>
int mesh_grade_size_2to1(const ptrdiff_t                     n_nodes,
                         const count_t *const SMESH_RESTRICT rowptr,
                         const idx_t *const SMESH_RESTRICT   colidx,
                         const int                           n_sweeps,
                         geom_t *const SMESH_RESTRICT        h);

template <typename idx_t>
void mesh_mark_nodes_from_edges(const ptrdiff_t                   n_nodes,
                                const ptrdiff_t                   n_edges,
                                const idx_t *const SMESH_RESTRICT e0,
                                const idx_t *const SMESH_RESTRICT e1,
                                uint8_t *const SMESH_RESTRICT     mark);

}  // namespace smesh

#endif
