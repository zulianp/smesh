#ifndef SMESH_SMOOTH_HPP
#define SMESH_SMOOTH_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

/// lock: 0 free (tangential Jacobi), 1 crease (1D along sharp edge), 2 pin.
/// `e0`/`e1` are sharp-edge node pairs (`n_sharp` of them). `n2n` is the full
/// undirected CRS. Jacobi uses two point buffers internally.
/// Optional `x0` + `max_abs_dev` / `max_normal_dev` clip surface (or all, if
/// `surface` is null) nodes after each write. `lock==2` copies `x0` when set.
/// `nx`/`ny`/`nz` are unit vertex normals; if null, no tangent projection.
/// Surface-marked nodes average only surface neighbors.
template <typename idx_t, typename count_t, typename geom_t>
int mesh_smooth_feature(const int                                               sdim,
                        const ptrdiff_t                                         n_nodes,
                        geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT      points,
                        const count_t *const SMESH_RESTRICT                     rowptr,
                        const idx_t *const SMESH_RESTRICT                       colidx,
                        const uint8_t *const SMESH_RESTRICT                     lock,
                        const ptrdiff_t                                         n_sharp,
                        const idx_t *const SMESH_RESTRICT                       e0,
                        const idx_t *const SMESH_RESTRICT                       e1,
                        const int                                               n_iters,
                        const geom_t                                            lambda,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT x0 = nullptr,
                        const uint8_t *const SMESH_RESTRICT                     surface = nullptr,
                        const geom_t                                            max_abs_dev = 0,
                        const geom_t                                            max_normal_dev = 0,
                        const geom_t *const SMESH_RESTRICT                     nx = nullptr,
                        const geom_t *const SMESH_RESTRICT                     ny = nullptr,
                        const geom_t *const SMESH_RESTRICT                     nz = nullptr);

}  // namespace smesh

#endif
