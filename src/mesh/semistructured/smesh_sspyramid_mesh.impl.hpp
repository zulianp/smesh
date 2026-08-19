#ifndef SMESH_SSPYRAMID_MESH_IMPL_HPP
#define SMESH_SSPYRAMID_MESH_IMPL_HPP

#include "smesh_sspyramid.hpp"
#include "smesh_sspyramid_mesh.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
int sspyramid_fill_points(const int                                                level,
                          const ptrdiff_t                                          nelements,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                          const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                          geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points) {
    const int corners[5] = {sspyramid_lidx(level, 0, 0, 0),
                            sspyramid_lidx(level, level, 0, 0),
                            sspyramid_lidx(level, level, level, 0),
                            sspyramid_lidx(level, 0, level, 0),
                            sspyramid_lidx(level, 0, 0, level)};

    for (int k = 0; k <= level; ++k) {
        const geom_t zeta = geom_t(k) / geom_t(level);
        const int    span = level - k;
        for (int j = 0; j <= span; ++j) {
            for (int i = 0; i <= span; ++i) {
                geom_t xi  = geom_t(0);
                geom_t eta = geom_t(0);
                if (span > 0) {
                    xi  = geom_t(i) / geom_t(span);
                    eta = geom_t(j) / geom_t(span);
                }
                const geom_t n0   = (geom_t(1) - xi) * (geom_t(1) - eta);
                const geom_t n1   = xi * (geom_t(1) - eta);
                const geom_t n2   = xi * eta;
                const geom_t n3   = (geom_t(1) - xi) * eta;
                const int    lidx = sspyramid_lidx(level, i, j, k);
                for (int d = 0; d < 3; d++) {
                    for (ptrdiff_t e = 0; e < nelements; e++) {
                        const geom_t base = n0 * macro_mesh_points[d][elements[corners[0]][e]] +
                                            n1 * macro_mesh_points[d][elements[corners[1]][e]] +
                                            n2 * macro_mesh_points[d][elements[corners[2]][e]] +
                                            n3 * macro_mesh_points[d][elements[corners[3]][e]];
                        const geom_t apex = macro_mesh_points[d][elements[corners[4]][e]];
                        points[d][elements[lidx][e]] = (geom_t(1) - zeta) * base + zeta * apex;
                    }
                }
            }
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
