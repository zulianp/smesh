#ifndef SMESH_SSWEDGE_MESH_IMPL_HPP
#define SMESH_SSWEDGE_MESH_IMPL_HPP

#include "smesh_sswedge.hpp"
#include "smesh_sswedge_mesh.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
int sswedge_fill_points(const int                                                level,
                        const ptrdiff_t                                          nelements,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                        geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points) {
    const int corners[6] = {sswedge_lidx(level, 0, 0, 0),
                            sswedge_lidx(level, level, 0, 0),
                            sswedge_lidx(level, 0, level, 0),
                            sswedge_lidx(level, 0, 0, level),
                            sswedge_lidx(level, level, 0, level),
                            sswedge_lidx(level, 0, level, level)};
    const geom_t h = geom_t(1) / geom_t(level);

    for (int z = 0; z <= level; ++z) {
        for (int y = 0; y <= level; ++y) {
            for (int x = 0; x <= level - y; ++x) {
                const geom_t fx   = geom_t(x) * h;
                const geom_t fy   = geom_t(y) * h;
                const geom_t fz   = geom_t(z) * h;
                const geom_t ftri = geom_t(1) - fx - fy;
                const int    lidx = sswedge_lidx(level, x, y, z);
                for (int d = 0; d < 3; d++) {
                    for (ptrdiff_t e = 0; e < nelements; e++) {
                        const geom_t btm = ftri * macro_mesh_points[d][elements[corners[0]][e]] +
                                           fx * macro_mesh_points[d][elements[corners[1]][e]] +
                                           fy * macro_mesh_points[d][elements[corners[2]][e]];
                        const geom_t top = ftri * macro_mesh_points[d][elements[corners[3]][e]] +
                                           fx * macro_mesh_points[d][elements[corners[4]][e]] +
                                           fy * macro_mesh_points[d][elements[corners[5]][e]];
                        points[d][elements[lidx][e]] = (geom_t(1) - fz) * btm + fz * top;
                    }
                }
            }
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
