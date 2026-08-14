#ifndef SMESH_SSTET4_MESH_IMPL_HPP
#define SMESH_SSTET4_MESH_IMPL_HPP

#include "smesh_sstet4.hpp"
#include "smesh_sstet4_mesh.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
int sstet4_fill_points(const int                                                level,
                       const ptrdiff_t                                          nelements,
                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                       const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                       geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points) {
    const int corners[4] = {sstet4_lidx(level, 0, 0, 0),
                            sstet4_lidx(level, level, 0, 0),
                            sstet4_lidx(level, 0, level, 0),
                            sstet4_lidx(level, 0, 0, level)};

    const geom_t h = geom_t(1) / geom_t(level);

    for (int z = 0; z <= level; ++z) {
        for (int y = 0; y <= level - z; ++y) {
            for (int x = 0; x <= level - z - y; ++x) {
                const geom_t fx = geom_t(x) * h;
                const geom_t fy = geom_t(y) * h;
                const geom_t fz = geom_t(z) * h;
                const geom_t f0 = geom_t(1) - fx - fy - fz;
                const geom_t f[4] = {f0, fx, fy, fz};
                const int    lidx = sstet4_lidx(level, x, y, z);

                for (int d = 0; d < 3; d++) {
                    for (ptrdiff_t e = 0; e < nelements; e++) {
                        geom_t acc = 0;
                        for (int c = 0; c < 4; c++) {
                            acc += macro_mesh_points[d][elements[corners[c]][e]] * f[c];
                        }
                        points[d][elements[lidx][e]] = acc;
                    }
                }
            }
        }
    }

    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif  // SMESH_SSTET4_MESH_IMPL_HPP
