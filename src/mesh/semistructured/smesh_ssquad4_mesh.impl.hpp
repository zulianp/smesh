#ifndef SMESH_SSQUAD4_MESH_IMPL_HPP
#define SMESH_SSQUAD4_MESH_IMPL_HPP

#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_mesh.hpp"

namespace smesh {

template <typename idx_t>
int ssquad4_to_standard_quad4_mesh(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements) {
  const int txe = ssquad4_txe(level);

  int lnode[4];
  for (int yi = 0; yi < level; yi++) {
    for (int xi = 0; xi < level; xi++) {
      lnode[0] = ssquad4_lidx(level, xi, yi);
      lnode[1] = ssquad4_lidx(level, xi + 1, yi);
      lnode[2] = ssquad4_lidx(level, xi + 1, yi + 1);
      lnode[3] = ssquad4_lidx(level, xi, yi + 1);

      int le = yi * level + xi;
      SMESH_ASSERT(le < txe);

      for (int l = 0; l < 4; l++) {
        for (ptrdiff_t e = 0; e < nelements; e++) {
          idx_t node = elements[lnode[l]][e];
          quad4_elements[l][e * txe + le] = node;
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int ssquad4_fill_points(const int                                                level,
                        const ptrdiff_t                                          nelements,
                        const int                                                n_dims,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                        geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points) {
    const int corners[4] = {ssquad4_lidx(level, 0, 0),
                            ssquad4_lidx(level, level, 0),
                            ssquad4_lidx(level, level, level),
                            ssquad4_lidx(level, 0, level)};

    const geom_t h = geom_t(1) / geom_t(level);

    for (int yi = 0; yi <= level; ++yi) {
        for (int xi = 0; xi <= level; ++xi) {
            const geom_t fx = geom_t(xi) * h;
            const geom_t fy = geom_t(yi) * h;
            const geom_t f[4] = {(geom_t(1) - fx) * (geom_t(1) - fy),
                                 fx * (geom_t(1) - fy),
                                 fx * fy,
                                 (geom_t(1) - fx) * fy};
            const int    lidx = ssquad4_lidx(level, xi, yi);

            for (int d = 0; d < n_dims; d++) {
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

    return SMESH_SUCCESS;
}

template <typename idx_t>
int ssedge_to_standard_edge2_mesh(const int                                               level,
                                  const ptrdiff_t                                         nelements,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                  idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             edge2_elements) {
    const int txe = ssedge_txe(level);
    for (int xi = 0; xi < level; ++xi) {
        const int n0 = ssedge_lidx(level, xi);
        const int n1 = ssedge_lidx(level, xi + 1);
        SMESH_ASSERT(xi < txe);
        for (ptrdiff_t e = 0; e < nelements; ++e) {
            edge2_elements[0][e * txe + xi] = elements[n0][e];
            edge2_elements[1][e * txe + xi] = elements[n1][e];
        }
    }
    return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_SSQUAD4_MESH_IMPL_HPP
