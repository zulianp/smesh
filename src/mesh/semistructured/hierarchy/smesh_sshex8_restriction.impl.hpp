#ifndef SMESH_SSHEX8_RESTRICTION_IMPL_HPP
#define SMESH_SSHEX8_RESTRICTION_IMPL_HPP

#include "smesh_sshex8_restriction.hpp"
#include "smesh_sshex8.hpp"

namespace smesh {

template <typename idx_t, typename real_t>
int sshex8_hierarchical_restriction(
    int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const real_t *const SMESH_RESTRICT from,
    real_t *const SMESH_RESTRICT to) {
#pragma omp parallel
  {
    const int nxe = sshex8_nxe(level);
    real_t **e_from = (real_t **)malloc(vec_size * sizeof(real_t *));
    real_t **e_to = (real_t **)malloc(vec_size * sizeof(real_t *));
    u16 *weight = (u16 *)malloc(nxe * sizeof(u16));

    for (int d = 0; d < vec_size; d++) {
      e_from[d] = (real_t *)malloc(nxe * sizeof(real_t));
      e_to[d] = (real_t *)malloc(8 * sizeof(real_t));
    }

    idx_t *ev = (idx_t *)malloc(nxe * sizeof(idx_t));

    const int corners[8] = {
        // Bottom
        sshex8_lidx(level, 0, 0, 0), sshex8_lidx(level, level, 0, 0),
        sshex8_lidx(level, level, level, 0), sshex8_lidx(level, 0, level, 0),
        // Top
        sshex8_lidx(level, 0, 0, level), sshex8_lidx(level, level, 0, level),
        sshex8_lidx(level, level, level, level),
        sshex8_lidx(level, 0, level, level)};

#pragma omp for
    for (ptrdiff_t e = 0; e < nelements; e++) {
      {
        // Gather elemental data
        for (int d = 0; d < nxe; d++) {
          ev[d] = elements[d][e];
        }

        for (int d = 0; d < vec_size; d++) {
          for (int v = 0; v < nxe; v++) {
            e_from[d][v] = from[ev[v] * vec_size + d];
            SMESH_ASSERT(e_from[d][v] == e_from[d][v]);
          }
        }

        for (int v = 0; v < nxe; v++) {
          weight[v] = element_to_node_incidence_count[ev[v]];
        }

        for (int d = 0; d < vec_size; d++) {
          for (int v = 0; v < 8; v++) {
            e_to[d][v] = 0;
          }
        }

        const real_t h = 1. / level;
        // Iterate over structrued grid (nodes)
        for (int zi = 0; zi < level + 1; zi++) {
          for (int yi = 0; yi < level + 1; yi++) {
            for (int xi = 0; xi < level + 1; xi++) {
              const int lidx = sshex8_lidx(level, xi, yi, zi);

              const real_t x = xi * h;
              const real_t y = yi * h;
              const real_t z = zi * h;

              // Evaluate Hex8 basis functions at x, y, z
              const real_t xm = (1 - x);
              const real_t ym = (1 - y);
              const real_t zm = (1 - z);

              real_t f[8];
              f[0] = xm * ym * zm; // (0, 0, 0)
              f[1] = x * ym * zm;  // (1, 0, 0)
              f[2] = x * y * zm;   // (1, 1, 0)
              f[3] = xm * y * zm;  // (0, 1, 0)
              f[4] = xm * ym * z;  // (0, 0, 1)
              f[5] = x * ym * z;   // (1, 0, 1)
              f[6] = x * y * z;    // (1, 1, 1)
              f[7] = xm * y * z;   // (0, 1, 1)

              for (int d = 0; d < vec_size; d++) {
                const real_t val = e_from[d][lidx] / weight[lidx];
                for (int i = 0; i < 8; i++) {
                  e_to[d][i] += f[i] * val;
                }
              }
            }
          }
        }

        for (int d = 0; d < vec_size; d++) {
          for (int i = 0; i < 8; i++) {
            const int idx = ev[corners[i]] * vec_size + d;
#pragma omp atomic update
            to[idx] += e_to[d][i];
          }
        }
      }
    }

    for (int d = 0; d < vec_size; d++) {
      free(e_to[d]);
      free(e_from[d]);
    }

    free(e_from);
    free(e_to);
    free(ev);
    free(weight);
  }

  return SMESH_SUCCESS;
}
} // namespace smesh

#endif // SMESH_SSHEX8_RESTRICTION_IMPL_HPP
