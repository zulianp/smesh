#ifndef SMESH_SSHEX8_RESTRICTION_IMPL_HPP
#define SMESH_SSHEX8_RESTRICTION_IMPL_HPP

#include "smesh_sshex8.hpp"
#include "smesh_alloc.hpp"
#include "smesh_sshex8_restriction.hpp"

#include <string.h>
#include <stdlib.h>

namespace smesh {

    template <typename idx_t, typename real_t>
    int sshex8_hierarchical_restriction(int                                                     level,
                                        const ptrdiff_t                                         nelements,
                                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                        const u16 *const SMESH_RESTRICT                         element_to_node_incidence_count,
                                        const int                                               vec_size,
                                        const real_t *const SMESH_RESTRICT                      from,
                                        real_t *const SMESH_RESTRICT                            to) {
#pragma omp parallel
        {
            const int nxe    = sshex8_nxe(level);
            real_t  **e_from = (real_t **)SMESH_ALLOC(vec_size * sizeof(real_t *));
            real_t  **e_to   = (real_t **)SMESH_ALLOC(vec_size * sizeof(real_t *));
            u16      *weight = (u16 *)SMESH_ALLOC(nxe * sizeof(u16));

            for (int d = 0; d < vec_size; d++) {
                e_from[d] = (real_t *)SMESH_ALLOC(nxe * sizeof(real_t));
                e_to[d]   = (real_t *)SMESH_ALLOC(8 * sizeof(real_t));
            }

            idx_t *ev = (idx_t *)SMESH_ALLOC(nxe * sizeof(idx_t));

            const int corners[8] = {// Bottom
                                    sshex8_lidx(level, 0, 0, 0),
                                    sshex8_lidx(level, level, 0, 0),
                                    sshex8_lidx(level, level, level, 0),
                                    sshex8_lidx(level, 0, level, 0),
                                    // Top
                                    sshex8_lidx(level, 0, 0, level),
                                    sshex8_lidx(level, level, 0, level),
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
                                f[0] = xm * ym * zm;  // (0, 0, 0)
                                f[1] = x * ym * zm;   // (1, 0, 0)
                                f[2] = x * y * zm;    // (1, 1, 0)
                                f[3] = xm * y * zm;   // (0, 1, 0)
                                f[4] = xm * ym * z;   // (0, 0, 1)
                                f[5] = x * ym * z;    // (1, 0, 1)
                                f[6] = x * y * z;     // (1, 1, 1)
                                f[7] = xm * y * z;    // (0, 1, 1)

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
                SMESH_FREE(e_to[d]);
                SMESH_FREE(e_from[d]);
            }

            SMESH_FREE(e_from);
            SMESH_FREE(e_to);
            SMESH_FREE(ev);
            SMESH_FREE(weight);
        }

        return SMESH_SUCCESS;
    }

    // OpenMP nested parallelization? https://docs.oracle.com/cd/E19205-01/819-5270/aewbc/index.html
    // https://ppc.cs.aalto.fi/ch3/nested/

    template <typename idx_t>
    int sshex8_element_node_incidence_count(const int                                               level,
                                            const int                                               stride,
                                            const ptrdiff_t                                         nelements,
                                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                            uint16_t *const SMESH_RESTRICT                          count) {
        for (int zi = 0; zi <= level; zi++) {
            for (int yi = 0; yi <= level; yi++) {
                for (int xi = 0; xi <= level; xi++) {
                    const int v = sshex8_lidx(level * stride, xi * stride, yi * stride, zi * stride);
#pragma omp parallel for
                    for (ptrdiff_t i = 0; i < nelements; ++i) {
#pragma omp atomic update
                        count[elements[v][i]]++;
                    }
                }
            }
        }
        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename real_t>
    int sshex8_restrict(const ptrdiff_t                                         nelements,
                        const int                                               from_level,
                        const int                                               from_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                        const uint16_t *const SMESH_RESTRICT                    from_element_to_node_incidence_count,
                        const int                                               to_level,
                        const int                                               to_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                        const int                                               vec_size,
                        const real_t *const SMESH_RESTRICT                      from,
                        real_t *const SMESH_RESTRICT                            to) {
        SMESH_ASSERT(from_level % to_level == 0);

        if (from_level % to_level != 0) {
            SMESH_ERROR("Nested meshes requirement: to_level must be divisible by from_level!");
            return SMESH_FAILURE;
        }

#pragma omp parallel
        {
            const int from_nxe    = sshex8_nxe(from_level);
            const int to_nxe      = sshex8_nxe(to_level);
            const int step_factor = from_level / to_level;

            scalar_t **from_coeffs = (scalar_t **)SMESH_ALLOC(vec_size * sizeof(scalar_t *));
            for (int d = 0; d < vec_size; d++) {
                from_coeffs[d] = (scalar_t *)SMESH_ALLOC(from_nxe * sizeof(scalar_t));
            }

            scalar_t **to_coeffs = (scalar_t **)SMESH_ALLOC(vec_size * sizeof(scalar_t *));
            for (int d = 0; d < vec_size; d++) {
                to_coeffs[d] = (scalar_t *)SMESH_ALLOC(to_nxe * sizeof(scalar_t));
            }

#pragma omp for
            for (ptrdiff_t e = 0; e < nelements; e++) {
                {  // Gather elemental data
                    for (int zi = 0; zi <= from_level; zi++) {
                        for (int yi = 0; yi <= from_level; yi++) {
                            for (int xi = 0; xi <= from_level; xi++) {
                                const int v         = sshex8_lidx(from_level, xi, yi, zi);
                                const int strided_v = sshex8_lidx(from_level * from_level_stride,
                                                                  xi * from_level_stride,
                                                                  yi * from_level_stride,
                                                                  zi * from_level_stride);

                                const ptrdiff_t gid = from_elements[strided_v][e];
                                for (int d = 0; d < vec_size; d++) {
                                    from_coeffs[d][v] = from[gid * vec_size + d] / from_element_to_node_incidence_count[gid];
                                }
                            }
                        }
                    }

                    for (int d = 0; d < vec_size; d++) {
                        memset(to_coeffs[d], 0, to_nxe * sizeof(scalar_t));
                    }
                }

                for (int d = 0; d < vec_size; d++) {
                    scalar_t *in  = from_coeffs[d];
                    scalar_t *out = to_coeffs[d];

                    for (int zi = 0; zi <= from_level; zi++) {
                        for (int yi = 0; yi <= from_level; yi++) {
                            for (int xi = 0; xi <= from_level; xi++) {
                                const int v = sshex8_lidx(from_level, xi, yi, zi);
                                // Floor
                                const int cxi = xi / step_factor;
                                const int cyi = yi / step_factor;
                                const int czi = zi / step_factor;

                                const int xinc = std::min(cxi + 1, to_level) - cxi;
                                const int yinc = std::min(cyi + 1, to_level) - cyi;
                                const int zinc = std::min(czi + 1, to_level) - czi;

                                const scalar_t lx = (xi - cxi * step_factor) / ((scalar_t)step_factor);
                                const scalar_t ly = (yi - cyi * step_factor) / ((scalar_t)step_factor);
                                const scalar_t lz = (zi - czi * step_factor) / ((scalar_t)step_factor);

                                SMESH_ASSERT(lx <= 1 + 1e-8);
                                SMESH_ASSERT(lx >= -1e-8);

                                const scalar_t phi0x[2] = {1 - lx, lx};
                                const scalar_t phi0y[2] = {1 - ly, ly};
                                const scalar_t phi0z[2] = {1 - lz, lz};

                                for (int kk = 0; kk <= zinc; kk++) {
                                    for (int jj = 0; jj <= yinc; jj++) {
                                        for (int ii = 0; ii <= xinc; ii++) {
                                            const scalar_t val = phi0x[ii] * phi0y[jj] * phi0z[kk] * in[v];
                                            out[sshex8_lidx(to_level, cxi + ii, cyi + jj, czi + kk)] += val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                for (int d = 0; d < vec_size; d++) {
                    for (int zi = 0; zi <= to_level; zi++) {
                        for (int yi = 0; yi <= to_level; yi++) {
                            for (int xi = 0; xi <= to_level; xi++) {
                                // Use top level stride
                                const int to_lidx =  // sshex8_lidx(to_level, xi, yi, zi);
                                        sshex8_lidx(to_level * to_level_stride,
                                                    xi * to_level_stride,
                                                    yi * to_level_stride,
                                                    zi * to_level_stride);

                                const idx_t idx = to_elements[to_lidx][e];
#pragma omp atomic update
                                to[idx * vec_size + d] += to_coeffs[d][sshex8_lidx(to_level, xi, yi, zi)];
                            }
                        }
                    }
                }
            }

            for (int d = 0; d < vec_size; d++) {
                SMESH_FREE(from_coeffs[d]);
            }

            SMESH_FREE(from_coeffs);

            for (int d = 0; d < vec_size; d++) {
                SMESH_FREE(to_coeffs[d]);
            }

            SMESH_FREE(to_coeffs);
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_SSHEX8_RESTRICTION_IMPL_HPP
