#ifndef SMESH_SSHEX8_PROLONGATION_IMPL_HPP
#define SMESH_SSHEX8_PROLONGATION_IMPL_HPP

#include "smesh_sshex8_prolongation.hpp"
#include "smesh_sshex8.hpp"

#include <math.h>
#include <stdio.h>
#include <string.h>

namespace smesh {

    template <typename idx_t, typename T>
    int sshex8_hierarchical_prolongation(int                                                     level,
                                         const ptrdiff_t                                         nelements,
                                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         const int                                               vec_size,
                                         const T *const SMESH_RESTRICT                           from,
                                         T *const SMESH_RESTRICT                                 to) {
#pragma omp parallel
        {
            const int nxe    = sshex8_nxe(level);
            T       **e_from = (T **)malloc(vec_size * sizeof(T *));
            T       **e_to   = (T **)malloc(vec_size * sizeof(T *));

            for (int d = 0; d < vec_size; d++) {
                e_from[d] = (T *)malloc(8 * sizeof(T));
                e_to[d]   = (T *)malloc(nxe * sizeof(T));
            }

            idx_t *ev = (idx_t *)malloc(nxe * sizeof(idx_t));

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
                            e_to[d][v] = 0;
                        }
                    }

                    for (int d = 0; d < vec_size; d++) {
                        for (int v = 0; v < 8; v++) {
                            e_from[d][v] = from[ev[corners[v]] * vec_size + d];
                            SMESH_ASSERT(e_from[d][v] == e_from[d][v]);
                        }
                    }

                    const T h = 1. / level;
                    // Iterate over structrued grid (nodes)
                    for (int zi = 0; zi < level + 1; zi++) {
                        for (int yi = 0; yi < level + 1; yi++) {
                            for (int xi = 0; xi < level + 1; xi++) {
                                int lidx = sshex8_lidx(level, xi, yi, zi);

                                const T x = xi * h;
                                const T y = yi * h;
                                const T z = zi * h;

                                // Evaluate Hex8 basis functions at x, y, z
                                const T xm = (1 - x);
                                const T ym = (1 - y);
                                const T zm = (1 - z);

                                T f[8];
                                f[0] = xm * ym * zm;  // (0, 0, 0)
                                f[1] = x * ym * zm;   // (1, 0, 0)
                                f[2] = x * y * zm;    // (1, 1, 0)
                                f[3] = xm * y * zm;   // (0, 1, 0)
                                f[4] = xm * ym * z;   // (0, 0, 1)
                                f[5] = x * ym * z;    // (1, 0, 1)
                                f[6] = x * y * z;     // (1, 1, 1)
                                f[7] = xm * y * z;    // (0, 1, 1)

                                for (int d = 0; d < vec_size; d++) {
                                    for (int v = 0; v < 8; v++) {
                                        e_to[d][lidx] += f[v] * e_from[d][v];
                                    }
                                }
                            }
                        }
                    }

                    for (int i = 0; i < nxe; i++) {
                        for (int d = 0; d < vec_size; d++) {
                            to[ev[i] * vec_size + d] = e_to[d][i];
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
        }

        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int sshex8_prolongate(const ptrdiff_t                                         nelements,
                          const int                                               from_level,
                          const int                                               from_level_stride,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                          const int                                               to_level,
                          const int                                               to_level_stride,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                          const int                                               vec_size,
                          const T *const SMESH_RESTRICT                           from,
                          T *const SMESH_RESTRICT                                 to) {
        SMESH_ASSERT(to_level % from_level == 0);

        if (to_level % from_level != 0) {
            SMESH_ERROR("Nested meshes requirement: to_level must be divisible by from_level!");
            return SMESH_FAILURE;
        }

#pragma omp parallel
        {
            const int from_nxe    = sshex8_nxe(from_level);
            const int to_nxe      = sshex8_nxe(to_level);
            const int step_factor = to_level / from_level;

            T **to_coeffs = (T **)malloc(vec_size * sizeof(T *));
            for (int d = 0; d < vec_size; d++) {
                to_coeffs[d] = (T *)malloc(to_nxe * sizeof(T));
            }

#pragma omp for
            for (ptrdiff_t e = 0; e < nelements; e++) {
                {
#ifndef NDEBUG
                    // Only for debugging
                    for (int d = 0; d < vec_size; d++) {
                        memset(to_coeffs[d], 0, to_nxe * sizeof(T));
                    }
#endif
                    // Gather elemental data
                    // Fill matching nodes with from data while gathering
                    for (int d = 0; d < vec_size; d++) {
                        for (int zi = 0; zi <= from_level; zi++) {
                            for (int yi = 0; yi <= from_level; yi++) {
                                for (int xi = 0; xi <= from_level; xi++) {
                                    // Use top level stride
                                    const int from_lidx = sshex8_lidx(from_level * from_level_stride,
                                                                      xi * from_level_stride,
                                                                      yi * from_level_stride,
                                                                      zi * from_level_stride);

                                    // Use stride to convert from "from" to "to" local indexing
                                    const int to_lidx =
                                            sshex8_lidx(to_level, xi * step_factor, yi * step_factor, zi * step_factor);

                                    const idx_t idx = from_elements[from_lidx][e];
                                    const T     val = from[idx * vec_size + d];

                                    to_coeffs[d][to_lidx] = val;
                                    SMESH_ASSERT(to_coeffs[d][to_lidx] == to_coeffs[d][to_lidx]);
                                }
                            }
                        }
                    }
                }

                const T to_h = from_level * 1. / to_level;
                for (int d = 0; d < vec_size; d++) {
                    T *c = to_coeffs[d];

                    // Interpolate the coefficients along the x-axis (edges)
                    for (int zi = 0; zi <= from_level; zi++) {
                        for (int yi = 0; yi <= from_level; yi++) {
                            for (int xi = 0; xi < from_level; xi++) {
                                const T c0 = c[sshex8_lidx(to_level, xi * step_factor, yi * step_factor, zi * step_factor)];
                                const T c1 = c[sshex8_lidx(to_level, (xi + 1) * step_factor, yi * step_factor, zi * step_factor)];

                                for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                                    const T   fl      = (1 - between_xi * to_h);
                                    const T   fr      = (between_xi * to_h);
                                    const int to_lidx = sshex8_lidx(
                                            to_level, xi * step_factor + between_xi, yi * step_factor, zi * step_factor);
                                    c[to_lidx] = fl * c0 + fr * c1;
                                }
                            }
                        }
                    }

                    // Interpolate the coefficients along the y-axis (edges)
                    for (int zi = 0; zi <= from_level; zi++) {
                        for (int yi = 0; yi < from_level; yi++) {
                            for (int xi = 0; xi <= from_level; xi++) {
                                const T c0 = c[sshex8_lidx(to_level, xi * step_factor, yi * step_factor, zi * step_factor)];
                                const T c1 = c[sshex8_lidx(to_level, xi * step_factor, (yi + 1) * step_factor, zi * step_factor)];

                                for (int between_yi = 1; between_yi < step_factor; between_yi++) {
                                    const T   fb      = (1 - between_yi * to_h);
                                    const T   ft      = (between_yi * to_h);
                                    const int to_lidx = sshex8_lidx(
                                            to_level, xi * step_factor, yi * step_factor + between_yi, zi * step_factor);
                                    c[to_lidx] = fb * c0 + ft * c1;
                                }
                            }
                        }
                    }

                    // Interpolate the coefficients along the z-axis (edges)
                    for (int zi = 0; zi < from_level; zi++) {
                        for (int yi = 0; yi <= from_level; yi++) {
                            for (int xi = 0; xi <= from_level; xi++) {
                                const T c0 = c[sshex8_lidx(to_level, xi * step_factor, yi * step_factor, zi * step_factor)];
                                const T c1 = c[sshex8_lidx(to_level, xi * step_factor, yi * step_factor, (zi + 1) * step_factor)];

                                for (int between_zi = 1; between_zi < step_factor; between_zi++) {
                                    const T   fb      = (1 - between_zi * to_h);
                                    const T   ft      = (between_zi * to_h);
                                    const int to_lidx = sshex8_lidx(
                                            to_level, xi * step_factor, yi * step_factor, zi * step_factor + between_zi);
                                    c[to_lidx] = fb * c0 + ft * c1;
                                }
                            }
                        }
                    }

                    // Interpolate the coefficients along the x-axis (center) in the x-y-planes
                    for (int zi = 0; zi <= from_level; zi++) {
                        for (int yi = 0; yi < from_level; yi++) {
                            for (int between_yi = 1; between_yi < step_factor; between_yi++) {
                                for (int xi = 0; xi < from_level; xi++) {
                                    const int zz    = zi * step_factor;
                                    const int yy    = yi * step_factor + between_yi;
                                    const int left  = sshex8_lidx(to_level, xi * step_factor, yy, zz);
                                    const int right = sshex8_lidx(to_level, (xi + 1) * step_factor, yy, zz);
                                    const T   cl    = c[left];
                                    const T   cr    = c[right];

                                    for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                                        const T fl = (1 - between_xi * to_h);
                                        const T fr = (between_xi * to_h);

                                        const int xx     = xi * step_factor + between_xi;
                                        const int center = sshex8_lidx(to_level, xx, yy, zz);
                                        c[center]        = fl * cl + fr * cr;
                                    }
                                }
                            }
                        }
                    }

                    // Interpolate the coefficients along the x-axis (center) in the x-z-planes
                    for (int zi = 0; zi < from_level; zi++) {
                        for (int between_zi = 1; between_zi < step_factor; between_zi++) {
                            for (int yi = 0; yi <= from_level; yi++) {
                                for (int xi = 0; xi < from_level; xi++) {
                                    const int yy    = yi * step_factor;
                                    const int zz    = zi * step_factor + between_zi;
                                    const int left  = sshex8_lidx(to_level, xi * step_factor, yy, zz);
                                    const int right = sshex8_lidx(to_level, (xi + 1) * step_factor, yy, zz);
                                    const T   cl    = c[left];
                                    const T   cr    = c[right];

                                    for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                                        const T fl = (1 - between_xi * to_h);
                                        const T fr = (between_xi * to_h);

                                        const int xx     = xi * step_factor + between_xi;
                                        const int center = sshex8_lidx(to_level, xx, yy, zz);
                                        c[center]        = fl * cl + fr * cr;
                                    }
                                }
                            }
                        }
                    }

                    // Interpolate the coefficients along the y-axis (center) in the y-z-planes
                    for (int zi = 0; zi < from_level; zi++) {
                        for (int between_zi = 1; between_zi < step_factor; between_zi++) {
                            for (int yi = 0; yi < from_level; yi++) {
                                for (int xi = 0; xi <= from_level; xi++) {
                                    const int xx    = xi * step_factor;
                                    const int zz    = zi * step_factor + between_zi;
                                    const int left  = sshex8_lidx(to_level, xx, yi * step_factor, zz);
                                    const int right = sshex8_lidx(to_level, xx, (yi + 1) * step_factor, zz);
                                    const T   cl    = c[left];
                                    const T   cr    = c[right];

                                    for (int between_yi = 1; between_yi < step_factor; between_yi++) {
                                        const T fl = (1 - between_yi * to_h);
                                        const T fr = (between_yi * to_h);

                                        const int yy     = yi * step_factor + between_yi;
                                        const int center = sshex8_lidx(to_level, xx, yy, zz);
                                        c[center]        = fl * cl + fr * cr;
                                    }
                                }
                            }
                        }
                    }

                    // Interpolate the coefficients along the x-axis (center) in the x-y-planes
                    for (int zi = 0; zi < from_level; zi++) {
                        for (int between_zi = 1; between_zi < step_factor; between_zi++) {
                            for (int yi = 0; yi < from_level; yi++) {
                                for (int between_yi = 1; between_yi < step_factor; between_yi++) {
                                    for (int xi = 0; xi < from_level; xi++) {
                                        const int zz    = zi * step_factor + between_zi;
                                        const int yy    = yi * step_factor + between_yi;
                                        const int left  = sshex8_lidx(to_level, xi * step_factor, yy, zz);
                                        const int right = sshex8_lidx(to_level, (xi + 1) * step_factor, yy, zz);
                                        const T   cl    = c[left];
                                        const T   cr    = c[right];

                                        for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                                            const T fl = (1 - between_xi * to_h);
                                            const T fr = (between_xi * to_h);

                                            const int xx     = xi * step_factor + between_xi;
                                            const int center = sshex8_lidx(to_level, xx, yy, zz);
                                            c[center]        = fl * cl + fr * cr;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Scatter elemental data
                for (int zi = 0; zi <= to_level; zi++) {
                    for (int yi = 0; yi <= to_level; yi++) {
                        for (int xi = 0; xi <= to_level; xi++) {
                            const int v         = sshex8_lidx(to_level, xi, yi, zi);
                            const int strided_v = sshex8_lidx(
                                    to_level * to_level_stride, xi * to_level_stride, yi * to_level_stride, zi * to_level_stride);
                            const ptrdiff_t gid = to_elements[strided_v][e];

                            for (int d = 0; d < vec_size; d++) {
                                to[gid * vec_size + d] = to_coeffs[d][v];
                            }
                        }
                    }
                }
            }

            for (int d = 0; d < vec_size; d++) {
                free(to_coeffs[d]);
            }

            free(to_coeffs);
        }

        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int sshex8_restrict(const ptrdiff_t                      nelements,
                        const int                            from_level,
                        const int                            from_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT         from_elements,
                        const uint16_t *const SMESH_RESTRICT from_element_to_node_incidence_count,
                        const int                            to_level,
                        const int                            to_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT         to_elements,
                        const int                            vec_size,
                        const T *const SMESH_RESTRICT        from,
                        T *const SMESH_RESTRICT              to) {
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

            T **from_coeffs = (T **)malloc(vec_size * sizeof(T *));
            for (int d = 0; d < vec_size; d++) {
                from_coeffs[d] = (T *)malloc(from_nxe * sizeof(T));
            }

            T **to_coeffs = (T **)malloc(vec_size * sizeof(T *));
            for (int d = 0; d < vec_size; d++) {
                to_coeffs[d] = (T *)malloc(to_nxe * sizeof(T));
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
                        memset(to_coeffs[d], 0, to_nxe * sizeof(T));
                    }
                }

                for (int d = 0; d < vec_size; d++) {
                    T *in  = from_coeffs[d];
                    T *out = to_coeffs[d];

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

                                const T lx = (xi - cxi * step_factor) / ((T)step_factor);
                                const T ly = (yi - cyi * step_factor) / ((T)step_factor);
                                const T lz = (zi - czi * step_factor) / ((T)step_factor);

                                SMESH_ASSERT(lx <= 1 + 1e-8);
                                SMESH_ASSERT(lx >= -1e-8);

                                const T phi0x[2] = {1 - lx, lx};
                                const T phi0y[2] = {1 - ly, ly};
                                const T phi0z[2] = {1 - lz, lz};

                                for (int kk = 0; kk <= zinc; kk++) {
                                    for (int jj = 0; jj <= yinc; jj++) {
                                        for (int ii = 0; ii <= xinc; ii++) {
                                            const T val = phi0x[ii] * phi0y[jj] * phi0z[kk] * in[v];
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
                free(from_coeffs[d]);
            }

            free(from_coeffs);

            for (int d = 0; d < vec_size; d++) {
                free(to_coeffs[d]);
            }

            free(to_coeffs);
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_SSHEX8_PROLONGATION_IMPL_HPP
