#ifndef SMESH_SSQUAD4_RESTRICTION_IMPL_HPP
#define SMESH_SSQUAD4_RESTRICTION_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_ssquad4_restriction.hpp"

#include "smesh_sort.hpp"
#include "smesh_ssquad4.hpp"

namespace smesh {

    template <typename idx_t>
    int ssquad4_element_node_incidence_count(const int                                             level,
                                             const int                                             stride,
                                             const ptrdiff_t                                       nelements,
                                             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                             uint16_t *const SMESH_RESTRICT                         count) {
        for (int yi = 0; yi <= level; yi++) {
            for (int xi = 0; xi <= level; xi++) {
                const int v = ssquad4_lidx(level * stride, xi * stride, yi * stride);
#pragma omp parallel for
                for (ptrdiff_t i = 0; i < nelements; ++i) {
#pragma omp atomic update
                    count[elements[v][i]]++;
                }
            }
        }

        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int ssquad4_hierarchical_restriction(int                                                   level,
                                         const ptrdiff_t                                       nelements,
                                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         const uint16_t *const SMESH_RESTRICT                   element_to_node_incidence_count,
                                         const int                                             vec_size,
                                         const T *const SMESH_RESTRICT                          from,
                                         T *const SMESH_RESTRICT                                to) {
#pragma omp parallel
        {
            const int nxe    = ssquad4_nxe(level);
            T       **e_from = (T **)malloc(vec_size * sizeof(T *));
            T       **e_to   = (T **)malloc(vec_size * sizeof(T *));
            uint16_t *weight = (uint16_t *)malloc(nxe * sizeof(uint16_t));

            for (int d = 0; d < vec_size; d++) {
                e_from[d] = (T *)malloc(nxe * sizeof(T));
                e_to[d]   = (T *)malloc(4 * sizeof(T));
            }

            idx_t *ev = (idx_t *)malloc(nxe * sizeof(idx_t));

            const int corners[4] = {ssquad4_lidx(level, 0, 0),
                                    ssquad4_lidx(level, level, 0),
                                    ssquad4_lidx(level, level, level),
                                    ssquad4_lidx(level, 0, level)};

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
                            assert(e_from[d][v] == e_from[d][v]);
                        }
                    }

                    for (int v = 0; v < nxe; v++) {
                        weight[v] = element_to_node_incidence_count[ev[v]];
                    }

                    for (int d = 0; d < vec_size; d++) {
                        for (int v = 0; v < 4; v++) {
                            e_to[d][v] = 0;
                        }
                    }

                    const T h = 1. / level;
                    // Iterate over structrued grid (nodes)

                    for (int yi = 0; yi < level + 1; yi++) {
                        for (int xi = 0; xi < level + 1; xi++) {
                            const int lidx = ssquad4_lidx(level, xi, yi);

                            const T x = xi * h;
                            const T y = yi * h;

                            // Evaluate Quad4 basis functions at x, y
                            const T xm = (1 - x);
                            const T ym = (1 - y);

                            T f[4];
                            f[0] = xm * ym;  // (0, 0, 0)
                            f[1] = x * ym;   // (1, 0, 0)
                            f[2] = x * y;    // (1, 1, 0)
                            f[3] = xm * y;   // (0, 1, 0)

                            for (int d = 0; d < vec_size; d++) {
                                const T val = e_from[d][lidx] / weight[lidx];
                                for (int i = 0; i < 4; i++) {
                                    e_to[d][i] += f[i] * val;
                                }
                            }
                        }
                    }

                    for (int d = 0; d < vec_size; d++) {
                        for (int i = 0; i < 4; i++) {
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

    template <typename idx_t, typename T>
    int ssquad4_restrict(const ptrdiff_t                     nelements,
                         const int                           from_level,
                         const int                           from_level_stride,
                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                         const uint16_t *const SMESH_RESTRICT from_element_to_node_incidence_count,
                         const int                           to_level,
                         const int                           to_level_stride,
                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                         const int                           vec_size,
                         const T *const SMESH_RESTRICT   from,
                         T *const SMESH_RESTRICT         to) {
        // FIXME this should be handled outside!!!
        // if (to_level == 1) {
        //     return ssquad4_hierarchical_restriction(
        //             from_level, nelements, from_elements, from_element_to_node_incidence_count, vec_size, from, to);
        // }

        assert(from_level % to_level == 0);

        if (from_level % to_level != 0) {
            SMESH_ERROR("Nested meshes requirement: from_level must be divisible by to_level!");
            return SMESH_FAILURE;
        }

#pragma omp parallel
        {
            const int from_nxe    = ssquad4_nxe(from_level);
            const int step_factor = from_level / to_level;

            T **from_coeffs = (T **)malloc(vec_size * sizeof(T *));
            for (int d = 0; d < vec_size; d++) {
                from_coeffs[d] = (T *)malloc(from_nxe * sizeof(T));
            }

#pragma omp for
            for (ptrdiff_t e = 0; e < nelements; e++) {
                {  // Gather elemental coefficients
                    for (int yi = 0; yi <= from_level; yi++) {
                        for (int xi = 0; xi <= from_level; xi++) {
                            const int v = ssquad4_lidx(from_level, xi, yi);
                            const int strided_v =
                                    ssquad4_lidx(from_level * from_level_stride, xi * from_level_stride, yi * from_level_stride);
                            const ptrdiff_t gid = from_elements[strided_v][e];

                            for (int d = 0; d < vec_size; d++) {
                                assert(from_element_to_node_incidence_count[gid] != 0);
                                from_coeffs[d][v] = from[gid * vec_size + d] / from_element_to_node_incidence_count[gid];
                            }
                        }
                    }
                }

                const T h = to_level * 1. / from_level;
                for (int d = 0; d < vec_size; d++) {
                    T *c = from_coeffs[d];

                    // Restict the coefficients along the x-axis (edges)
                    for (int yi = 0; yi <= to_level; yi++) {
                        for (int xi = 0; xi < to_level; xi++) {
                            for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                                const T   fl          = (1 - between_xi * h);
                                const T   fr          = (between_xi * h);
                                const int between_idx = ssquad4_lidx(from_level, xi * step_factor + between_xi, yi * step_factor);
                                c[ssquad4_lidx(from_level, xi * step_factor, yi * step_factor)] += fl * c[between_idx];
                                c[ssquad4_lidx(from_level, (xi + 1) * step_factor, yi * step_factor)] += fr * c[between_idx];
                            }
                        }
                    }

                    // Restrict the coefficients along the x-axis (center)
                    for (int yi = 0; yi < to_level; yi++) {
                        for (int xi = 0; xi < to_level; xi++) {
                            for (int between_yi = 1; between_yi < step_factor; between_yi++) {
                                const int yy    = yi * step_factor + between_yi;
                                const int left  = ssquad4_lidx(from_level, xi * step_factor, yy);
                                const int right = ssquad4_lidx(from_level, (xi + 1) * step_factor, yy);

                                for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                                    const T fl = (1 - between_xi * h);
                                    const T fr = (between_xi * h);

                                    const int xx     = xi * step_factor + between_xi;
                                    const int center = ssquad4_lidx(from_level, xx, yy);

                                    c[left] += fl * c[center];
                                    c[right] += fr * c[center];
                                }
                            }
                        }
                    }

                    // Restrict the coefficients along the y-axis (edges)
                    for (int yi = 0; yi < to_level; yi++) {
                        for (int xi = 0; xi <= to_level; xi++) {
                            for (int between_yi = 1; between_yi < step_factor; between_yi++) {
                                const T   fb      = (1 - between_yi * h);
                                const T   ft      = (between_yi * h);
                                const int to_lidx = ssquad4_lidx(from_level, xi * step_factor, yi * step_factor + between_yi);
                                c[ssquad4_lidx(from_level, xi * step_factor, yi * step_factor)] += fb * c[to_lidx];
                                c[ssquad4_lidx(from_level, xi * step_factor, (yi + 1) * step_factor)] += ft * c[to_lidx];
                            }
                        }
                    }
                }

                // Scatter elemental data
                // Extract coarse coeffs and discard rest
                for (int d = 0; d < vec_size; d++) {
                    for (int yi = 0; yi <= to_level; yi++) {
                        for (int xi = 0; xi <= to_level; xi++) {
                            // Use top level stride
                            const int to_lidx =
                                    ssquad4_lidx(to_level * to_level_stride, xi * to_level_stride, yi * to_level_stride);

                            // Use stride to convert from "from" to "to" local indexing
                            const int from_lidx = ssquad4_lidx(from_level, xi * step_factor, yi * step_factor);

                            const idx_t idx = to_elements[to_lidx][e];
                            assert(from_coeffs[d][from_lidx] == from_coeffs[d][from_lidx]);
#pragma omp atomic update
                            to[idx * vec_size + d] += from_coeffs[d][from_lidx];
                        }
                    }
                }
            }

            for (int d = 0; d < vec_size; d++) {
                free(from_coeffs[d]);
            }

            free(from_coeffs);
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif // SMESH_SSQUAD4_RESTRICTION_IMPL_HPP
