#ifndef SMESH_SSEDGE_RESTRICTION_IMPL_HPP
#define SMESH_SSEDGE_RESTRICTION_IMPL_HPP

#include "smesh_ssedge_restriction.hpp"
#include "smesh_alloc.hpp"
#include "smesh_ssquad4.hpp"

namespace smesh {

    template <typename idx_t>
    int ssedge_element_node_incidence_count(const int                                               level,
                                            const int                                               stride,
                                            const ptrdiff_t                                         nelements,
                                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                            uint16_t *const SMESH_RESTRICT                          count) {
        for (int x = 0; x <= level; ++x) {
            const int v = ssedge_lidx(level * stride, x * stride);
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < nelements; ++i) {
#pragma omp atomic update
                count[elements[v][i]]++;
            }
        }

        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int ssedge_restrict(const ptrdiff_t                                         nelements,
                        const int                                               from_level,
                        const int                                               from_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                        const uint16_t *const SMESH_RESTRICT                    from_element_to_node_incidence_count,
                        const int                                               to_level,
                        const int                                               to_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                        const int                                               vec_size,
                        const T *const SMESH_RESTRICT                           from,
                        T *const SMESH_RESTRICT                                 to) {
        if (from_level % to_level != 0) {
            SMESH_ERROR("Nested meshes requirement: from_level must be divisible by to_level!");
            return SMESH_FAILURE;
        }

#pragma omp parallel
        {
            const int from_nxe    = ssedge_nxe(from_level);
            const int step_factor = from_level / to_level;

            T **from_coeffs = (T **)SMESH_ALLOC(vec_size * sizeof(T *));
            for (int d = 0; d < vec_size; d++) {
                from_coeffs[d] = (T *)SMESH_ALLOC(from_nxe * sizeof(T));
            }

#pragma omp for
            for (ptrdiff_t e = 0; e < nelements; e++) {
                for (int xi = 0; xi <= from_level; xi++) {
                    const int v =
                            ssedge_lidx(from_level * from_level_stride, xi * from_level_stride);
                    const ptrdiff_t gid = from_elements[v][e];
                    for (int d = 0; d < vec_size; d++) {
                        from_coeffs[d][xi] = from[gid * vec_size + d] / from_element_to_node_incidence_count[gid];
                    }
                }

                const T h = to_level * T(1) / T(from_level);
                for (int d = 0; d < vec_size; d++) {
                    T *c = from_coeffs[d];
                    for (int xi = 0; xi < to_level; xi++) {
                        for (int between_xi = 1; between_xi < step_factor; between_xi++) {
                            const T fl          = (T(1) - between_xi * h);
                            const T fr          = (between_xi * h);
                            const int between_idx = xi * step_factor + between_xi;
                            c[xi * step_factor] += fl * c[between_idx];
                            c[(xi + 1) * step_factor] += fr * c[between_idx];
                        }
                    }
                }

                for (int d = 0; d < vec_size; d++) {
                    for (int xi = 0; xi <= to_level; xi++) {
                        const int to_lidx   = ssedge_lidx(to_level * to_level_stride, xi * to_level_stride);
                        const int from_lidx = xi * step_factor;
                        const idx_t idx     = to_elements[to_lidx][e];
#pragma omp atomic update
                        to[idx * vec_size + d] += from_coeffs[d][from_lidx];
                    }
                }
            }

            for (int d = 0; d < vec_size; d++) {
                SMESH_FREE(from_coeffs[d]);
            }
            SMESH_FREE(from_coeffs);
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_SSEDGE_RESTRICTION_IMPL_HPP
