#ifndef SMESH_SSTET4_RESTRICTION_IMPL_HPP
#define SMESH_SSTET4_RESTRICTION_IMPL_HPP

#include "smesh_sstet4_restriction.hpp"
#include "smesh_alloc.hpp"
#include "smesh_sstet4_transfer.impl.hpp"

namespace smesh {

    template <typename idx_t>
    int sstet4_element_node_incidence_count(const int                                               level,
                                            const int                                               stride,
                                            const ptrdiff_t                                         nelements,
                                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                            uint16_t *const SMESH_RESTRICT                          count) {
        for (int z = 0; z <= level; ++z) {
            for (int y = 0; y <= level - z; ++y) {
                for (int x = 0; x <= level - z - y; ++x) {
                    const int v = sstet4_lidx(level * stride, x * stride, y * stride, z * stride);
#pragma omp parallel for schedule(static)
                    for (ptrdiff_t e = 0; e < nelements; ++e) {
#pragma omp atomic update
                        count[elements[v][e]]++;
                    }
                }
            }
        }

        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int sstet4_hierarchical_restriction(int                                                     level,
                                        const ptrdiff_t                                         nelements,
                                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                        const uint16_t *const SMESH_RESTRICT                    element_to_node_incidence_count,
                                        const int                                               vec_size,
                                        const T *const SMESH_RESTRICT                           from,
                                        T *const SMESH_RESTRICT                                 to) {
        const int corners[4] = {sstet4_lidx(level, 0, 0, 0),
                                sstet4_lidx(level, level, 0, 0),
                                sstet4_lidx(level, 0, level, 0),
                                sstet4_lidx(level, 0, 0, level)};

#pragma omp parallel for schedule(static)
        for (ptrdiff_t e = 0; e < nelements; ++e) {
            for (int z = 0; z <= level; ++z) {
                for (int y = 0; y <= level - z; ++y) {
                    for (int x = 0; x <= level - z - y; ++x) {
                        const int   lidx = sstet4_lidx(level, x, y, z);
                        const T     h    = T(1) / T(level);
                        const T     w[4] = {T(1) - T(x + y + z) * h, T(x) * h, T(y) * h, T(z) * h};
                        const idx_t gid  = elements[lidx][e];
                        const T     inv_count = T(1) / T(element_to_node_incidence_count[gid]);

                        for (int d = 0; d < vec_size; ++d) {
                            const T val = from[gid * vec_size + d] * inv_count;
                            for (int v = 0; v < 4; ++v) {
                                const idx_t coarse_gid = elements[corners[v]][e];
#pragma omp atomic update
                                to[coarse_gid * vec_size + d] += w[v] * val;
                            }
                        }
                    }
                }
            }
        }

        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int sstet4_restrict(const ptrdiff_t                                         nelements,
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
        SMESH_ASSERT(from_level % to_level == 0);
        if (from_level % to_level != 0) {
            SMESH_ERROR("Nested meshes requirement: from_level must be divisible by to_level!");
            return SMESH_FAILURE;
        }

        const int to_nxe = sstet4_nxe(to_level);
        int       invalid = 0;

#pragma omp parallel reduction(| : invalid)
        {
            int *const cx = (int *)SMESH_ALLOC(static_cast<size_t>(to_nxe) * sizeof(int));
            int *const cy = (int *)SMESH_ALLOC(static_cast<size_t>(to_nxe) * sizeof(int));
            int *const cz = (int *)SMESH_ALLOC(static_cast<size_t>(to_nxe) * sizeof(int));
            sstet4_transfer::fill_lattice_coords(to_level, cx, cy, cz);

#pragma omp for schedule(static)
            for (ptrdiff_t e = 0; e < nelements; ++e) {
                for (int z = 0; z <= from_level; ++z) {
                    for (int y = 0; y <= from_level - z; ++y) {
                        for (int x = 0; x <= from_level - z - y; ++x) {
                            int coarse_lidx[4];
                            T   w[4];
                            const int err = sstet4_transfer::find_weights(to_level,
                                                                          T(x) * T(to_level) / T(from_level),
                                                                          T(y) * T(to_level) / T(from_level),
                                                                          T(z) * T(to_level) / T(from_level),
                                                                          cx,
                                                                          cy,
                                                                          cz,
                                                                          coarse_lidx,
                                                                          w);
                            if (err != SMESH_SUCCESS) {
                                invalid = 1;
                                continue;
                            }

                            const int from_lidx = sstet4_lidx(from_level * from_level_stride,
                                                               x * from_level_stride,
                                                               y * from_level_stride,
                                                               z * from_level_stride);
                            const idx_t from_gid = from_elements[from_lidx][e];
                            const T     inv_count = T(1) / T(from_element_to_node_incidence_count[from_gid]);

                            for (int d = 0; d < vec_size; ++d) {
                                const T val = from[from_gid * vec_size + d] * inv_count;
                                for (int v = 0; v < 4; ++v) {
                                    const int to_lidx = sstet4_lidx(to_level * to_level_stride,
                                                                    cx[coarse_lidx[v]] * to_level_stride,
                                                                    cy[coarse_lidx[v]] * to_level_stride,
                                                                    cz[coarse_lidx[v]] * to_level_stride);
                                    const idx_t to_gid = to_elements[to_lidx][e];
#pragma omp atomic update
                                    to[to_gid * vec_size + d] += w[v] * val;
                                }
                            }
                        }
                    }
                }
            }

            SMESH_FREE(cx);
            SMESH_FREE(cy);
            SMESH_FREE(cz);
        }

        return invalid ? SMESH_FAILURE : SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_SSTET4_RESTRICTION_IMPL_HPP
