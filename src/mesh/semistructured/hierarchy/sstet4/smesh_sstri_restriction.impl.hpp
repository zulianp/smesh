#ifndef SMESH_SSTRI_RESTRICTION_IMPL_HPP
#define SMESH_SSTRI_RESTRICTION_IMPL_HPP

#include "smesh_sstri_restriction.hpp"
#include "smesh_alloc.hpp"
#include "smesh_sstet4.hpp"

namespace smesh {

    template <typename idx_t>
    int sstri_element_node_incidence_count(const int                                               level,
                                           const int                                               stride,
                                           const ptrdiff_t                                         nelements,
                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                           uint16_t *const SMESH_RESTRICT                          count) {
        for (int y = 0; y <= level; ++y) {
            for (int x = 0; x <= level - y; ++x) {
                const int v = sstri_lidx(level * stride, x * stride, y * stride);
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
    int sstri_restrict(const ptrdiff_t                                         nelements,
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

        const int step = from_level / to_level;
        const T   h    = T(1) / T(step);

#pragma omp parallel for
        for (ptrdiff_t e = 0; e < nelements; ++e) {
            for (int y = 0; y <= from_level; ++y) {
                for (int x = 0; x <= from_level - y; ++x) {
                    const int from_lidx =
                            sstri_lidx(from_level * from_level_stride, x * from_level_stride, y * from_level_stride);
                    const idx_t from_gid  = from_elements[from_lidx][e];
                    const T     inv_count = T(1) / T(from_element_to_node_incidence_count[from_gid]);

                    const int X  = x / step;
                    const int Y  = y / step;
                    const int lx = x - X * step;
                    const int ly = y - Y * step;

                    int coarse_xy[3][2];
                    T   w[3];
                    int n_w = 0;

                    if (lx == 0 && ly == 0) {
                        coarse_xy[0][0] = X;
                        coarse_xy[0][1] = Y;
                        w[0]            = T(1);
                        n_w             = 1;
                    } else if (lx + ly <= step) {
                        coarse_xy[0][0] = X;
                        coarse_xy[0][1] = Y;
                        coarse_xy[1][0] = X + 1;
                        coarse_xy[1][1] = Y;
                        coarse_xy[2][0] = X;
                        coarse_xy[2][1] = Y + 1;
                        w[0]            = T(1) - T(lx) * h - T(ly) * h;
                        w[1]            = T(lx) * h;
                        w[2]            = T(ly) * h;
                        n_w             = 3;
                    } else {
                        coarse_xy[0][0] = X + 1;
                        coarse_xy[0][1] = Y;
                        coarse_xy[1][0] = X + 1;
                        coarse_xy[1][1] = Y + 1;
                        coarse_xy[2][0] = X;
                        coarse_xy[2][1] = Y + 1;
                        w[0]            = T(1) - T(ly) * h;
                        w[1]            = T(lx) * h + T(ly) * h - T(1);
                        w[2]            = T(1) - T(lx) * h;
                        n_w             = 3;
                    }

                    for (int d = 0; d < vec_size; ++d) {
                        const T val = from[from_gid * vec_size + d] * inv_count;
                        for (int v = 0; v < n_w; ++v) {
                            const int to_lidx = sstri_lidx(to_level * to_level_stride,
                                                           coarse_xy[v][0] * to_level_stride,
                                                           coarse_xy[v][1] * to_level_stride);
                            const idx_t to_gid = to_elements[to_lidx][e];
#pragma omp atomic update
                            to[to_gid * vec_size + d] += w[v] * val;
                        }
                    }
                }
            }
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_SSTRI_RESTRICTION_IMPL_HPP
