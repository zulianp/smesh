#ifndef SMESH_CURVATURE_IMPL_HPP
#define SMESH_CURVATURE_IMPL_HPP

#include "smesh_curvature.hpp"

#include "smesh_alloc.hpp"
#include "smesh_common.hpp"
#include "smesh_elem_type.hpp"

#include <cmath>
#include <stdio.h>
#include <string.h>

namespace smesh {

namespace {

template <typename geom_t>
SMESH_INLINE void cross3_raw(geom_t ax,
                             geom_t ay,
                             geom_t az,
                             geom_t bx,
                             geom_t by,
                             geom_t bz,
                             geom_t *cx,
                             geom_t *cy,
                             geom_t *cz) {
    *cx = ay * bz - az * by;
    *cy = az * bx - ax * bz;
    *cz = ax * by - ay * bx;
}

template <typename idx_t, typename geom_t>
void accumulate_triangle(const idx_t                                              i0,
                         const idx_t                                              i1,
                         const idx_t                                              i2,
                         const int                                                sdim,
                         const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                         geom_t *const SMESH_RESTRICT                             kx,
                         geom_t *const SMESH_RESTRICT                             ky,
                         geom_t *const SMESH_RESTRICT                             kz,
                         geom_t *const SMESH_RESTRICT                             area) {
    const geom_t x0 = points[0][i0], y0 = points[1][i0];
    const geom_t x1 = points[0][i1], y1 = points[1][i1];
    const geom_t x2 = points[0][i2], y2 = points[1][i2];
    const geom_t z0 = sdim >= 3 ? points[2][i0] : static_cast<geom_t>(0);
    const geom_t z1 = sdim >= 3 ? points[2][i1] : static_cast<geom_t>(0);
    const geom_t z2 = sdim >= 3 ? points[2][i2] : static_cast<geom_t>(0);

    const geom_t u01x = x1 - x0, u01y = y1 - y0, u01z = z1 - z0;
    const geom_t u02x = x2 - x0, u02y = y2 - y0, u02z = z2 - z0;
    const geom_t u12x = x2 - x1, u12y = y2 - y1, u12z = z2 - z1;

    geom_t nx, ny, nz;
    cross3_raw(u01x, u01y, u01z, u02x, u02y, u02z, &nx, &ny, &nz);
    const geom_t twice_a = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (!(twice_a > static_cast<geom_t>(0))) {
        return;
    }
    const geom_t a = static_cast<geom_t>(0.5) * twice_a;
    const geom_t a3 = a / static_cast<geom_t>(3);

    auto cot_at = [&](geom_t ux, geom_t uy, geom_t uz, geom_t vx, geom_t vy, geom_t vz) -> geom_t {
        geom_t cx, cy, cz;
        cross3_raw(ux, uy, uz, vx, vy, vz, &cx, &cy, &cz);
        const geom_t den = std::sqrt(cx * cx + cy * cy + cz * cz);
        if (!(den > static_cast<geom_t>(0))) {
            return static_cast<geom_t>(0);
        }
        return (ux * vx + uy * vy + uz * vz) / den;
    };

    const geom_t cot0 = cot_at(u01x, u01y, u01z, u02x, u02y, u02z);
    const geom_t cot1 = cot_at(-u01x, -u01y, -u01z, u12x, u12y, u12z);
    const geom_t cot2 = cot_at(-u02x, -u02y, -u02z, -u12x, -u12y, -u12z);

#pragma omp atomic
    area[i0] += a3;
#pragma omp atomic
    area[i1] += a3;
#pragma omp atomic
    area[i2] += a3;

    // Edge opposite vertex k contributes cot_k (x_i - x_j).
#pragma omp atomic
    kx[i0] += cot2 * (x0 - x1);
#pragma omp atomic
    ky[i0] += cot2 * (y0 - y1);
#pragma omp atomic
    kz[i0] += cot2 * (z0 - z1);
#pragma omp atomic
    kx[i1] += cot2 * (x1 - x0);
#pragma omp atomic
    ky[i1] += cot2 * (y1 - y0);
#pragma omp atomic
    kz[i1] += cot2 * (z1 - z0);

#pragma omp atomic
    kx[i1] += cot0 * (x1 - x2);
#pragma omp atomic
    ky[i1] += cot0 * (y1 - y2);
#pragma omp atomic
    kz[i1] += cot0 * (z1 - z2);
#pragma omp atomic
    kx[i2] += cot0 * (x2 - x1);
#pragma omp atomic
    ky[i2] += cot0 * (y2 - y1);
#pragma omp atomic
    kz[i2] += cot0 * (z2 - z1);

#pragma omp atomic
    kx[i2] += cot1 * (x2 - x0);
#pragma omp atomic
    ky[i2] += cot1 * (y2 - y0);
#pragma omp atomic
    kz[i2] += cot1 * (z2 - z0);
#pragma omp atomic
    kx[i0] += cot1 * (x0 - x2);
#pragma omp atomic
    ky[i0] += cot1 * (y0 - y2);
#pragma omp atomic
    kz[i0] += cot1 * (z0 - z2);
}

}  // namespace

template <typename idx_t, typename geom_t>
int mesh_node_curvature(const enum ElemType                                     element_type,
                        const ptrdiff_t                                         n_elements,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                        const int                                               sdim,
                        const ptrdiff_t                                         n_nodes,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                        const uint8_t *const SMESH_RESTRICT                     sharp_node,
                        geom_t *const SMESH_RESTRICT                            kappa) {
    if (!elements || !points || !kappa || n_nodes < 0 || n_elements < 0 || sdim < 2) {
        return SMESH_FAILURE;
    }
    const bool tri =
            element_type == TRI3 || element_type == TRISHELL3;
    const bool quad =
            element_type == QUAD4 || element_type == QUADSHELL4;
    if (!tri && !quad) {
        fprintf(stderr,
                "mesh_node_curvature: unsupported type %s\n",
                type_to_string(element_type));
        return SMESH_FAILURE;
    }

#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        kappa[i] = static_cast<geom_t>(0);
    }
    if (n_nodes == 0 || n_elements == 0) {
        return SMESH_SUCCESS;
    }

    geom_t *kx = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t *ky = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t *kz = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t *area = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    if (!kx || !ky || !kz || !area) {
        SMESH_FREE(kx);
        SMESH_FREE(ky);
        SMESH_FREE(kz);
        SMESH_FREE(area);
        return SMESH_FAILURE;
    }

    if (sdim >= 3) {
#pragma omp parallel for schedule(static)
        for (ptrdiff_t e = 0; e < n_elements; ++e) {
            const idx_t a = elements[0][e];
            const idx_t b = elements[1][e];
            const idx_t c = elements[2][e];
            accumulate_triangle(a, b, c, sdim, points, kx, ky, kz, area);
            if (quad) {
                const idx_t d = elements[3][e];
                accumulate_triangle(a, c, d, sdim, points, kx, ky, kz, area);
            }
        }

#pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            if (sharp_node && sharp_node[i]) {
                kappa[i] = static_cast<geom_t>(0);
                continue;
            }
            if (!(area[i] > static_cast<geom_t>(0))) {
                kappa[i] = static_cast<geom_t>(0);
                continue;
            }
            const geom_t inv = static_cast<geom_t>(0.5) / area[i];
            const geom_t hx  = kx[i] * inv;
            const geom_t hy  = ky[i] * inv;
            const geom_t hz  = kz[i] * inv;
            kappa[i]         = std::sqrt(hx * hx + hy * hy + hz * hz);
        }
    } else {
        // Planar: turning angle on the boundary polygon; interior κ = 0.
        const int nxe = tri ? 3 : 4;
        int *valence = (int *)SMESH_CALLOC((size_t)n_nodes, sizeof(int));
        geom_t *tx   = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
        geom_t *ty   = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
        if (!valence || !tx || !ty) {
            SMESH_FREE(kx);
            SMESH_FREE(ky);
            SMESH_FREE(kz);
            SMESH_FREE(area);
            SMESH_FREE(valence);
            SMESH_FREE(tx);
            SMESH_FREE(ty);
            return SMESH_FAILURE;
        }
        for (ptrdiff_t e = 0; e < n_elements; ++e) {
            for (int k = 0; k < nxe; ++k) {
                const idx_t a = elements[k][e];
                const idx_t b = elements[(k + 1) % nxe][e];
                const geom_t dx = points[0][b] - points[0][a];
                const geom_t dy = points[1][b] - points[1][a];
                tx[a] += dx;
                ty[a] += dy;
                tx[b] -= dx;
                ty[b] -= dy;
                valence[a] += 1;
                valence[b] += 1;
            }
        }
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            if (sharp_node && sharp_node[i]) {
                kappa[i] = static_cast<geom_t>(0);
                continue;
            }
            const geom_t n2 = tx[i] * tx[i] + ty[i] * ty[i];
            if (valence[i] <= 0 || !(n2 > static_cast<geom_t>(0))) {
                kappa[i] = static_cast<geom_t>(0);
                continue;
            }
            // |sum of unit-edge residuals| / mean incident length ~ turning / ē.
            kappa[i] = std::sqrt(n2);
        }
        SMESH_FREE(valence);
        SMESH_FREE(tx);
        SMESH_FREE(ty);
    }

    SMESH_FREE(kx);
    SMESH_FREE(ky);
    SMESH_FREE(kz);
    SMESH_FREE(area);
    return SMESH_SUCCESS;
}

template <typename geom_t>
int mesh_size_from_curvature(const ptrdiff_t              n_nodes,
                             const geom_t *const SMESH_RESTRICT kappa,
                             const geom_t                     cells_per_radius,
                             const geom_t                     h_min,
                             const geom_t                     h_max,
                             const geom_t                     bbox_diag,
                             geom_t *const SMESH_RESTRICT     h) {
    if (!kappa || !h || n_nodes < 0) {
        return SMESH_FAILURE;
    }
    const geom_t diag = bbox_diag > static_cast<geom_t>(0) ? bbox_diag : static_cast<geom_t>(1);
    const geom_t hlo =
            h_min > static_cast<geom_t>(0) ? h_min : diag * static_cast<geom_t>(1e-4);
    const geom_t hhi = h_max > static_cast<geom_t>(0) ? h_max : diag;
    const geom_t twopi = static_cast<geom_t>(6.28318530717958647692);
    const geom_t ncpr =
            cells_per_radius > static_cast<geom_t>(0) ? cells_per_radius : static_cast<geom_t>(8);
#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const geom_t k = kappa[i];
        geom_t       hi = hhi;
        if (k > static_cast<geom_t>(0)) {
            hi = twopi / (k * ncpr);
        }
        if (hi < hlo) {
            hi = hlo;
        }
        if (hi > hhi) {
            hi = hhi;
        }
        h[i] = hi;
    }
    return SMESH_SUCCESS;
}

template <typename geom_t>
int mesh_size_from_curvature_error(const ptrdiff_t              n_nodes,
                                   const geom_t *const SMESH_RESTRICT kappa,
                                   const geom_t                     cells_per_radius,
                                   const geom_t                     geom_error,
                                   const geom_t                     h_min,
                                   const geom_t                     h_max,
                                   const geom_t                     bbox_diag,
                                   geom_t *const SMESH_RESTRICT     h) {
    if (!kappa || !h || n_nodes < 0) {
        return SMESH_FAILURE;
    }
    const geom_t diag = bbox_diag > static_cast<geom_t>(0) ? bbox_diag : static_cast<geom_t>(1);
    const geom_t hlo =
            h_min > static_cast<geom_t>(0) ? h_min : diag * static_cast<geom_t>(1e-4);
    const geom_t huge = diag * static_cast<geom_t>(1e6);
    const geom_t ncpr =
            cells_per_radius > static_cast<geom_t>(0) ? cells_per_radius : static_cast<geom_t>(8);
    const geom_t pi = static_cast<geom_t>(3.14159265358979323846);
    geom_t       eps = geom_error;
    if (!(eps > static_cast<geom_t>(0))) {
        eps = (pi * pi) * diag / (static_cast<geom_t>(2) * ncpr * ncpr);
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const geom_t k = kappa[i];
        if (k * diag > static_cast<geom_t>(0.5)) {
            geom_t hi = std::sqrt(static_cast<geom_t>(8) * eps / k);
            if (hi < hlo) {
                hi = hlo;
            }
            if (h_max > static_cast<geom_t>(0) && hi > h_max) {
                hi = h_max;
            }
            h[i] = hi;
        } else {
            h[i] = huge;
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename geom_t>
int mesh_grade_size_2to1(const ptrdiff_t                     n_nodes,
                         const count_t *const SMESH_RESTRICT rowptr,
                         const idx_t *const SMESH_RESTRICT   colidx,
                         const int                           n_sweeps,
                         geom_t *const SMESH_RESTRICT        h) {
    if (!rowptr || !colidx || !h || n_nodes < 0) {
        return SMESH_FAILURE;
    }
    if (n_sweeps <= 0 || n_nodes == 0) {
        return SMESH_SUCCESS;
    }
    geom_t *tmp = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    if (!tmp) {
        return SMESH_FAILURE;
    }
    const geom_t two = static_cast<geom_t>(2);
    memcpy(tmp, h, (size_t)n_nodes * sizeof(geom_t));
    geom_t hmin = static_cast<geom_t>(0);
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        if (h[i] > static_cast<geom_t>(0) && (!(hmin > static_cast<geom_t>(0)) || h[i] < hmin)) {
            hmin = h[i];
        }
    }
    const geom_t hfree = (hmin > static_cast<geom_t>(0)) ? hmin * static_cast<geom_t>(1e4)
                                                         : static_cast<geom_t>(0);
    for (int s = 0; s < n_sweeps; ++s) {
        geom_t *src = (s % 2 == 0) ? h : tmp;
        geom_t *dst = (s % 2 == 0) ? tmp : h;
#pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            geom_t hi = src[i];
            const count_t b = rowptr[i];
            const count_t e = rowptr[i + 1];
            for (count_t k = b; k < e; ++k) {
                const geom_t hj = src[colidx[k]];
                if (!(hj > static_cast<geom_t>(0)) || (hfree > static_cast<geom_t>(0) && hj >= hfree)) {
                    continue;
                }
                const geom_t lo = hj / two;
                const geom_t hi2 = hj * two;
                if (hi < lo) {
                    hi = lo;
                }
                if (hi > hi2) {
                    hi = hi2;
                }
            }
            dst[i] = hi;
        }
    }
    if (n_sweeps % 2 == 1) {
        memcpy(h, tmp, (size_t)n_nodes * sizeof(geom_t));
    }
    SMESH_FREE(tmp);
    return SMESH_SUCCESS;
}

template <typename idx_t>
void mesh_mark_nodes_from_edges(const ptrdiff_t                   n_nodes,
                                const ptrdiff_t                   n_edges,
                                const idx_t *const SMESH_RESTRICT e0,
                                const idx_t *const SMESH_RESTRICT e1,
                                uint8_t *const SMESH_RESTRICT     mark) {
#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        mark[i] = 0;
    }
    for (ptrdiff_t e = 0; e < n_edges; ++e) {
        const idx_t a = e0[e];
        const idx_t b = e1[e];
        if (a >= 0 && (ptrdiff_t)a < n_nodes) {
            mark[a] = 1;
        }
        if (b >= 0 && (ptrdiff_t)b < n_nodes) {
            mark[b] = 1;
        }
    }
}

}  // namespace smesh

#endif
