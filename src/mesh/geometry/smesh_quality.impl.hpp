#ifndef SMESH_QUALITY_IMPL_HPP
#define SMESH_QUALITY_IMPL_HPP

#include "smesh_quality.hpp"

#include "smesh_alloc.hpp"
#include "smesh_elem_type.hpp"

#include <cmath>
#include <string.h>

namespace smesh {

namespace {

template <typename real_t>
SMESH_INLINE real_t tri_mean_ratio(real_t x0,
                                   real_t y0,
                                   real_t z0,
                                   real_t x1,
                                   real_t y1,
                                   real_t z1,
                                   real_t x2,
                                   real_t y2,
                                   real_t z2,
                                   const int sdim,
                                   const int signed_2d) {
    const real_t ux = x1 - x0, uy = y1 - y0, uz = z1 - z0;
    const real_t vx = x2 - x0, vy = y2 - y0, vz = z2 - z0;
    const real_t wx = x2 - x1, wy = y2 - y1, wz = z2 - z1;
    const real_t s2 = ux * ux + uy * uy + uz * uz + vx * vx + vy * vy + vz * vz + wx * wx + wy * wy +
                      wz * wz;
    if (!(s2 > static_cast<real_t>(0))) {
        return static_cast<real_t>(0);
    }
    real_t nx = uy * vz - uz * vy;
    real_t ny = uz * vx - ux * vz;
    real_t nz = ux * vy - uy * vx;
    if (sdim < 3) {
        nx = static_cast<real_t>(0);
        ny = static_cast<real_t>(0);
        nz = ux * vy - uy * vx;
        if (signed_2d && !(nz > static_cast<real_t>(0))) {
            return static_cast<real_t>(0);
        }
    }
    const real_t twice = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (!(twice > static_cast<real_t>(0))) {
        return static_cast<real_t>(0);
    }
    const real_t area       = static_cast<real_t>(0.5) * twice;
    const real_t four_sqrt3 = static_cast<real_t>(6.928203230275509);
    real_t       q          = four_sqrt3 * area / s2;
    if (q > static_cast<real_t>(1)) {
        q = static_cast<real_t>(1);
    }
    return q;
}

template <typename real_t>
SMESH_INLINE real_t tet_mean_ratio(real_t x0,
                                   real_t y0,
                                   real_t z0,
                                   real_t x1,
                                   real_t y1,
                                   real_t z1,
                                   real_t x2,
                                   real_t y2,
                                   real_t z2,
                                   real_t x3,
                                   real_t y3,
                                   real_t z3) {
    const real_t ux = x1 - x0, uy = y1 - y0, uz = z1 - z0;
    const real_t vx = x2 - x0, vy = y2 - y0, vz = z2 - z0;
    const real_t wx = x3 - x0, wy = y3 - y0, wz = z3 - z0;
    const real_t vol6 =
            ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx);
    if (!(vol6 > static_cast<real_t>(0))) {
        return static_cast<real_t>(0);
    }
    const real_t vol = vol6 / static_cast<real_t>(6);
    auto         l2  = [](real_t ax, real_t ay, real_t az, real_t bx, real_t by, real_t bz) {
        const real_t dx = bx - ax, dy = by - ay, dz = bz - az;
        return dx * dx + dy * dy + dz * dz;
    };
    const real_t s2 = l2(x0, y0, z0, x1, y1, z1) + l2(x1, y1, z1, x2, y2, z2) +
                      l2(x2, y2, z2, x0, y0, z0) + l2(x0, y0, z0, x3, y3, z3) +
                      l2(x1, y1, z1, x3, y3, z3) + l2(x2, y2, z2, x3, y3, z3);
    if (!(s2 > static_cast<real_t>(0))) {
        return static_cast<real_t>(0);
    }
    const real_t t = static_cast<real_t>(3) * vol;
    const real_t q =
            static_cast<real_t>(12) * std::pow(t, static_cast<real_t>(2) / static_cast<real_t>(3)) / s2;
    return q > static_cast<real_t>(1) ? static_cast<real_t>(1) : q;
}

}  // namespace

template <typename idx_t, typename real_t>
real_t mesh_elem_mean_ratio(const enum ElemType                                     element_type,
                            const int                                               sdim,
                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                            const real_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                            const ptrdiff_t                                         e) {
    auto P = [&](int d, idx_t i) -> real_t {
        if (d == 2 && sdim < 3) {
            return static_cast<real_t>(0);
        }
        return points[d][i];
    };
    const bool tri  = element_type == TRI3 || element_type == TRISHELL3;
    const bool quad = element_type == QUAD4 || element_type == QUADSHELL4;
    const bool tet  = element_type == TET4;
    if (tri) {
        const idx_t a = elements[0][e], b = elements[1][e], c = elements[2][e];
        return tri_mean_ratio(P(0, a), P(1, a), P(2, a), P(0, b), P(1, b), P(2, b), P(0, c), P(1, c),
                              P(2, c), sdim, sdim < 3);
    }
    if (tet) {
        if (sdim < 3) {
            return static_cast<real_t>(0);
        }
        const idx_t a = elements[0][e], b = elements[1][e], c = elements[2][e], d = elements[3][e];
        return tet_mean_ratio(P(0, a), P(1, a), P(2, a), P(0, b), P(1, b), P(2, b), P(0, c), P(1, c),
                              P(2, c), P(0, d), P(1, d), P(2, d));
    }
    if (quad) {
        const idx_t a = elements[0][e], b = elements[1][e], c = elements[2][e], d = elements[3][e];
        const real_t q02a = tri_mean_ratio(P(0, a), P(1, a), P(2, a), P(0, b), P(1, b), P(2, b),
                                           P(0, c), P(1, c), P(2, c), sdim, 0);
        const real_t q02b = tri_mean_ratio(P(0, a), P(1, a), P(2, a), P(0, c), P(1, c), P(2, c),
                                           P(0, d), P(1, d), P(2, d), sdim, 0);
        const real_t q13a = tri_mean_ratio(P(0, a), P(1, a), P(2, a), P(0, b), P(1, b), P(2, b),
                                           P(0, d), P(1, d), P(2, d), sdim, 0);
        const real_t q13b = tri_mean_ratio(P(0, b), P(1, b), P(2, b), P(0, c), P(1, c), P(2, c),
                                           P(0, d), P(1, d), P(2, d), sdim, 0);
        const real_t m02 = q02a < q02b ? q02a : q02b;
        const real_t m13 = q13a < q13b ? q13a : q13b;
        return m02 > m13 ? m02 : m13;
    }
    return static_cast<real_t>(0);
}

template <typename idx_t, typename real_t>
int mesh_element_quality(const enum ElemType                                     element_type,
                         const ptrdiff_t                                         n_elements,
                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                         const int                                               sdim,
                         const real_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                         real_t *const SMESH_RESTRICT                            q) {
    if (!elements || !points || !q || n_elements < 0 || sdim < 2) {
        return SMESH_FAILURE;
    }
    const bool ok = element_type == TRI3 || element_type == TRISHELL3 || element_type == QUAD4 ||
                    element_type == QUADSHELL4 || element_type == TET4;
    if (!ok) {
        return SMESH_FAILURE;
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
        q[e] = mesh_elem_mean_ratio<idx_t, real_t>(element_type, sdim, elements, points, e);
    }
    return SMESH_SUCCESS;
}

template <typename real_t>
real_t mesh_quality_min(const ptrdiff_t n_elements, const real_t *const SMESH_RESTRICT q) {
    if (!q || n_elements <= 0) {
        return static_cast<real_t>(0);
    }
    real_t m = q[0];
    for (ptrdiff_t e = 1; e < n_elements; ++e) {
        if (q[e] < m) {
            m = q[e];
        }
    }
    return m;
}

template <typename real_t>
void mesh_clamp_to_band(const int    sdim,
                        real_t      *x,
                        real_t      *y,
                        real_t      *z,
                        const real_t x0,
                        const real_t y0,
                        const real_t z0,
                        const real_t nx,
                        const real_t ny,
                        const real_t nz,
                        const real_t max_abs_dev,
                        const real_t max_normal_dev) {
    real_t dx = *x - x0;
    real_t dy = *y - y0;
    real_t dz = (sdim >= 3 && z) ? (*z - z0) : static_cast<real_t>(0);
    if (sdim >= 3 && max_normal_dev > static_cast<real_t>(0)) {
        const real_t dn = dx * nx + dy * ny + dz * nz;
        if (dn > max_normal_dev) {
            const real_t t = dn - max_normal_dev;
            dx -= t * nx;
            dy -= t * ny;
            dz -= t * nz;
        } else if (dn < -max_normal_dev) {
            const real_t t = dn + max_normal_dev;
            dx -= t * nx;
            dy -= t * ny;
            dz -= t * nz;
        }
    }
    if (max_abs_dev > static_cast<real_t>(0)) {
        const real_t r2 = dx * dx + dy * dy + dz * dz;
        const real_t m2 = max_abs_dev * max_abs_dev;
        if (r2 > m2 && r2 > static_cast<real_t>(0)) {
            const real_t s = max_abs_dev / std::sqrt(r2);
            dx *= s;
            dy *= s;
            dz *= s;
        }
    }
    *x = x0 + dx;
    *y = y0 + dy;
    if (sdim >= 3 && z) {
        *z = z0 + dz;
    }
}

namespace {

template <typename real_t>
SMESH_INLINE void closest_on_seg(const real_t px,
                                 const real_t py,
                                 const real_t pz,
                                 const real_t ax,
                                 const real_t ay,
                                 const real_t az,
                                 const real_t bx,
                                 const real_t by,
                                 const real_t bz,
                                 real_t      *qx,
                                 real_t      *qy,
                                 real_t      *qz,
                                 real_t      *d2) {
    const real_t ux = bx - ax, uy = by - ay, uz = bz - az;
    const real_t wx = px - ax, wy = py - ay, wz = pz - az;
    const real_t uu = ux * ux + uy * uy + uz * uz;
    real_t       t  = 0;
    if (uu > static_cast<real_t>(0)) {
        t = (wx * ux + wy * uy + wz * uz) / uu;
        if (t < static_cast<real_t>(0)) {
            t = static_cast<real_t>(0);
        } else if (t > static_cast<real_t>(1)) {
            t = static_cast<real_t>(1);
        }
    }
    *qx         = ax + t * ux;
    *qy         = ay + t * uy;
    *qz         = az + t * uz;
    const real_t dx = px - *qx, dy = py - *qy, dz = pz - *qz;
    *d2 = dx * dx + dy * dy + dz * dz;
}

template <typename real_t>
SMESH_INLINE void closest_on_tri(const real_t px,
                                 const real_t py,
                                 const real_t pz,
                                 const real_t ax,
                                 const real_t ay,
                                 const real_t az,
                                 const real_t bx,
                                 const real_t by,
                                 const real_t bz,
                                 const real_t cx,
                                 const real_t cy,
                                 const real_t cz,
                                 real_t      *qx,
                                 real_t      *qy,
                                 real_t      *qz,
                                 real_t      *d2) {
    const real_t abx = bx - ax, aby = by - ay, abz = bz - az;
    const real_t acx = cx - ax, acy = cy - ay, acz = cz - az;
    const real_t apx = px - ax, apy = py - ay, apz = pz - az;
    const real_t d1 = abx * apx + aby * apy + abz * apz;
    const real_t d2v = acx * apx + acy * apy + acz * apz;
    if (d1 <= static_cast<real_t>(0) && d2v <= static_cast<real_t>(0)) {
        *qx = ax;
        *qy = ay;
        *qz = az;
        *d2 = apx * apx + apy * apy + apz * apz;
        return;
    }
    const real_t bpx = px - bx, bpy = py - by, bpz = pz - bz;
    const real_t d3 = abx * bpx + aby * bpy + abz * bpz;
    const real_t d4 = acx * bpx + acy * bpy + acz * bpz;
    if (d3 >= static_cast<real_t>(0) && d4 <= d3) {
        *qx = bx;
        *qy = by;
        *qz = bz;
        *d2 = bpx * bpx + bpy * bpy + bpz * bpz;
        return;
    }
    const real_t vc = d1 * d4 - d3 * d2v;
    if (vc <= static_cast<real_t>(0) && d1 >= static_cast<real_t>(0) && d3 <= static_cast<real_t>(0)) {
        const real_t den = d1 - d3;
        const real_t v   = den > static_cast<real_t>(0) ? d1 / den : static_cast<real_t>(0);
        *qx              = ax + v * abx;
        *qy              = ay + v * aby;
        *qz              = az + v * abz;
        const real_t dx = px - *qx, dy = py - *qy, dz = pz - *qz;
        *d2 = dx * dx + dy * dy + dz * dz;
        return;
    }
    const real_t cpx = px - cx, cpy = py - cy, cpz = pz - cz;
    const real_t d5 = abx * cpx + aby * cpy + abz * cpz;
    const real_t d6 = acx * cpx + acy * cpy + acz * cpz;
    if (d6 >= static_cast<real_t>(0) && d5 <= d6) {
        *qx = cx;
        *qy = cy;
        *qz = cz;
        *d2 = cpx * cpx + cpy * cpy + cpz * cpz;
        return;
    }
    const real_t vb = d5 * d2v - d1 * d6;
    if (vb <= static_cast<real_t>(0) && d2v >= static_cast<real_t>(0) && d6 <= static_cast<real_t>(0)) {
        const real_t den = d2v - d6;
        const real_t w   = den > static_cast<real_t>(0) ? d2v / den : static_cast<real_t>(0);
        *qx              = ax + w * acx;
        *qy              = ay + w * acy;
        *qz              = az + w * acz;
        const real_t dx = px - *qx, dy = py - *qy, dz = pz - *qz;
        *d2 = dx * dx + dy * dy + dz * dz;
        return;
    }
    const real_t va = d3 * d6 - d5 * d4;
    if (va <= static_cast<real_t>(0) && (d4 - d3) >= static_cast<real_t>(0) &&
        (d5 - d6) >= static_cast<real_t>(0)) {
        const real_t den = (d4 - d3) + (d5 - d6);
        const real_t w =
                den > static_cast<real_t>(0) ? (d4 - d3) / den : static_cast<real_t>(0);
        *qx = bx + w * (cx - bx);
        *qy = by + w * (cy - by);
        *qz = bz + w * (cz - bz);
        const real_t dx = px - *qx, dy = py - *qy, dz = pz - *qz;
        *d2 = dx * dx + dy * dy + dz * dz;
        return;
    }
    const real_t den = va + vb + vc;
    const real_t v   = den > static_cast<real_t>(0) ? vb / den : static_cast<real_t>(0);
    const real_t w   = den > static_cast<real_t>(0) ? vc / den : static_cast<real_t>(0);
    *qx              = ax + abx * v + acx * w;
    *qy              = ay + aby * v + acy * w;
    *qz              = az + abz * v + acz * w;
    const real_t dx = px - *qx, dy = py - *qy, dz = pz - *qz;
    *d2 = dx * dx + dy * dy + dz * dz;
}

template <typename idx_t, typename geom_t>
SMESH_INLINE void accum_tri_n(const geom_t ax,
                              const geom_t ay,
                              const geom_t az,
                              const geom_t bx,
                              const geom_t by,
                              const geom_t bz,
                              const geom_t cx,
                              const geom_t cy,
                              const geom_t cz,
                              geom_t      *nx,
                              geom_t      *ny,
                              geom_t      *nz,
                              const idx_t  i0,
                              const idx_t  i1,
                              const idx_t  i2) {
    const geom_t ux = bx - ax, uy = by - ay, uz = bz - az;
    const geom_t vx = cx - ax, vy = cy - ay, vz = cz - az;
    const geom_t nnx = uy * vz - uz * vy;
    const geom_t nny = uz * vx - ux * vz;
    const geom_t nnz = ux * vy - uy * vx;
#pragma omp atomic
    nx[i0] += nnx;
#pragma omp atomic
    ny[i0] += nny;
#pragma omp atomic
    nz[i0] += nnz;
#pragma omp atomic
    nx[i1] += nnx;
#pragma omp atomic
    ny[i1] += nny;
#pragma omp atomic
    nz[i1] += nnz;
#pragma omp atomic
    nx[i2] += nnx;
#pragma omp atomic
    ny[i2] += nny;
#pragma omp atomic
    nz[i2] += nnz;
}

}  // namespace

template <typename idx_t, typename geom_t>
int mesh_vertex_normals_from_faces(const enum ElemType                                     element_type,
                                   const ptrdiff_t                                         n_elements,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                   const int                                               sdim,
                                   const ptrdiff_t                                         n_nodes,
                                   const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                   geom_t *const SMESH_RESTRICT                            nx,
                                   geom_t *const SMESH_RESTRICT                            ny,
                                   geom_t *const SMESH_RESTRICT                            nz) {
    if (!elements || !points || !nx || !ny || n_nodes < 0 || n_elements < 0 || sdim < 2) {
        return SMESH_FAILURE;
    }
    const bool tri  = element_type == TRI3 || element_type == TRISHELL3;
    const bool quad = element_type == QUAD4 || element_type == QUADSHELL4;
    if (!tri && !quad) {
        return SMESH_FAILURE;
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        nx[i] = static_cast<geom_t>(0);
        ny[i] = static_cast<geom_t>(0);
        if (nz) {
            nz[i] = static_cast<geom_t>(0);
        }
    }
    if (n_elements == 0 || n_nodes == 0) {
        return SMESH_SUCCESS;
    }
    geom_t *nz_w = nz;
    geom_t *nz_tmp = nullptr;
    if (!nz_w) {
        nz_tmp = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
        if (!nz_tmp) {
            return SMESH_FAILURE;
        }
        nz_w = nz_tmp;
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
        const idx_t a = elements[0][e];
        const idx_t b = elements[1][e];
        const idx_t c = elements[2][e];
        const geom_t az = sdim >= 3 ? points[2][a] : static_cast<geom_t>(0);
        const geom_t bz = sdim >= 3 ? points[2][b] : static_cast<geom_t>(0);
        const geom_t cz = sdim >= 3 ? points[2][c] : static_cast<geom_t>(0);
        accum_tri_n<idx_t, geom_t>(points[0][a],
                                   points[1][a],
                                   az,
                                   points[0][b],
                                   points[1][b],
                                   bz,
                                   points[0][c],
                                   points[1][c],
                                   cz,
                                   nx,
                                   ny,
                                   nz_w,
                                   a,
                                   b,
                                   c);
        if (quad) {
            const idx_t d  = elements[3][e];
            const geom_t dz = sdim >= 3 ? points[2][d] : static_cast<geom_t>(0);
            accum_tri_n<idx_t, geom_t>(points[0][a],
                                       points[1][a],
                                       az,
                                       points[0][c],
                                       points[1][c],
                                       cz,
                                       points[0][d],
                                       points[1][d],
                                       dz,
                                       nx,
                                       ny,
                                       nz_w,
                                       a,
                                       c,
                                       d);
        }
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const geom_t n2 = nx[i] * nx[i] + ny[i] * ny[i] + nz_w[i] * nz_w[i];
        if (n2 > static_cast<geom_t>(0)) {
            const geom_t inv = static_cast<geom_t>(1) / std::sqrt(n2);
            nx[i] *= inv;
            ny[i] *= inv;
            nz_w[i] *= inv;
        }
    }
    SMESH_FREE(nz_tmp);
    return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int mesh_project_to_surface(const enum ElemType                                     element_type,
                            const ptrdiff_t                                         n_elements,
                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                            const int                                               sdim,
                            const ptrdiff_t                                         n_surf_nodes,
                            const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT surf_points,
                            const ptrdiff_t                                         n_query,
                            const uint8_t *const SMESH_RESTRICT                     query_mask,
                            geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT      query_points,
                            const geom_t                                            max_move) {
    if (!elements || !surf_points || !query_points || n_elements < 0 || n_surf_nodes < 0 ||
        n_query < 0 || sdim < 2) {
        return SMESH_FAILURE;
    }
    const bool tri  = element_type == TRI3 || element_type == TRISHELL3;
    const bool quad = element_type == QUAD4 || element_type == QUADSHELL4;
    if (!tri && !quad) {
        return SMESH_FAILURE;
    }
    if (n_query == 0 || n_elements == 0) {
        return SMESH_SUCCESS;
    }
    if (n_surf_nodes == 0) {
        return SMESH_FAILURE;
    }

    const ptrdiff_t n_tri = quad ? (n_elements * 2) : n_elements;
    idx_t *t0 = (idx_t *)SMESH_ALLOC((size_t)n_tri * sizeof(idx_t));
    idx_t *t1 = (idx_t *)SMESH_ALLOC((size_t)n_tri * sizeof(idx_t));
    idx_t *t2 = (idx_t *)SMESH_ALLOC((size_t)n_tri * sizeof(idx_t));
    if (!t0 || !t1 || !t2) {
        SMESH_FREE(t0);
        SMESH_FREE(t1);
        SMESH_FREE(t2);
        return SMESH_FAILURE;
    }
    ptrdiff_t nt = 0;
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
        t0[nt] = elements[0][e];
        t1[nt] = elements[1][e];
        t2[nt] = elements[2][e];
        ++nt;
        if (quad) {
            t0[nt] = elements[0][e];
            t1[nt] = elements[2][e];
            t2[nt] = elements[3][e];
            ++nt;
        }
    }

    geom_t lo[3] = {surf_points[0][0], surf_points[1][0],
                    sdim >= 3 ? surf_points[2][0] : static_cast<geom_t>(0)};
    geom_t hi[3] = {lo[0], lo[1], lo[2]};
    for (ptrdiff_t i = 1; i < n_surf_nodes; ++i) {
        for (int d = 0; d < sdim; ++d) {
            const geom_t v = surf_points[d][i];
            if (v < lo[d]) {
                lo[d] = v;
            }
            if (v > hi[d]) {
                hi[d] = v;
            }
        }
    }
    geom_t ext[3];
    geom_t cell[3];
    int    res[3];
    int    nbin = 1;
    {
        ptrdiff_t r = 8;
        while (r * r * r < n_tri && r < 48) {
            r += 8;
        }
        for (int d = 0; d < 3; ++d) {
            ext[d] = hi[d] - lo[d];
            if (!(ext[d] > static_cast<geom_t>(0))) {
                ext[d] = static_cast<geom_t>(1);
            }
            res[d]  = (d < sdim) ? (int)r : 1;
            cell[d] = ext[d] / static_cast<geom_t>(res[d]);
            nbin *= res[d];
        }
    }
    ptrdiff_t *bptr = (ptrdiff_t *)SMESH_CALLOC((size_t)(nbin + 1), sizeof(ptrdiff_t));
    if (!bptr) {
        SMESH_FREE(t0);
        SMESH_FREE(t1);
        SMESH_FREE(t2);
        return SMESH_FAILURE;
    }

    auto bin_of = [&](geom_t x, geom_t y, geom_t z, int *ix, int *iy, int *iz) {
        int a = (int)((x - lo[0]) / cell[0]);
        int b = (int)((y - lo[1]) / cell[1]);
        int c = (int)((z - lo[2]) / cell[2]);
        if (a < 0) {
            a = 0;
        } else if (a >= res[0]) {
            a = res[0] - 1;
        }
        if (b < 0) {
            b = 0;
        } else if (b >= res[1]) {
            b = res[1] - 1;
        }
        if (c < 0) {
            c = 0;
        } else if (c >= res[2]) {
            c = res[2] - 1;
        }
        *ix = a;
        *iy = b;
        *iz = c;
        return (a * res[1] + b) * res[2] + c;
    };

    for (ptrdiff_t t = 0; t < nt; ++t) {
        const idx_t a = t0[t], b = t1[t], c = t2[t];
        geom_t      x0 = surf_points[0][a], x1 = surf_points[0][b], x2 = surf_points[0][c];
        geom_t      y0 = surf_points[1][a], y1 = surf_points[1][b], y2 = surf_points[1][c];
        geom_t      z0 = sdim >= 3 ? surf_points[2][a] : static_cast<geom_t>(0);
        geom_t      z1 = sdim >= 3 ? surf_points[2][b] : static_cast<geom_t>(0);
        geom_t      z2 = sdim >= 3 ? surf_points[2][c] : static_cast<geom_t>(0);
        geom_t      mn[3] = {x0, y0, z0}, mx[3] = {x0, y0, z0};
        const geom_t xs[3] = {x0, x1, x2}, ys[3] = {y0, y1, y2}, zs[3] = {z0, z1, z2};
        for (int k = 1; k < 3; ++k) {
            if (xs[k] < mn[0]) {
                mn[0] = xs[k];
            }
            if (xs[k] > mx[0]) {
                mx[0] = xs[k];
            }
            if (ys[k] < mn[1]) {
                mn[1] = ys[k];
            }
            if (ys[k] > mx[1]) {
                mx[1] = ys[k];
            }
            if (zs[k] < mn[2]) {
                mn[2] = zs[k];
            }
            if (zs[k] > mx[2]) {
                mx[2] = zs[k];
            }
        }
        int i0, j0, k0, i1, j1, k1;
        bin_of(mn[0], mn[1], mn[2], &i0, &j0, &k0);
        bin_of(mx[0], mx[1], mx[2], &i1, &j1, &k1);
        for (int ix = i0; ix <= i1; ++ix) {
            for (int iy = j0; iy <= j1; ++iy) {
                for (int iz = k0; iz <= k1; ++iz) {
                    const int id = (ix * res[1] + iy) * res[2] + iz;
                    bptr[id + 1] += 1;
                }
            }
        }
    }
    for (int i = 0; i < nbin; ++i) {
        bptr[i + 1] += bptr[i];
    }
    const ptrdiff_t nnz = bptr[nbin];
    ptrdiff_t *bind = (ptrdiff_t *)SMESH_ALLOC((size_t)nnz * sizeof(ptrdiff_t));
    ptrdiff_t *fill = (ptrdiff_t *)SMESH_ALLOC((size_t)(nbin + 1) * sizeof(ptrdiff_t));
    if (!bind || !fill) {
        SMESH_FREE(t0);
        SMESH_FREE(t1);
        SMESH_FREE(t2);
        SMESH_FREE(bptr);
        SMESH_FREE(bind);
        SMESH_FREE(fill);
        return SMESH_FAILURE;
    }
    memcpy(fill, bptr, (size_t)(nbin + 1) * sizeof(ptrdiff_t));
    for (ptrdiff_t t = 0; t < nt; ++t) {
        const idx_t a = t0[t], b = t1[t], c = t2[t];
        geom_t      x0 = surf_points[0][a], x1 = surf_points[0][b], x2 = surf_points[0][c];
        geom_t      y0 = surf_points[1][a], y1 = surf_points[1][b], y2 = surf_points[1][c];
        geom_t      z0 = sdim >= 3 ? surf_points[2][a] : static_cast<geom_t>(0);
        geom_t      z1 = sdim >= 3 ? surf_points[2][b] : static_cast<geom_t>(0);
        geom_t      z2 = sdim >= 3 ? surf_points[2][c] : static_cast<geom_t>(0);
        geom_t      mn[3] = {x0, y0, z0}, mxv[3] = {x0, y0, z0};
        const geom_t xs[3] = {x0, x1, x2}, ys[3] = {y0, y1, y2}, zs[3] = {z0, z1, z2};
        for (int k = 1; k < 3; ++k) {
            if (xs[k] < mn[0]) {
                mn[0] = xs[k];
            }
            if (xs[k] > mxv[0]) {
                mxv[0] = xs[k];
            }
            if (ys[k] < mn[1]) {
                mn[1] = ys[k];
            }
            if (ys[k] > mxv[1]) {
                mxv[1] = ys[k];
            }
            if (zs[k] < mn[2]) {
                mn[2] = zs[k];
            }
            if (zs[k] > mxv[2]) {
                mxv[2] = zs[k];
            }
        }
        int i0, j0, k0, i1, j1, k1;
        bin_of(mn[0], mn[1], mn[2], &i0, &j0, &k0);
        bin_of(mxv[0], mxv[1], mxv[2], &i1, &j1, &k1);
        for (int ix = i0; ix <= i1; ++ix) {
            for (int iy = j0; iy <= j1; ++iy) {
                for (int iz = k0; iz <= k1; ++iz) {
                    const int id = (ix * res[1] + iy) * res[2] + iz;
                    bind[fill[id]++] = t;
                }
            }
        }
    }
    SMESH_FREE(fill);

    const geom_t huge = static_cast<geom_t>(1e30);
#pragma omp parallel for schedule(static)
    for (ptrdiff_t q = 0; q < n_query; ++q) {
        if (query_mask && !query_mask[q]) {
            continue;
        }
        const geom_t px = query_points[0][q];
        const geom_t py = query_points[1][q];
        const geom_t pz = sdim >= 3 ? query_points[2][q] : static_cast<geom_t>(0);
        geom_t       best = huge, bx = px, by = py, bz = pz;
        int          ix, iy, iz;
        bin_of(px, py, pz, &ix, &iy, &iz);
        const int rmax = res[0] > res[1] ? (res[0] > res[2] ? res[0] : res[2])
                                         : (res[1] > res[2] ? res[1] : res[2]);
        int found = 0;
        for (int r = 0; r <= rmax; ++r) {
            const int x0 = ix - r < 0 ? 0 : ix - r;
            const int x1 = ix + r >= res[0] ? res[0] - 1 : ix + r;
            const int y0 = iy - r < 0 ? 0 : iy - r;
            const int y1 = iy + r >= res[1] ? res[1] - 1 : iy + r;
            const int z0 = iz - r < 0 ? 0 : iz - r;
            const int z1 = iz + r >= res[2] ? res[2] - 1 : iz + r;
            for (int ax = x0; ax <= x1; ++ax) {
                for (int ay = y0; ay <= y1; ++ay) {
                    for (int az = z0; az <= z1; ++az) {
                        if (r > 0) {
                            const int on = (ax == x0 || ax == x1 || ay == y0 || ay == y1 ||
                                            az == z0 || az == z1);
                            if (!on) {
                                continue;
                            }
                        }
                        const int id = (ax * res[1] + ay) * res[2] + az;
                        for (ptrdiff_t k = bptr[id]; k < bptr[id + 1]; ++k) {
                            const ptrdiff_t t = bind[k];
                            const idx_t     a = t0[t], b = t1[t], c = t2[t];
                            geom_t          qx, qy, qz, d2;
                            closest_on_tri(px,
                                           py,
                                           pz,
                                           surf_points[0][a],
                                           surf_points[1][a],
                                           sdim >= 3 ? surf_points[2][a] : static_cast<geom_t>(0),
                                           surf_points[0][b],
                                           surf_points[1][b],
                                           sdim >= 3 ? surf_points[2][b] : static_cast<geom_t>(0),
                                           surf_points[0][c],
                                           surf_points[1][c],
                                           sdim >= 3 ? surf_points[2][c] : static_cast<geom_t>(0),
                                           &qx,
                                           &qy,
                                           &qz,
                                           &d2);
                            if (d2 < best) {
                                best = d2;
                                bx   = qx;
                                by   = qy;
                                bz   = qz;
                                found = 1;
                            }
                        }
                    }
                }
            }
            if (found) {
                const geom_t reach = static_cast<geom_t>(r + 1);
                geom_t       cs    = cell[0];
                if (cell[1] < cs) {
                    cs = cell[1];
                }
                if (sdim >= 3 && cell[2] < cs) {
                    cs = cell[2];
                }
                if (best <= reach * reach * cs * cs) {
                    break;
                }
            }
        }
        if (!found || !(best < huge)) {
            for (ptrdiff_t t = 0; t < nt; ++t) {
                const idx_t a = t0[t], b = t1[t], c = t2[t];
                geom_t      qx, qy, qz, d2;
                closest_on_tri(px,
                               py,
                               pz,
                               surf_points[0][a],
                               surf_points[1][a],
                               sdim >= 3 ? surf_points[2][a] : static_cast<geom_t>(0),
                               surf_points[0][b],
                               surf_points[1][b],
                               sdim >= 3 ? surf_points[2][b] : static_cast<geom_t>(0),
                               surf_points[0][c],
                               surf_points[1][c],
                               sdim >= 3 ? surf_points[2][c] : static_cast<geom_t>(0),
                               &qx,
                               &qy,
                               &qz,
                               &d2);
                if (d2 < best) {
                    best = d2;
                    bx   = qx;
                    by   = qy;
                    bz   = qz;
                }
            }
        }
        if (max_move > static_cast<geom_t>(0) && best > max_move * max_move) {
            continue;
        }
        query_points[0][q] = bx;
        query_points[1][q] = by;
        if (sdim >= 3) {
            query_points[2][q] = bz;
        }
    }

    SMESH_FREE(t0);
    SMESH_FREE(t1);
    SMESH_FREE(t2);
    SMESH_FREE(bptr);
    SMESH_FREE(bind);
    return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int mesh_project_to_segments(const int                                               sdim,
                             const ptrdiff_t                                         n_seg,
                             const idx_t *const SMESH_RESTRICT                       e0,
                             const idx_t *const SMESH_RESTRICT                       e1,
                             const ptrdiff_t                                         n_seg_nodes,
                             const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT seg_points,
                             const ptrdiff_t                                         n_query,
                             const uint8_t *const SMESH_RESTRICT                     query_mask,
                             geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT      query_points,
                             const geom_t                                            max_move) {
    if (!seg_points || !query_points || n_seg < 0 || n_seg_nodes < 0 || n_query < 0 || sdim < 2) {
        return SMESH_FAILURE;
    }
    if (n_seg > 0 && (!e0 || !e1)) {
        return SMESH_FAILURE;
    }
    if (n_query == 0 || n_seg == 0) {
        return SMESH_SUCCESS;
    }
    const geom_t huge = static_cast<geom_t>(1e30);
#pragma omp parallel for schedule(static)
    for (ptrdiff_t q = 0; q < n_query; ++q) {
        if (query_mask && !query_mask[q]) {
            continue;
        }
        const geom_t px = query_points[0][q];
        const geom_t py = query_points[1][q];
        const geom_t pz = sdim >= 3 ? query_points[2][q] : static_cast<geom_t>(0);
        geom_t       best = huge, bx = px, by = py, bz = pz;
        for (ptrdiff_t s = 0; s < n_seg; ++s) {
            const idx_t a = e0[s], b = e1[s];
            if (a < 0 || b < 0 || (ptrdiff_t)a >= n_seg_nodes || (ptrdiff_t)b >= n_seg_nodes) {
                continue;
            }
            geom_t qx, qy, qz, d2;
            closest_on_seg(px,
                           py,
                           pz,
                           seg_points[0][a],
                           seg_points[1][a],
                           sdim >= 3 ? seg_points[2][a] : static_cast<geom_t>(0),
                           seg_points[0][b],
                           seg_points[1][b],
                           sdim >= 3 ? seg_points[2][b] : static_cast<geom_t>(0),
                           &qx,
                           &qy,
                           &qz,
                           &d2);
            if (d2 < best) {
                best = d2;
                bx   = qx;
                by   = qy;
                bz   = qz;
            }
        }
        if (best < huge) {
            if (max_move > static_cast<geom_t>(0) && best > max_move * max_move) {
                continue;
            }
            query_points[0][q] = bx;
            query_points[1][q] = by;
            if (sdim >= 3) {
                query_points[2][q] = bz;
            }
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
