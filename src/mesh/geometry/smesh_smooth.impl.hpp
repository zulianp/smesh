#ifndef SMESH_SMOOTH_IMPL_HPP
#define SMESH_SMOOTH_IMPL_HPP

#include "smesh_smooth.hpp"

#include "smesh_alloc.hpp"
#include "smesh_common.hpp"
#include "smesh_quality.hpp"

#include <cmath>
#include <string.h>

namespace smesh {

template <typename idx_t, typename count_t, typename geom_t>
int mesh_smooth_feature(const int                                               sdim,
                        const ptrdiff_t                                         n_nodes,
                        geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT      points,
                        const count_t *const SMESH_RESTRICT                     rowptr,
                        const idx_t *const SMESH_RESTRICT                       colidx,
                        const uint8_t *const SMESH_RESTRICT                     lock,
                        const ptrdiff_t                                         n_sharp,
                        const idx_t *const SMESH_RESTRICT                       e0,
                        const idx_t *const SMESH_RESTRICT                       e1,
                        const int                                               n_iters,
                        const geom_t                                            lambda,
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT x0,
                        const uint8_t *const SMESH_RESTRICT                     surface,
                        const geom_t                                            max_abs_dev,
                        const geom_t                                            max_normal_dev,
                        const geom_t *const SMESH_RESTRICT                     vnx,
                        const geom_t *const SMESH_RESTRICT                     vny,
                        const geom_t *const SMESH_RESTRICT                     vnz) {
    if (!points || !rowptr || !colidx || !lock || sdim < 2 || n_nodes < 0) {
        return SMESH_FAILURE;
    }
    if (n_iters <= 0 || n_nodes == 0) {
        return SMESH_SUCCESS;
    }
    const geom_t lam = lambda > static_cast<geom_t>(0) ? lambda : static_cast<geom_t>(0.5);

    geom_t *nx = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t *ny = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t *nz = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t *ax = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    geom_t *ay = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    geom_t *az = sdim >= 3 ? (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t)) : nullptr;
    geom_t *bx = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    geom_t *by = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    geom_t *bz = sdim >= 3 ? (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t)) : nullptr;
    idx_t  *c0 = (idx_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(idx_t));
    idx_t  *c1 = (idx_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(idx_t));
    int    *cn = (int *)SMESH_CALLOC((size_t)n_nodes, sizeof(int));
    if (!nx || !ny || !nz || !ax || !ay || !bx || !by || !c0 || !c1 || !cn ||
        (sdim >= 3 && (!az || !bz))) {
        SMESH_FREE(nx);
        SMESH_FREE(ny);
        SMESH_FREE(nz);
        SMESH_FREE(ax);
        SMESH_FREE(ay);
        SMESH_FREE(az);
        SMESH_FREE(bx);
        SMESH_FREE(by);
        SMESH_FREE(bz);
        SMESH_FREE(c0);
        SMESH_FREE(c1);
        SMESH_FREE(cn);
        return SMESH_FAILURE;
    }

#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        c0[i] = static_cast<idx_t>(-1);
        c1[i] = static_cast<idx_t>(-1);
        ax[i] = points[0][i];
        ay[i] = points[1][i];
        if (sdim >= 3) {
            az[i] = points[2][i];
        }
        nx[i] = vnx ? vnx[i] : static_cast<geom_t>(0);
        ny[i] = vny ? vny[i] : static_cast<geom_t>(0);
        nz[i] = vnz ? vnz[i] : static_cast<geom_t>(0);
    }
    for (ptrdiff_t e = 0; e < n_sharp; ++e) {
        const idx_t u = e0[e];
        const idx_t v = e1[e];
        if (u < 0 || v < 0 || (ptrdiff_t)u >= n_nodes || (ptrdiff_t)v >= n_nodes) {
            continue;
        }
        if (cn[u] == 0) {
            c0[u] = v;
        } else if (cn[u] == 1) {
            c1[u] = v;
        }
        cn[u] += 1;
        if (cn[v] == 0) {
            c0[v] = u;
        } else if (cn[v] == 1) {
            c1[v] = u;
        }
        cn[v] += 1;
    }

    geom_t *srcx = ax, *srcy = ay, *srcz = az;
    geom_t *dstx = bx, *dsty = by, *dstz = bz;
    for (int it = 0; it < n_iters; ++it) {
#pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            if (lock[i] == 2) {
                if (x0) {
                    dstx[i] = x0[0][i];
                    dsty[i] = x0[1][i];
                    if (sdim >= 3) {
                        dstz[i] = x0[2][i];
                    }
                } else {
                    dstx[i] = srcx[i];
                    dsty[i] = srcy[i];
                    if (sdim >= 3) {
                        dstz[i] = srcz[i];
                    }
                }
                continue;
            }
            if (surface && !surface[i]) {
                dstx[i] = srcx[i];
                dsty[i] = srcy[i];
                if (sdim >= 3) {
                    dstz[i] = srcz[i];
                }
                continue;
            }
            if (lock[i] == 1 && cn[i] == 2 && c0[i] >= 0 && c1[i] >= 0) {
                dstx[i] = srcx[i] + lam * (static_cast<geom_t>(0.5) * (srcx[c0[i]] + srcx[c1[i]]) - srcx[i]);
                dsty[i] = srcy[i] + lam * (static_cast<geom_t>(0.5) * (srcy[c0[i]] + srcy[c1[i]]) - srcy[i]);
                if (sdim >= 3) {
                    dstz[i] =
                            srcz[i] + lam * (static_cast<geom_t>(0.5) * (srcz[c0[i]] + srcz[c1[i]]) - srcz[i]);
                }
            } else {
                const count_t b = rowptr[i];
                const count_t e = rowptr[i + 1];
                const int     on_surf = !surface || surface[i];
                geom_t        mx = 0, my = 0, mz = 0;
                count_t       nacc = 0;
                for (count_t k = b; k < e; ++k) {
                    const idx_t j = colidx[k];
                    if (on_surf && surface && !surface[j]) {
                        continue;
                    }
                    mx += srcx[j];
                    my += srcy[j];
                    if (sdim >= 3) {
                        mz += srcz[j];
                    }
                    nacc += 1;
                }
                if (nacc == 0) {
                    dstx[i] = srcx[i];
                    dsty[i] = srcy[i];
                    if (sdim >= 3) {
                        dstz[i] = srcz[i];
                    }
                } else {
                    const geom_t w = static_cast<geom_t>(1) / static_cast<geom_t>(nacc);
                    mx *= w;
                    my *= w;
                    mz *= w;
                    geom_t dx = mx - srcx[i];
                    geom_t dy = my - srcy[i];
                    geom_t dz = sdim >= 3 ? (mz - srcz[i]) : static_cast<geom_t>(0);
                    if (sdim >= 3 && on_surf) {
                        const geom_t n2 = nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i];
                        if (n2 > static_cast<geom_t>(0)) {
                            const geom_t dn = dx * nx[i] + dy * ny[i] + dz * nz[i];
                            dx -= dn * nx[i];
                            dy -= dn * ny[i];
                            dz -= dn * nz[i];
                        }
                    }
                    dstx[i] = srcx[i] + lam * dx;
                    dsty[i] = srcy[i] + lam * dy;
                    if (sdim >= 3) {
                        dstz[i] = srcz[i] + lam * dz;
                    }
                }
            }
            const int clip = x0 && (max_abs_dev > static_cast<geom_t>(0) ||
                                    max_normal_dev > static_cast<geom_t>(0)) &&
                             (!surface || surface[i]);
            if (clip) {
                geom_t zz = sdim >= 3 ? dstz[i] : static_cast<geom_t>(0);
                mesh_clamp_to_band(sdim,
                                   &dstx[i],
                                   &dsty[i],
                                   sdim >= 3 ? &zz : nullptr,
                                   x0[0][i],
                                   x0[1][i],
                                   sdim >= 3 ? x0[2][i] : static_cast<geom_t>(0),
                                   nx[i],
                                   ny[i],
                                   nz[i],
                                   max_abs_dev,
                                   max_normal_dev);
                if (sdim >= 3) {
                    dstz[i] = zz;
                }
            }
        }
        geom_t *tx = srcx;
        srcx       = dstx;
        dstx       = tx;
        geom_t *ty = srcy;
        srcy       = dsty;
        dsty       = ty;
        geom_t *tz = srcz;
        srcz       = dstz;
        dstz       = tz;
    }

#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        points[0][i] = srcx[i];
        points[1][i] = srcy[i];
        if (sdim >= 3) {
            points[2][i] = srcz[i];
        }
    }

    SMESH_FREE(nx);
    SMESH_FREE(ny);
    SMESH_FREE(nz);
    SMESH_FREE(ax);
    SMESH_FREE(ay);
    SMESH_FREE(az);
    SMESH_FREE(bx);
    SMESH_FREE(by);
    SMESH_FREE(bz);
    SMESH_FREE(c0);
    SMESH_FREE(c1);
    SMESH_FREE(cn);
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
