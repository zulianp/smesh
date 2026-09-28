#ifndef SMESH_SMOOTH_HPP
#define SMESH_SMOOTH_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <stddef.h>
#include <cmath>

namespace smesh {

/// Local feature through `o,a,b`: a straight line when the three points are
/// colinear, otherwise their circumcircle. `tx` is a unit tangent at `o`.
/// `circle` is 1 when `c`/`r2`/`n` describe that circle (`r2` is radius squared).
template <typename geom_t>
inline int feature_curve_frame(geom_t ox, geom_t oy, geom_t oz,
                               geom_t ax, geom_t ay, geom_t az,
                               geom_t bx, geom_t by, geom_t bz,
                               geom_t *tx, geom_t *ty, geom_t *tz,
                               geom_t *cx, geom_t *cy, geom_t *cz,
                               geom_t *r2,
                               geom_t *nx, geom_t *ny, geom_t *nz,
                               int    *circle) {
    const geom_t aox = ax - ox, aoy = ay - oy, aoz = az - oz;
    const geom_t box = bx - ox, boy = by - oy, boz = bz - oz;
    const geom_t nnx = aoy * boz - aoz * boy;
    const geom_t nny = aoz * box - aox * boz;
    const geom_t nnz = aox * boy - aoy * box;
    const geom_t n2  = nnx * nnx + nny * nny + nnz * nnz;
    const geom_t oa2 = aox * aox + aoy * aoy + aoz * aoz;
    const geom_t ob2 = box * box + boy * boy + boz * boz;
    if (!(oa2 > static_cast<geom_t>(0)) || !(ob2 > static_cast<geom_t>(0))) {
        return 0;
    }
    geom_t ux, uy, uz;
    int    circ = 0;
    geom_t ccx = 0, ccy = 0, ccz = 0, rad2 = 0, unx = 0, uny = 0, unz = 0;
    if (n2 > static_cast<geom_t>(1e-12) * oa2 * ob2) {
        const geom_t abx = bx - ax, aby = by - ay, abz = bz - az;
        const geom_t la  = abx * abx + aby * aby + abz * abz;
        const geom_t wa  = la * (ob2 + oa2 - la);
        const geom_t wb  = ob2 * (oa2 + la - ob2);
        const geom_t wc  = oa2 * (la + ob2 - oa2);
        const geom_t w   = wa + wb + wc;
        const geom_t w2  = w * w;
        const geom_t scale = la + ob2 + oa2;
        if (!(w2 > static_cast<geom_t>(1e-24) * scale * scale)) {
            return 0;
        }
        const geom_t invw = static_cast<geom_t>(1) / w;
        ccx = (wa * ox + wb * ax + wc * bx) * invw;
        ccy = (wa * oy + wb * ay + wc * by) * invw;
        ccz = (wa * oz + wb * az + wc * bz) * invw;
        const geom_t rx = ox - ccx, ry = oy - ccy, rz = oz - ccz;
        rad2 = rx * rx + ry * ry + rz * rz;
        if (!(rad2 > static_cast<geom_t>(0))) {
            return 0;
        }
        const geom_t invn = static_cast<geom_t>(1) / std::sqrt(n2);
        unx = nnx * invn;
        uny = nny * invn;
        unz = nnz * invn;
        ux  = uny * rz - unz * ry;
        uy  = unz * rx - unx * rz;
        uz  = unx * ry - uny * rx;
        circ = 1;
    } else {
        ux = bx - ax;
        uy = by - ay;
        uz = bz - az;
    }
    const geom_t t2 = ux * ux + uy * uy + uz * uz;
    if (!(t2 > static_cast<geom_t>(0))) {
        return 0;
    }
    const geom_t invt = static_cast<geom_t>(1) / std::sqrt(t2);
    *tx = ux * invt;
    *ty = uy * invt;
    *tz = uz * invt;
    *cx = ccx;
    *cy = ccy;
    *cz = ccz;
    *r2 = rad2;
    *nx = unx;
    *ny = uny;
    *nz = unz;
    *circle = circ;
    return 1;
}

/// Move `o` by the part of `(dx,dy,dz)` tangent to the feature through `o,a,b`,
/// then snap back onto that feature. A step longer than about half the shorter
/// incident feature edge is clamped so the node cannot jump to the far arc.
template <typename geom_t>
inline int feature_curve_point(geom_t ox, geom_t oy, geom_t oz,
                               geom_t ax, geom_t ay, geom_t az,
                               geom_t bx, geom_t by, geom_t bz,
                               geom_t dx, geom_t dy, geom_t dz,
                               geom_t *x, geom_t *y, geom_t *z) {
    geom_t tx, ty, tz, cx, cy, cz, r2, nx, ny, nz;
    int    circle = 0;
    if (!feature_curve_frame(ox, oy, oz, ax, ay, az, bx, by, bz,
                             &tx, &ty, &tz, &cx, &cy, &cz, &r2, &nx, &ny, &nz, &circle)) {
        return 0;
    }
    geom_t s = dx * tx + dy * ty + dz * tz;
    const geom_t oa2 = (ax - ox) * (ax - ox) + (ay - oy) * (ay - oy) + (az - oz) * (az - oz);
    const geom_t ob2 = (bx - ox) * (bx - ox) + (by - oy) * (by - oy) + (bz - oz) * (bz - oz);
    const geom_t olen = oa2 < ob2 ? oa2 : ob2;
    const geom_t lim  = static_cast<geom_t>(0.49) * std::sqrt(olen);
    if (s > lim) {
        s = lim;
    } else if (s < -lim) {
        s = -lim;
    }
    if (!(s > static_cast<geom_t>(1e-8) * lim) && !(s < static_cast<geom_t>(-1e-8) * lim)) {
        return 0;
    }
    geom_t qx = ox + s * tx, qy = oy + s * ty, qz = oz + s * tz;
    if (circle) {
        geom_t vx = qx - cx, vy = qy - cy, vz = qz - cz;
        const geom_t dn = vx * nx + vy * ny + vz * nz;
        vx -= dn * nx;
        vy -= dn * ny;
        vz -= dn * nz;
        const geom_t v2 = vx * vx + vy * vy + vz * vz;
        if (!(v2 > static_cast<geom_t>(0))) {
            return 0;
        }
        const geom_t scale = std::sqrt(r2 / v2);
        qx = cx + vx * scale;
        qy = cy + vy * scale;
        qz = cz + vz * scale;
    }
    *x = qx;
    *y = qy;
    *z = qz;
    return 1;
}

/// Drop `q` onto the feature through `o,a,b` (line, or the near point of the circumcircle).
template <typename geom_t>
inline int feature_curve_project(geom_t qx, geom_t qy, geom_t qz,
                                 geom_t ox, geom_t oy, geom_t oz,
                                 geom_t ax, geom_t ay, geom_t az,
                                 geom_t bx, geom_t by, geom_t bz,
                                 geom_t *x, geom_t *y, geom_t *z) {
    geom_t tx, ty, tz, cx, cy, cz, r2, nx, ny, nz;
    int    circle = 0;
    if (!feature_curve_frame(ox, oy, oz, ax, ay, az, bx, by, bz,
                             &tx, &ty, &tz, &cx, &cy, &cz, &r2, &nx, &ny, &nz, &circle)) {
        return 0;
    }
    if (!circle) {
        const geom_t s = (qx - ox) * tx + (qy - oy) * ty + (qz - oz) * tz;
        *x = ox + s * tx;
        *y = oy + s * ty;
        *z = oz + s * tz;
        return 1;
    }
    geom_t vx = qx - cx, vy = qy - cy, vz = qz - cz;
    const geom_t dn = vx * nx + vy * ny + vz * nz;
    vx -= dn * nx;
    vy -= dn * ny;
    vz -= dn * nz;
    const geom_t v2 = vx * vx + vy * vy + vz * vz;
    if (!(v2 > static_cast<geom_t>(0))) {
        return 0;
    }
    const geom_t scale = std::sqrt(r2 / v2);
    *x = cx + vx * scale;
    *y = cy + vy * scale;
    *z = cz + vz * scale;
    return 1;
}

/// lock: 0 free (tangential Jacobi), 1 crease (1D along sharp edge), 2 pin.
/// `e0`/`e1` are sharp-edge node pairs (`n_sharp` of them). `n2n` is the full
/// undirected CRS. Jacobi uses two point buffers internally.
/// Optional `x0` + `max_abs_dev` / `max_normal_dev` clip surface (or all, if
/// `surface` is null) nodes after each write. `lock==2` copies `x0` when set.
/// `nx`/`ny`/`nz` are unit vertex normals; if null, no tangent projection.
/// Surface-marked nodes average only surface neighbors.
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
                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT x0 = nullptr,
                        const uint8_t *const SMESH_RESTRICT                     surface = nullptr,
                        const geom_t                                            max_abs_dev = 0,
                        const geom_t                                            max_normal_dev = 0,
                        const geom_t *const SMESH_RESTRICT                     nx = nullptr,
                        const geom_t *const SMESH_RESTRICT                     ny = nullptr,
                        const geom_t *const SMESH_RESTRICT                     nz = nullptr);

}  // namespace smesh

#endif
