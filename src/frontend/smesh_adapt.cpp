#include "smesh_mesh.hpp"

#include "smesh_adapt_refine.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_common.hpp"
#include "smesh_curvature.hpp"
#include "smesh_edgeset.hpp"
#include "smesh_extractions.hpp"
#include "smesh_geom_map.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_parametrization.hpp"
#include "smesh_sideset.hpp"
#include "smesh_improve.hpp"
#include "smesh_quality.hpp"
#include "smesh_smooth.hpp"
#include "smesh_tracer.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

namespace smesh {
namespace {

geom_t tet4_relvol(geom_t **p, idx_t a, idx_t b, idx_t c, idx_t d) {
    const geom_t abx = p[0][b] - p[0][a], aby = p[1][b] - p[1][a], abz = p[2][b] - p[2][a];
    const geom_t acx = p[0][c] - p[0][a], acy = p[1][c] - p[1][a], acz = p[2][c] - p[2][a];
    const geom_t adx = p[0][d] - p[0][a], ady = p[1][d] - p[1][a], adz = p[2][d] - p[2][a];
    const geom_t v =
            abx * (acy * adz - acz * ady) + aby * (acz * adx - acx * adz) + abz * (acx * ady - acy * adx);
    const geom_t s = std::fabs(abx) + std::fabs(aby) + std::fabs(abz) + std::fabs(acx) + std::fabs(acy) +
                     std::fabs(acz) + std::fabs(adx) + std::fabs(ady) + std::fabs(adz);
    return v / (s * s * s + std::numeric_limits<geom_t>::min());
}

void pull_mids_of_inverted_tets(Mesh &mesh, const idx_t *na, const idx_t *nb) {
    if (mesh.element_type(0) != TET4 || !na || !nb || mesh.spatial_dimension() < 3) {
        return;
    }
    idx_t **const   el = mesh.elements(0)->data();
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t ne = mesh.n_elements(0);
    for (int it = 0; it < 32; ++it) {
        int any = 0;
        for (ptrdiff_t e = 0; e < ne; ++e) {
            const idx_t a = el[0][e], b = el[1][e], c = el[2][e], d = el[3][e];
            if (tet4_relvol(p, a, b, c, d) >= static_cast<geom_t>(-1e-4)) {
                continue;
            }
            any = 1;
            const idx_t q[4] = {a, b, c, d};
            for (int k = 0; k < 4; ++k) {
                const idx_t i = q[k];
                if (na[i] == nb[i]) {
                    continue;
                }
                for (int dim = 0; dim < 3; ++dim) {
                    p[dim][i] = static_cast<geom_t>(0.5) * (p[dim][na[i]] + p[dim][nb[i]]);
                }
            }
        }
        if (!any) {
            return;
        }
    }
}

geom_t incident_qmin(const idx_t                                              i,
                     const count_t                                           *n2eptr,
                     const element_idx_t                                     *elindex,
                     idx_t **const                                            el,
                     geom_t **const                                           p) {
    geom_t mn  = 1;
    int    any = 0;
    for (count_t k = n2eptr[i]; k < n2eptr[i + 1]; ++k) {
        const ptrdiff_t e = (ptrdiff_t)elindex[k];
        const geom_t    q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e);
        if (!any || q < mn) {
            mn = q;
        }
        any = 1;
    }
    return any ? mn : static_cast<geom_t>(1);
}

void quality_laplace_tets(Mesh                    &mesh,
                          const uint8_t           *lock,
                          const count_t           *rowptr,
                          const idx_t             *colidx,
                          const count_t           *n2eptr,
                          const element_idx_t     *elindex,
                          const geom_t            *nx,
                          const geom_t            *ny,
                          const geom_t            *nz,
                          const int                n_iters,
                          const geom_t             lambda) {
    if (!lock || !rowptr || !colidx || !n2eptr || !elindex || n_iters <= 0) {
        return;
    }
    idx_t **const   el = mesh.elements(0)->data();
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t nn = mesh.n_nodes();
    const geom_t    lam =
            lambda > static_cast<geom_t>(0) ? lambda : static_cast<geom_t>(0.5);
    const geom_t step[3] = {lam, static_cast<geom_t>(0.5) * lam, static_cast<geom_t>(0.25) * lam};
    for (int it = 0; it < n_iters; ++it) {
        for (ptrdiff_t i = 0; i < nn; ++i) {
            if (lock[i] >= 2 || (nx && lock[i] >= 1)) {
                continue;
            }
            geom_t  mx = 0, my = 0, mz = 0;
            count_t nacc = 0;
            for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
                const idx_t j = colidx[k];
                if (j == (idx_t)i) {
                    continue;
                }
                mx += p[0][j];
                my += p[1][j];
                mz += p[2][j];
                nacc += 1;
            }
            if (nacc == 0) {
                continue;
            }
            const geom_t inv = static_cast<geom_t>(1) / (geom_t)nacc;
            geom_t       dx = mx * inv - p[0][i];
            geom_t       dy = my * inv - p[1][i];
            geom_t       dz = mz * inv - p[2][i];
            if (nx && ny && nz) {
                const geom_t nl2 = nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i];
                if (nl2 > static_cast<geom_t>(0)) {
                    const geom_t invn = static_cast<geom_t>(1) / std::sqrt(nl2);
                    const geom_t dtn  = (dx * nx[i] + dy * ny[i] + dz * nz[i]) * invn;
                    dx -= dtn * nx[i] * invn;
                    dy -= dtn * ny[i] * invn;
                    dz -= dtn * nz[i] * invn;
                }
            }
            const geom_t ox = p[0][i], oy = p[1][i], oz = p[2][i];
            const geom_t qold = incident_qmin((idx_t)i, n2eptr, elindex, el, p);
            geom_t       qbest = qold;
            geom_t       bx = ox, by = oy, bz = oz;
            for (int s = 0; s < 3; ++s) {
                p[0][i] = ox + step[s] * dx;
                p[1][i] = oy + step[s] * dy;
                p[2][i] = oz + step[s] * dz;
                const geom_t qn = incident_qmin((idx_t)i, n2eptr, elindex, el, p);
                if (qn > qbest) {
                    qbest = qn;
                    bx    = p[0][i];
                    by    = p[1][i];
                    bz    = p[2][i];
                }
            }
            p[0][i] = bx;
            p[1][i] = by;
            p[2][i] = bz;
        }
    }
}

void sliver_exorcise_tets(Mesh                    &mesh,
                          const uint8_t           *lock,
                          const count_t           *n2eptr,
                          const element_idx_t     *elindex,
                          const geom_t             qbar,
                          const geom_t            *tnx,
                          const geom_t            *tny,
                          const geom_t            *tnz) {
    if (!lock || !n2eptr || !elindex) {
        return;
    }
    idx_t **const   el = mesh.elements(0)->data();
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t ne = mesh.n_elements(0);
    const int       faces[4][4] = {{0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {1, 2, 3, 0}};
    for (int sweep = 0; sweep < 16; ++sweep) {
        int any = 0;
        for (ptrdiff_t e = 0; e < ne; ++e) {
            const geom_t qe = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e);
            if (!(qe < qbar)) {
                continue;
            }
            const idx_t v[4] = {el[0][e], el[1][e], el[2][e], el[3][e]};
            for (int f = 0; f < 4; ++f) {
                const idx_t iv = v[faces[f][3]];
                if (lock[iv] >= 2) {
                    continue;
                }
                const idx_t ia = v[faces[f][0]], ib = v[faces[f][1]], ic = v[faces[f][2]];
                geom_t nx = (p[1][ib] - p[1][ia]) * (p[2][ic] - p[2][ia]) -
                            (p[2][ib] - p[2][ia]) * (p[1][ic] - p[1][ia]);
                geom_t ny = (p[2][ib] - p[2][ia]) * (p[0][ic] - p[0][ia]) -
                            (p[0][ib] - p[0][ia]) * (p[2][ic] - p[2][ia]);
                geom_t nz = (p[0][ib] - p[0][ia]) * (p[1][ic] - p[1][ia]) -
                            (p[1][ib] - p[1][ia]) * (p[0][ic] - p[0][ia]);
                const geom_t nl = std::sqrt(nx * nx + ny * ny + nz * nz);
                if (!(nl > static_cast<geom_t>(0))) {
                    continue;
                }
                nx /= nl;
                ny /= nl;
                nz /= nl;
                if (tnx && tny && tnz && lock[iv] < 2) {
                    const geom_t nl2 = tnx[iv] * tnx[iv] + tny[iv] * tny[iv] + tnz[iv] * tnz[iv];
                    if (nl2 > static_cast<geom_t>(0)) {
                        const geom_t invn = static_cast<geom_t>(1) / std::sqrt(nl2);
                        const geom_t dtn  = (nx * tnx[iv] + ny * tny[iv] + nz * tnz[iv]) * invn;
                        nx -= dtn * tnx[iv] * invn;
                        ny -= dtn * tny[iv] * invn;
                        nz -= dtn * tnz[iv] * invn;
                        const geom_t tl = std::sqrt(nx * nx + ny * ny + nz * nz);
                        if (!(tl > static_cast<geom_t>(0))) {
                            continue;
                        }
                        nx /= tl;
                        ny /= tl;
                        nz /= tl;
                    }
                }
                const geom_t h    = static_cast<geom_t>(0.15) * nl;
                const geom_t ox   = p[0][iv], oy = p[1][iv], oz = p[2][iv];
                const geom_t qold = incident_qmin(iv, n2eptr, elindex, el, p);
                geom_t       qbest = qold;
                geom_t       bx = ox, by = oy, bz = oz;
                const geom_t frac[6] = {static_cast<geom_t>(0.05),
                                        static_cast<geom_t>(0.1),
                                        static_cast<geom_t>(0.25),
                                        static_cast<geom_t>(0.5),
                                        static_cast<geom_t>(1),
                                        static_cast<geom_t>(-0.25)};
                for (int s = 0; s < 6; ++s) {
                    p[0][iv] = ox + frac[s] * h * nx;
                    p[1][iv] = oy + frac[s] * h * ny;
                    p[2][iv] = oz + frac[s] * h * nz;
                    const geom_t qn = incident_qmin(iv, n2eptr, elindex, el, p);
                    if (qn > qbest) {
                        qbest = qn;
                        bx    = p[0][iv];
                        by    = p[1][iv];
                        bz    = p[2][iv];
                        any   = 1;
                    }
                }
                p[0][iv] = bx;
                p[1][iv] = by;
                p[2][iv] = bz;
            }
            {
                geom_t elen = 0;
                int    nee  = 0;
                const int ed[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                for (int k = 0; k < 6; ++k) {
                    const idx_t a = v[ed[k][0]], b = v[ed[k][1]];
                    const geom_t dx = p[0][b] - p[0][a], dy = p[1][b] - p[1][a], dz = p[2][b] - p[2][a];
                    elen += std::sqrt(dx * dx + dy * dy + dz * dz);
                    nee += 1;
                }
                const geom_t h = static_cast<geom_t>(0.04) * (elen / (geom_t)nee);
                const geom_t dir[6][3] = {{h, 0, 0},
                                          {-h, 0, 0},
                                          {0, h, 0},
                                          {0, -h, 0},
                                          {0, 0, h},
                                          {0, 0, -h}};
                for (int k = 0; k < 4; ++k) {
                    const idx_t iv = v[k];
                    if (lock[iv] >= 2) {
                        continue;
                    }
                    const geom_t ox = p[0][iv], oy = p[1][iv], oz = p[2][iv];
                    const geom_t qold = incident_qmin(iv, n2eptr, elindex, el, p);
                    geom_t       qbest = qold;
                    geom_t       bx = ox, by = oy, bz = oz;
                    for (int s = 0; s < 6; ++s) {
                        geom_t dx = dir[s][0], dy = dir[s][1], dz = dir[s][2];
                        if (tnx && tny && tnz) {
                            const geom_t nl2 =
                                    tnx[iv] * tnx[iv] + tny[iv] * tny[iv] + tnz[iv] * tnz[iv];
                            if (nl2 > static_cast<geom_t>(0)) {
                                const geom_t invn = static_cast<geom_t>(1) / std::sqrt(nl2);
                                const geom_t dtn  = (dx * tnx[iv] + dy * tny[iv] + dz * tnz[iv]) * invn;
                                dx -= dtn * tnx[iv] * invn;
                                dy -= dtn * tny[iv] * invn;
                                dz -= dtn * tnz[iv] * invn;
                            }
                        }
                        p[0][iv] = ox + dx;
                        p[1][iv] = oy + dy;
                        p[2][iv] = oz + dz;
                        const geom_t qn = incident_qmin(iv, n2eptr, elindex, el, p);
                        if (qn > qbest) {
                            qbest = qn;
                            bx    = p[0][iv];
                            by    = p[1][iv];
                            bz    = p[2][iv];
                            any   = 1;
                        }
                    }
                    p[0][iv] = bx;
                    p[1][iv] = by;
                    p[2][iv] = bz;
                }
            }
        }
        if (!any) {
            break;
        }
    }
}

void quality_gated_param_apply(Mesh                &mesh,
                               const count_t       *n2eptr,
                               const element_idx_t *elindex) {
    if (mesh.parametrizations().empty()) {
        return;
    }
    const ptrdiff_t nn = mesh.n_nodes();
    geom_t **const  p  = mesh.points()->data();
    geom_t         *pre[3] = {nullptr, nullptr, nullptr};
    geom_t         *snp[3] = {nullptr, nullptr, nullptr};
    int             ok     = 1;
    for (int d = 0; d < 3; ++d) {
        pre[d] = (geom_t *)SMESH_ALLOC((size_t)nn * sizeof(geom_t));
        snp[d] = (geom_t *)SMESH_ALLOC((size_t)nn * sizeof(geom_t));
        if (!pre[d] || !snp[d]) {
            ok = 0;
            break;
        }
        std::memcpy(pre[d], p[d], (size_t)nn * sizeof(geom_t));
    }
    if (!ok) {
        for (int d = 0; d < 3; ++d) {
            SMESH_FREE(pre[d]);
            SMESH_FREE(snp[d]);
        }
        for (const auto &kv : mesh.parametrizations()) {
            if (kv.second) {
                kv.second->apply(mesh);
            }
        }
        return;
    }
    for (const auto &kv : mesh.parametrizations()) {
        if (kv.second) {
            kv.second->apply(mesh);
        }
    }
    for (int d = 0; d < 3; ++d) {
        std::memcpy(snp[d], p[d], (size_t)nn * sizeof(geom_t));
        std::memcpy(p[d], pre[d], (size_t)nn * sizeof(geom_t));
    }
    if (mesh.element_type(0) == TET4 && n2eptr && elindex) {
        idx_t **const el = mesh.elements(0)->data();
        const geom_t  frac[4] = {static_cast<geom_t>(0.25),
                                static_cast<geom_t>(0.5),
                                static_cast<geom_t>(0.75),
                                static_cast<geom_t>(1)};
        for (const auto &kv : mesh.parametrizations()) {
            if (!kv.second || !kv.second->nodeset()) {
                continue;
            }
            auto         ns  = kv.second->nodeset();
            const idx_t *ids = ns->nodes()->data();
            for (ptrdiff_t k = 0; k < ns->size(); ++k) {
                const idx_t i = ids[k];
                if (i < 0 || (ptrdiff_t)i >= nn) {
                    continue;
                }
                const geom_t ox = p[0][i], oy = p[1][i], oz = p[2][i];
                const geom_t sx = snp[0][i], sy = snp[1][i], sz = snp[2][i];
                const geom_t q0 = incident_qmin(i, n2eptr, elindex, el, p);
                geom_t       bf = 0;
                for (int s = 0; s < 4; ++s) {
                    p[0][i] = ox + frac[s] * (sx - ox);
                    p[1][i] = oy + frac[s] * (sy - oy);
                    p[2][i] = oz + frac[s] * (sz - oz);
                    if (incident_qmin(i, n2eptr, elindex, el, p) >= q0) {
                        bf = frac[s];
                    }
                }
                p[0][i] = ox + bf * (sx - ox);
                p[1][i] = oy + bf * (sy - oy);
                p[2][i] = oz + bf * (sz - oz);
            }
        }
    } else {
        for (int d = 0; d < 3; ++d) {
            std::memcpy(p[d], snp[d], (size_t)nn * sizeof(geom_t));
        }
    }
    for (int d = 0; d < 3; ++d) {
        SMESH_FREE(pre[d]);
        SMESH_FREE(snp[d]);
    }
}

void restore_positive_tets(Mesh &mesh, geom_t **x0) {
    if (mesh.element_type(0) != TET4 || !x0 || mesh.spatial_dimension() < 3) {
        return;
    }
    idx_t **const   el = mesh.elements(0)->data();
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t ne = mesh.n_elements(0);
    for (int it = 0; it < 12; ++it) {
        int any = 0;
        for (ptrdiff_t e = 0; e < ne; ++e) {
            if (mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e) > static_cast<geom_t>(0)) {
                continue;
            }
            any = 1;
            const idx_t q[4] = {el[0][e], el[1][e], el[2][e], el[3][e]};
            int         kmax = 0;
            geom_t      dmax = -1;
            for (int k = 0; k < 4; ++k) {
                const idx_t i = q[k];
                geom_t      s = 0;
                for (int dim = 0; dim < 3; ++dim) {
                    const geom_t t = p[dim][i] - x0[dim][i];
                    s += t * t;
                }
                if (s > dmax) {
                    dmax = s;
                    kmax = k;
                }
            }
            const idx_t i = q[kmax];
            geom_t      cx[3], ox[3];
            for (int dim = 0; dim < 3; ++dim) {
                cx[dim] = p[dim][i];
                ox[dim] = x0[dim][i];
            }
            geom_t lo = 0, hi = 1;
            for (int bs = 0; bs < 10; ++bs) {
                const geom_t mid = static_cast<geom_t>(0.5) * (lo + hi);
                for (int dim = 0; dim < 3; ++dim) {
                    p[dim][i] = ox[dim] + mid * (cx[dim] - ox[dim]);
                }
                if (mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e) > static_cast<geom_t>(0)) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            for (int dim = 0; dim < 3; ++dim) {
                p[dim][i] = ox[dim] + lo * (cx[dim] - ox[dim]);
            }
            if (!(mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e) > static_cast<geom_t>(0))) {
                for (int k = 0; k < 4; ++k) {
                    const idx_t v = q[k];
                    for (int dim = 0; dim < 3; ++dim) {
                        p[dim][v] = x0[dim][v];
                    }
                }
            }
        }
        if (!any) {
            return;
        }
    }
}

geom_t mesh_tet_qmin(const Mesh &mesh) {
    if (mesh.element_type(0) != TET4) {
        return static_cast<geom_t>(1);
    }
    idx_t **const   el = mesh.elements(0)->data();
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t ne = mesh.n_elements(0);
    geom_t          mn = 1;
    for (ptrdiff_t e = 0; e < ne; ++e) {
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e);
        if (q < mn) {
            mn = q;
        }
    }
    return mn;
}

ptrdiff_t mesh_tet_nbelow(const Mesh &mesh, const geom_t qbar) {
    if (mesh.element_type(0) != TET4) {
        return 0;
    }
    idx_t **const   el = mesh.elements(0)->data();
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t ne = mesh.n_elements(0);
    ptrdiff_t       n  = 0;
    for (ptrdiff_t e = 0; e < ne; ++e) {
        n += mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e) < qbar;
    }
    return n;
}

geom_t bbox_diag(const Mesh &mesh) {
    const int       sdim = mesh.spatial_dimension();
    const ptrdiff_t n    = mesh.n_nodes();
    geom_t        **p    = mesh.points()->data();
    geom_t          s    = 0;
    for (int d = 0; d < sdim; ++d) {
        geom_t lo = 0, hi = 0;
        minmax(n, p[d], &lo, &hi);
        const geom_t e = hi - lo;
        s += e * e;
    }
    return std::sqrt(s);
}

void tighten_h_from_param_error(Mesh &surf, const Mesh &src, geom_t *h, const geom_t eps) {
    if (!h || !(eps > static_cast<geom_t>(0)) || src.parametrizations().empty()) {
        return;
    }
    auto n2n = surf.node_to_node_graph();
    if (!n2n) {
        return;
    }
    const count_t *const rp = n2n->rowptr()->data();
    const idx_t   *const ci = n2n->colidx()->data();
    const ptrdiff_t      nn = surf.n_nodes();
    const int            sdim = surf.spatial_dimension();
    geom_t **const       p    = surf.points()->data();
    if (!rp || !ci || !p || nn <= 0 || sdim < 2) {
        return;
    }

    for (const auto &kv : src.parametrizations()) {
        if (!kv.second) {
            continue;
        }
        auto ns = kv.second->nodeset();
        if (!ns || ns->size() <= 0) {
            continue;
        }
        std::vector<uint8_t> on((size_t)nn, 0);
        const idx_t *ids = ns->nodes()->data();
        for (ptrdiff_t i = 0; i < ns->size(); ++i) {
            const idx_t id = ids[i];
            if (id >= 0 && (ptrdiff_t)id < nn) {
                on[(size_t)id] = 1;
            }
        }

        ptrdiff_t n_mids = 0;
        for (ptrdiff_t i = 0; i < nn; ++i) {
            for (count_t k = rp[i]; k < rp[i + 1]; ++k) {
                const idx_t j = ci[k];
                if ((ptrdiff_t)j > i && on[(size_t)i] && on[(size_t)j]) {
                    ++n_mids;
                }
            }
        }
        if (n_mids == 0) {
            continue;
        }

        auto pbuf = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_mids);
        auto mid_ids = create_host_buffer<idx_t>((size_t)n_mids);
        std::vector<idx_t> ea((size_t)n_mids), eb((size_t)n_mids);
        std::vector<geom_t> chord((size_t)n_mids * (size_t)sdim);
        ptrdiff_t t = 0;
        for (ptrdiff_t i = 0; i < nn; ++i) {
            for (count_t k = rp[i]; k < rp[i + 1]; ++k) {
                const idx_t j = ci[k];
                if ((ptrdiff_t)j <= i || !on[(size_t)i] || !on[(size_t)j]) {
                    continue;
                }
                for (int d = 0; d < sdim; ++d) {
                    const geom_t m = static_cast<geom_t>(0.5) * (p[d][i] + p[d][j]);
                    pbuf->data()[d][t]              = m;
                    chord[(size_t)t * (size_t)sdim + (size_t)d] = m;
                }
                mid_ids->data()[t] = (idx_t)t;
                ea[(size_t)t]      = (idx_t)i;
                eb[(size_t)t]      = j;
                ++t;
            }
        }

        auto ebuf = create_host_buffer<idx_t>(3, 1);
        ebuf->data()[0][0] = 0;
        ebuf->data()[1][0] = n_mids > 1 ? 1 : 0;
        ebuf->data()[2][0] = n_mids > 2 ? 2 : 0;
        Mesh dummy(surf.comm(), TRI3, ebuf, pbuf);
        auto nsm = Nodeset::create(surf.comm(), mid_ids);
        auto np  = kv.second->with_nodeset(nsm);
        if (!np || np->apply(dummy) != SMESH_SUCCESS) {
            continue;
        }
        geom_t **const q = dummy.points()->data();
        for (ptrdiff_t m = 0; m < n_mids; ++m) {
            const idx_t a = ea[(size_t)m];
            const idx_t b = eb[(size_t)m];
            geom_t      err2 = 0;
            geom_t      l2   = 0;
            for (int d = 0; d < sdim; ++d) {
                const geom_t de = q[d][m] - chord[(size_t)m * (size_t)sdim + (size_t)d];
                err2 += de * de;
                const geom_t dl = p[d][b] - p[d][a];
                l2 += dl * dl;
            }
            if (err2 <= eps * eps || !(l2 > static_cast<geom_t>(0))) {
                continue;
            }
            const geom_t hi = static_cast<geom_t>(0.49) * std::sqrt(l2);
            if (hi < h[a]) {
                h[a] = hi;
            }
            if (hi < h[b]) {
                h[b] = hi;
            }
        }
    }
}

int fill_lock_from_sharp(Mesh     &surf,
                         uint8_t  *lock,
                         ptrdiff_t n_nodes,
                         geom_t    cos_th,
                         idx_t   **e0_out,
                         idx_t   **e1_out,
                         ptrdiff_t *n_sharp_out) {
    std::memset(lock, 0, (size_t)n_nodes);
    auto se = extract_sharp_edges(surf, cos_th);
    if (!se) {
        *e0_out      = nullptr;
        *e1_out      = nullptr;
        *n_sharp_out = 0;
        return SMESH_SUCCESS;
    }
    auto corners = extract_sharp_corners(surf, se, false);
    LocalEdgeTable let;
    auto           blk = surf.block(0);
    if (!blk || let.fill(blk->element_type()) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    const ptrdiff_t nse = se->size();
    idx_t          *e0  = nse > 0 ? (idx_t *)SMESH_ALLOC((size_t)nse * sizeof(idx_t)) : nullptr;
    idx_t          *e1  = nse > 0 ? (idx_t *)SMESH_ALLOC((size_t)nse * sizeof(idx_t)) : nullptr;
    idx_t **const   el  = blk->elements()->data();
    const element_idx_t *par = nse > 0 ? se->parent()->data() : nullptr;
    const i16           *lei = nse > 0 ? se->lei()->data() : nullptr;
    for (ptrdiff_t i = 0; i < nse; ++i) {
        const ptrdiff_t e = (ptrdiff_t)par[i];
        const int       s = (int)lei[i];
        e0[i]             = el[let(s, 0)][e];
        e1[i]             = el[let(s, 1)][e];
        if (e0[i] >= 0 && (ptrdiff_t)e0[i] < n_nodes && lock[e0[i]] < 1) {
            lock[e0[i]] = 1;
        }
        if (e1[i] >= 0 && (ptrdiff_t)e1[i] < n_nodes && lock[e1[i]] < 1) {
            lock[e1[i]] = 1;
        }
    }
    if (corners && corners->size() > 0) {
        const idx_t *ids = corners->nodes()->data();
        for (ptrdiff_t i = 0; i < corners->size(); ++i) {
            const idx_t id = ids[i];
            if (id >= 0 && (ptrdiff_t)id < n_nodes) {
                lock[id] = 2;
            }
        }
    }
    *e0_out      = e0;
    *e1_out      = e1;
    *n_sharp_out = nse;
    return SMESH_SUCCESS;
}

std::shared_ptr<Mesh> wrap_adapt_mesh(const std::shared_ptr<Mesh> &coarse,
                                      enum ElemType                et,
                                      const int                    nxe,
                                      const int                    sdim,
                                      ptrdiff_t                    n_elem,
                                      idx_t                      **elems,
                                      ptrdiff_t                    n_nodes,
                                      geom_t                     **pts) {
    auto ebuf = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_elem);
    auto pbuf = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_nodes);
    for (int d = 0; d < nxe; ++d) {
        std::memcpy(ebuf->data()[d], elems[d], (size_t)n_elem * sizeof(idx_t));
    }
    for (int d = 0; d < sdim; ++d) {
        std::memcpy(pbuf->data()[d], pts[d], (size_t)n_nodes * sizeof(geom_t));
    }
    std::string name = "default";
    if (coarse->n_blocks() > 0) {
        name = coarse->block(0)->name();
    }
    auto block = std::make_shared<Mesh::Block>(name, et, ebuf);
    const enum GeomMap gm =
            detect_geom_map(et, n_elem, ebuf->data(), sdim, pbuf->data());
    block->set_geom_map(gm);
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(block);
    return std::make_shared<Mesh>(coarse->comm(), blocks, pbuf);
}

void remap_nodeset_mids(const ptrdiff_t                   n_coarse_nodes,
                        const ptrdiff_t                   n_nodes,
                        const idx_t *const SMESH_RESTRICT node_a,
                        const idx_t *const SMESH_RESTRICT node_b,
                        uint8_t *const                    in) {
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const idx_t a = node_a[i];
        const idx_t b = node_b[i];
        if (a >= 0 && b >= 0 && in[a] && in[b]) {
            in[i] = 1;
        }
    }
    (void)n_coarse_nodes;
}

std::shared_ptr<Nodeset> expand_nodeset(const std::shared_ptr<Mesh>    &coarse,
                                        const std::shared_ptr<Nodeset> &ns,
                                        const std::shared_ptr<Mesh>    &fine,
                                        const idx_t                    *node_a,
                                        const idx_t                    *node_b) {
    if (!ns) {
        return nullptr;
    }
    const ptrdiff_t n_c = coarse->n_nodes();
    const ptrdiff_t n_f = fine->n_nodes();
    std::vector<uint8_t> in((size_t)n_f, 0);
    if (ns->size() > 0) {
        const idx_t *ids = ns->nodes()->data();
        for (ptrdiff_t i = 0; i < ns->size(); ++i) {
            const idx_t id = ids[i];
            if (id >= 0 && (ptrdiff_t)id < n_c) {
                in[(size_t)id] = 1;
            }
        }
    }
    remap_nodeset_mids(n_c, n_f, node_a, node_b, in.data());
    ptrdiff_t count = 0;
    for (ptrdiff_t i = 0; i < n_f; ++i) {
        count += in[(size_t)i];
    }
    auto buf = create_host_buffer<idx_t>((size_t)count);
    ptrdiff_t k = 0;
    for (ptrdiff_t i = 0; i < n_f; ++i) {
        if (in[(size_t)i]) {
            buf->data()[k++] = (idx_t)i;
        }
    }
    return Nodeset::create(fine->comm(), buf);
}

void remap_sets_through_adapt(const std::shared_ptr<Mesh> &coarse,
                              const std::shared_ptr<Mesh> &fine,
                              const ptrdiff_t             *parent_elem,
                              const idx_t                 *node_a,
                              const idx_t                 *node_b) {
    const enum ElemType et   = coarse->element_type(0);
    const int           nxe  = elem_num_nodes(et);
    LocalSideTable      lst;
    LocalEdgeTable      let;
    lst.fill(et);
    let.fill(et);
    const int n_sides = elem_num_sides(et);
    const int n_lei   = elem_num_edges(et);
    const ptrdiff_t n_cnode = coarse->n_nodes();
    const ptrdiff_t n_fnode = fine->n_nodes();
    const ptrdiff_t n_fe    = fine->n_elements(0);
    idx_t **const   fel     = fine->elements(0)->data();
    idx_t **const   cel     = coarse->elements(0)->data();

    auto node_on_side = [&](ptrdiff_t coarse_e, int lfi, idx_t node) -> int {
        const int nv = lst.nnxs_side[lfi];
        std::vector<uint8_t> on((size_t)n_fnode, 0);
        for (int k = 0; k < nv; ++k) {
            const idx_t v = cel[lst(lfi, k)][coarse_e];
            if (v >= 0 && (ptrdiff_t)v < n_cnode) {
                on[(size_t)v] = 1;
            }
        }
        for (ptrdiff_t i = 0; i < n_fnode; ++i) {
            const idx_t a = node_a[i];
            const idx_t b = node_b[i];
            if (a >= 0 && b >= 0 && on[(size_t)a] && on[(size_t)b]) {
                on[(size_t)i] = 1;
            }
        }
        return node >= 0 && (ptrdiff_t)node < n_fnode ? on[(size_t)node] : 0;
    };

    for (const auto &kv : coarse->sidesets()) {
        auto ss = kv.second;
        if (!ss || ss->block_id() != 0) {
            continue;
        }
        std::vector<element_idx_t> par;
        std::vector<i16>           lfi;
        const ptrdiff_t            nss = ss->size();
        const element_idx_t       *sp  = nss > 0 ? ss->parent()->data() : nullptr;
        const i16                 *sl  = nss > 0 ? ss->lfi()->data() : nullptr;
        for (ptrdiff_t i = 0; i < nss; ++i) {
            const ptrdiff_t pe  = (ptrdiff_t)sp[i];
            const int       slf = (int)sl[i];
            if (slf < 0 || slf >= n_sides) {
                continue;
            }
            const int nv = lst.nnxs_side[slf];
            for (ptrdiff_t e = 0; e < n_fe; ++e) {
                if (parent_elem[e] != pe) {
                    continue;
                }
                for (int s = 0; s < n_sides; ++s) {
                    int ok = 1;
                    for (int k = 0; k < nv && k < lst.nnxs_side[s]; ++k) {
                        if (!node_on_side(pe, slf, fel[lst(s, k)][e])) {
                            ok = 0;
                            break;
                        }
                    }
                    if (ok && lst.nnxs_side[s] == nv) {
                        par.push_back((element_idx_t)e);
                        lfi.push_back((i16)s);
                    }
                }
            }
        }
        auto pb = create_host_buffer<element_idx_t>(par.size());
        auto lb = create_host_buffer<i16>(lfi.size());
        for (size_t i = 0; i < par.size(); ++i) {
            pb->data()[i] = par[i];
            lb->data()[i] = lfi[i];
        }
        fine->add_sideset(kv.first, Sideset::create(fine->comm(), pb, lb, 0));
    }

    for (const auto &kv : coarse->edgesets()) {
        auto es = kv.second;
        if (!es || es->block_id() != 0) {
            continue;
        }
        std::vector<element_idx_t> par;
        std::vector<i16>           lei;
        const ptrdiff_t            nes = es->size();
        const element_idx_t       *ep  = nes > 0 ? es->parent()->data() : nullptr;
        const i16                 *el  = nes > 0 ? es->lei()->data() : nullptr;
        for (ptrdiff_t i = 0; i < nes; ++i) {
            const ptrdiff_t pe  = (ptrdiff_t)ep[i];
            const int       sli = (int)el[i];
            if (sli < 0 || sli >= n_lei) {
                continue;
            }
            std::vector<uint8_t> on((size_t)n_fnode, 0);
            {
                const idx_t u = cel[let(sli, 0)][pe];
                const idx_t v = cel[let(sli, 1)][pe];
                if (u >= 0 && (ptrdiff_t)u < n_cnode) {
                    on[(size_t)u] = 1;
                }
                if (v >= 0 && (ptrdiff_t)v < n_cnode) {
                    on[(size_t)v] = 1;
                }
                for (ptrdiff_t k = 0; k < n_fnode; ++k) {
                    const idx_t a = node_a[k];
                    const idx_t b = node_b[k];
                    if (a >= 0 && b >= 0 && on[(size_t)a] && on[(size_t)b]) {
                        on[(size_t)k] = 1;
                    }
                }
            }
            for (ptrdiff_t e = 0; e < n_fe; ++e) {
                if (parent_elem[e] != pe) {
                    continue;
                }
                for (int s = 0; s < n_lei; ++s) {
                    const idx_t a = fel[let(s, 0)][e];
                    const idx_t b = fel[let(s, 1)][e];
                    if (a >= 0 && b >= 0 && on[(size_t)a] && on[(size_t)b]) {
                        par.push_back((element_idx_t)e);
                        lei.push_back((i16)s);
                    }
                }
            }
        }
        auto pb = create_host_buffer<element_idx_t>(par.size());
        auto lb = create_host_buffer<i16>(lei.size());
        for (size_t i = 0; i < par.size(); ++i) {
            pb->data()[i] = par[i];
            lb->data()[i] = lei[i];
        }
        fine->add_edgeset(kv.first, Edgeset::create(fine->comm(), pb, lb, 0));
    }

    for (const auto &kv : coarse->nodesets()) {
        auto mapped = expand_nodeset(coarse, kv.second, fine, node_a, node_b);
        if (mapped) {
            fine->add_nodeset(kv.first, mapped);
        }
    }
    (void)nxe;
}

struct SurfaceSnap {
    enum ElemType et      = TRI3;
    int           nxe     = 0;
    int           sdim    = 0;
    ptrdiff_t     n_elem  = 0;
    ptrdiff_t     n_nodes = 0;
    idx_t       **el      = nullptr;
    geom_t      **p       = nullptr;
};

void surface_snap_free(SurfaceSnap &s) {
    if (s.el) {
        for (int d = 0; d < s.nxe; ++d) {
            SMESH_FREE(s.el[d]);
        }
        SMESH_FREE(s.el);
    }
    if (s.p) {
        for (int d = 0; d < s.sdim; ++d) {
            SMESH_FREE(s.p[d]);
        }
        SMESH_FREE(s.p);
    }
    s.el      = nullptr;
    s.p       = nullptr;
    s.n_elem  = 0;
    s.n_nodes = 0;
    s.nxe     = 0;
    s.sdim    = 0;
}

int surface_snap_copy(const enum ElemType                                      et,
                      const ptrdiff_t                                          n_elem,
                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  el,
                      const int                                                sdim,
                      const ptrdiff_t                                          n_nodes,
                      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT pts,
                      SurfaceSnap                                             *s) {
    if (!s || !el || !pts || n_elem < 0 || n_nodes < 0 || sdim < 2) {
        return SMESH_FAILURE;
    }
    surface_snap_free(*s);
    s->et      = et;
    s->nxe     = elem_num_nodes(et);
    s->sdim    = sdim;
    s->n_elem  = n_elem;
    s->n_nodes = n_nodes;
    s->el      = (idx_t **)SMESH_ALLOC((size_t)s->nxe * sizeof(idx_t *));
    s->p       = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    if (!s->el || !s->p) {
        surface_snap_free(*s);
        return SMESH_FAILURE;
    }
    for (int d = 0; d < s->nxe; ++d) {
        s->el[d] = (idx_t *)SMESH_ALLOC((size_t)n_elem * sizeof(idx_t));
        if (!s->el[d]) {
            s->nxe = d;
            surface_snap_free(*s);
            return SMESH_FAILURE;
        }
        if (n_elem > 0) {
            std::memcpy(s->el[d], el[d], (size_t)n_elem * sizeof(idx_t));
        }
    }
    for (int d = 0; d < sdim; ++d) {
        s->p[d] = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
        if (!s->p[d]) {
            s->sdim = d;
            surface_snap_free(*s);
            return SMESH_FAILURE;
        }
        if (n_nodes > 0) {
            std::memcpy(s->p[d], pts[d], (size_t)n_nodes * sizeof(geom_t));
        }
    }
    return SMESH_SUCCESS;
}

int surface_snap_from_mesh(Mesh &mesh, SurfaceSnap *s) {
    if (!s || mesh.n_blocks() != 1) {
        return SMESH_FAILURE;
    }
    const enum ElemType et   = mesh.element_type(0);
    const int           sdim = mesh.spatial_dimension();
    const ptrdiff_t     nn   = mesh.n_nodes();
    geom_t **const      pts  = mesh.points()->data();
    if (et == TET4) {
        auto view = std::make_shared<Mesh>(mesh.comm(), et, mesh.elements(0), mesh.points());
        auto ss   = skin_sideset(view);
        if (!ss) {
            return SMESH_FAILURE;
        }
        auto created = create_surface_from_sideset(view, ss);
        if (!created.second) {
            return SMESH_FAILURE;
        }
        return surface_snap_copy(created.first,
                                 (ptrdiff_t)created.second->extent(1),
                                 created.second->data(),
                                 sdim,
                                 nn,
                                 pts,
                                 s);
    }
    return surface_snap_copy(et, mesh.n_elements(0), mesh.elements(0)->data(), sdim, nn, pts, s);
}

void mark_nodes_from_faces(const ptrdiff_t                                        n_elem,
                           const int                                              nxe,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT el,
                           const ptrdiff_t                                        n_nodes,
                           uint8_t *const                                         mark) {
    std::memset(mark, 0, (size_t)n_nodes);
    for (ptrdiff_t e = 0; e < n_elem; ++e) {
        for (int d = 0; d < nxe; ++d) {
            const idx_t v = el[d][e];
            if (v >= 0 && (ptrdiff_t)v < n_nodes) {
                mark[v] = 1;
            }
        }
    }
}

int project_onto_snap(Mesh                &mesh,
                      const SurfaceSnap   &snap,
                      const uint8_t       *surf,
                      const uint8_t       *lock,
                      const ptrdiff_t      n_crease,
                      const idx_t         *c0,
                      const idx_t         *c1,
                      const geom_t         max_move) {
    const ptrdiff_t n_nodes = mesh.n_nodes();
    geom_t **const  pts     = mesh.points()->data();
    uint8_t        *pmask   = (uint8_t *)SMESH_ALLOC((size_t)n_nodes);
    uint8_t        *cmask   = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    if (!pmask || !cmask) {
        SMESH_FREE(pmask);
        SMESH_FREE(cmask);
        return SMESH_FAILURE;
    }
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const uint8_t lk = lock ? lock[i] : 0;
        pmask[i]         = (uint8_t)((!surf || surf[i]) && lk == 0);
        cmask[i]         = (uint8_t)(lk == 1);
    }
    int err = mesh_project_to_surface<idx_t, geom_t>(snap.et,
                                                     snap.n_elem,
                                                     snap.el,
                                                     snap.sdim,
                                                     snap.n_nodes,
                                                     snap.p,
                                                     n_nodes,
                                                     pmask,
                                                     pts,
                                                     max_move);
    if (err == SMESH_SUCCESS && n_crease > 0 && c0 && c1) {
        err = mesh_project_to_segments<idx_t, geom_t>(snap.sdim,
                                                      n_crease,
                                                      c0,
                                                      c1,
                                                      snap.n_nodes,
                                                      snap.p,
                                                      n_nodes,
                                                      cmask,
                                                      pts,
                                                      max_move);
    }
    if (err == SMESH_SUCCESS && lock) {
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            if (lock[i] == 2 && (ptrdiff_t)i < snap.n_nodes) {
                for (int d = 0; d < snap.sdim; ++d) {
                    pts[d][i] = snap.p[d][i];
                }
            }
        }
    }
    SMESH_FREE(pmask);
    SMESH_FREE(cmask);
    return err;
}

int constrained_smooth(Mesh              &mesh,
                       const int          n_iters,
                       const geom_t       lambda,
                       const geom_t       sharp_cos,
                       const geom_t       max_abs,
                       const geom_t       max_n,
                       const bool         use_param,
                       const SurfaceSnap *frozen,
                       const ptrdiff_t    n_crease,
                       const idx_t       *crease0,
                       const idx_t       *crease1,
                       const ptrdiff_t    n_pin,
                       const geom_t       q_min) {
    SMESH_TRACE_SCOPE("constrained_smooth");
    const geom_t qbar = q_min > static_cast<geom_t>(0) ? q_min : static_cast<geom_t>(0.5);
    if (n_iters <= 0) {
        return SMESH_SUCCESS;
    }
    SurfaceSnap local;
    const SurfaceSnap *snap = frozen;
    if (!snap) {
        if (surface_snap_from_mesh(mesh, &local) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        snap = &local;
    }

    const int       sdim    = mesh.spatial_dimension();
    const ptrdiff_t n_nodes = mesh.n_nodes();
    geom_t **const  pts     = mesh.points()->data();
    uint8_t        *surf    = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    uint8_t        *lock    = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    geom_t         *nx      = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t         *ny      = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t         *nz      = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    geom_t        **x0      = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    if (!surf || !lock || !nx || !ny || !nz || !x0) {
        SMESH_FREE(surf);
        SMESH_FREE(lock);
        SMESH_FREE(nx);
        SMESH_FREE(ny);
        SMESH_FREE(nz);
        SMESH_FREE(x0);
        surface_snap_free(local);
        return SMESH_FAILURE;
    }
    for (int d = 0; d < sdim; ++d) {
        x0[d] = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
        if (!x0[d]) {
            for (int k = 0; k < d; ++k) {
                SMESH_FREE(x0[k]);
            }
            SMESH_FREE(x0);
            SMESH_FREE(surf);
            SMESH_FREE(lock);
            SMESH_FREE(nx);
            SMESH_FREE(ny);
            SMESH_FREE(nz);
            surface_snap_free(local);
            return SMESH_FAILURE;
        }
        std::memcpy(x0[d], pts[d], (size_t)n_nodes * sizeof(geom_t));
    }

    SurfaceSnap fine_s;
    if (surface_snap_from_mesh(mesh, &fine_s) != SMESH_SUCCESS) {
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(x0[d]);
        }
        SMESH_FREE(x0);
        SMESH_FREE(surf);
        SMESH_FREE(lock);
        SMESH_FREE(nx);
        SMESH_FREE(ny);
        SMESH_FREE(nz);
        surface_snap_free(local);
        return SMESH_FAILURE;
    }
    mark_nodes_from_faces(fine_s.n_elem, fine_s.nxe, fine_s.el, n_nodes, surf);
    idx_t    *se0 = nullptr;
    idx_t    *se1 = nullptr;
    ptrdiff_t nsh = 0;
    {
        auto ebuf = create_host_buffer<idx_t>((size_t)fine_s.nxe, (size_t)fine_s.n_elem);
        for (int d = 0; d < fine_s.nxe; ++d) {
            if (fine_s.n_elem > 0) {
                std::memcpy(ebuf->data()[d], fine_s.el[d], (size_t)fine_s.n_elem * sizeof(idx_t));
            }
        }
        auto sm = std::make_shared<Mesh>(mesh.comm(), fine_s.et, ebuf, mesh.points());
        if (fill_lock_from_sharp(*sm, lock, n_nodes, sharp_cos, &se0, &se1, &nsh) != SMESH_SUCCESS) {
            for (int d = 0; d < sdim; ++d) {
                SMESH_FREE(x0[d]);
            }
            SMESH_FREE(x0);
            SMESH_FREE(surf);
            SMESH_FREE(lock);
            SMESH_FREE(nx);
            SMESH_FREE(ny);
            SMESH_FREE(nz);
            SMESH_FREE(se0);
            SMESH_FREE(se1);
            surface_snap_free(fine_s);
            surface_snap_free(local);
            return SMESH_FAILURE;
        }
        if (mesh_vertex_normals_from_faces<idx_t, geom_t>(
                    fine_s.et, fine_s.n_elem, fine_s.el, sdim, n_nodes, pts, nx, ny, nz) !=
            SMESH_SUCCESS) {
            for (int d = 0; d < sdim; ++d) {
                SMESH_FREE(x0[d]);
            }
            SMESH_FREE(x0);
            SMESH_FREE(surf);
            SMESH_FREE(lock);
            SMESH_FREE(nx);
            SMESH_FREE(ny);
            SMESH_FREE(nz);
            SMESH_FREE(se0);
            SMESH_FREE(se1);
            surface_snap_free(fine_s);
            surface_snap_free(local);
            return SMESH_FAILURE;
        }
    }

    const ptrdiff_t npin = n_pin < n_nodes ? n_pin : n_nodes;
    for (ptrdiff_t i = 0; i < npin; ++i) {
        lock[i] = 2;
    }
    if (n_crease > 0 && crease0 && crease1) {
        for (ptrdiff_t s = 0; s < n_crease; ++s) {
            const idx_t a = crease0[s], b = crease1[s];
            if (a >= 0 && (ptrdiff_t)a < n_nodes && lock[a] < 1) {
                lock[a] = 1;
            }
            if (b >= 0 && (ptrdiff_t)b < n_nodes && lock[b] < 1) {
                lock[b] = 1;
            }
        }
    }

    auto n2n_vol = mesh.node_to_node_graph();
    auto n2n     = n2n_vol;
    auto n2e     = mesh.element_type(0) == TET4 ? mesh.node_to_element_graph() : nullptr;
    uint8_t *lk_int = nullptr;
    if (mesh.element_type(0) == TET4 && n2n && n2e) {
        lk_int = (uint8_t *)SMESH_ALLOC((size_t)n_nodes);
        if (lk_int) {
            std::memcpy(lk_int, lock, (size_t)n_nodes);
            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                if (surf[i]) {
                    lk_int[i] = 2;
                }
            }
            quality_laplace_tets(mesh,
                                 lk_int,
                                 n2n->rowptr()->data(),
                                 n2n->colidx()->data(),
                                 n2e->rowptr()->data(),
                                 n2e->colidx()->data(),
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 n_iters,
                                 lambda);
            sliver_exorcise_tets(mesh,
                                 lk_int,
                                 n2e->rowptr()->data(),
                                 n2e->colidx()->data(),
                                 qbar,
                                 nullptr,
                                 nullptr,
                                 nullptr);
        }
        auto vol = std::make_shared<Mesh>(
                mesh.comm(), mesh.element_type(0), mesh.elements(0), mesh.points());
        auto ss = skin_sideset(vol);
        if (ss) {
            auto created = create_surface_from_sideset(vol, ss);
            if (created.second) {
                auto sm = std::make_shared<Mesh>(
                        mesh.comm(), created.first, created.second, mesh.points());
                n2n = sm->node_to_node_graph();
            }
        }
    }
    int err = SMESH_SUCCESS;
    if (mesh.element_type(0) == TET4 && n2e && n2n) {
        quality_laplace_tets(mesh,
                             lock,
                             n2n->rowptr()->data(),
                             n2n->colidx()->data(),
                             n2e->rowptr()->data(),
                             n2e->colidx()->data(),
                             nx,
                             ny,
                             nz,
                             n_iters,
                             lambda);
        sliver_exorcise_tets(mesh,
                             lock,
                             n2e->rowptr()->data(),
                             n2e->colidx()->data(),
                             qbar,
                             nx,
                             ny,
                             nz);
    } else {
        err = mesh_smooth_feature<idx_t, count_t, geom_t>(sdim,
                                                          n_nodes,
                                                          pts,
                                                          n2n->rowptr()->data(),
                                                          n2n->colidx()->data(),
                                                          lock,
                                                          nsh,
                                                          se0,
                                                          se1,
                                                          n_iters,
                                                          lambda,
                                                          x0,
                                                          surf,
                                                          max_abs,
                                                          max_n,
                                                          nx,
                                                          ny,
                                                          nz);
    }
    SMESH_FREE(nx);
    SMESH_FREE(ny);
    SMESH_FREE(nz);
    if (err != SMESH_SUCCESS) {
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(x0[d]);
        }
        SMESH_FREE(x0);
        SMESH_FREE(se0);
        SMESH_FREE(se1);
        SMESH_FREE(surf);
        SMESH_FREE(lock);
        SMESH_FREE(lk_int);
        surface_snap_free(fine_s);
        surface_snap_free(local);
        return err;
    }

    const idx_t    *pc0 = crease0;
    const idx_t    *pc1 = crease1;
    ptrdiff_t       ncr = (n_crease > 0 && crease0 && crease1) ? n_crease : 0;
    if (ncr == 0 && nsh > 0 && se0 && se1) {
        ncr = nsh;
        pc0 = se0;
        pc1 = se1;
    }
    const geom_t reach =
            max_abs > static_cast<geom_t>(0) ? max_abs * static_cast<geom_t>(2) : static_cast<geom_t>(0);
    const bool snap_pl = !use_param || mesh.parametrizations().empty();
    const int  perr    = snap_pl ? project_onto_snap(mesh, *snap, surf, lock, ncr, pc0, pc1, reach)
                                 : SMESH_SUCCESS;
    SMESH_FREE(se0);
    SMESH_FREE(se1);
    SMESH_FREE(surf);
    SMESH_FREE(lock);
    surface_snap_free(fine_s);
    surface_snap_free(local);
    if (perr != SMESH_SUCCESS) {
        SMESH_FREE(lk_int);
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(x0[d]);
        }
        SMESH_FREE(x0);
        return perr;
    }
    if (use_param) {
        if (n2e) {
            quality_gated_param_apply(mesh, n2e->rowptr()->data(), n2e->colidx()->data());
        } else {
            for (const auto &kv : mesh.parametrizations()) {
                if (kv.second && kv.second->apply(mesh) != SMESH_SUCCESS) {
                    SMESH_FREE(lk_int);
                    for (int d = 0; d < sdim; ++d) {
                        SMESH_FREE(x0[d]);
                    }
                    SMESH_FREE(x0);
                    return SMESH_FAILURE;
                }
            }
        }
    }
    if (lk_int && n2e && n2n_vol) {
        quality_laplace_tets(mesh,
                             lk_int,
                             n2n_vol->rowptr()->data(),
                             n2n_vol->colidx()->data(),
                             n2e->rowptr()->data(),
                             n2e->colidx()->data(),
                             nullptr,
                             nullptr,
                             nullptr,
                             n_iters,
                             lambda);
        sliver_exorcise_tets(mesh,
                             lk_int,
                             n2e->rowptr()->data(),
                             n2e->colidx()->data(),
                             qbar,
                             nullptr,
                             nullptr,
                             nullptr);
    }
    restore_positive_tets(mesh, x0);
    SMESH_FREE(lk_int);
    for (int d = 0; d < sdim; ++d) {
        SMESH_FREE(x0[d]);
    }
    SMESH_FREE(x0);
    return SMESH_SUCCESS;
}

int rebind_params_to_surface(Mesh &mesh) {
    const auto stored = mesh.parametrizations();
    if (stored.empty()) {
        return SMESH_SUCCESS;
    }
    auto view = std::make_shared<Mesh>(
            mesh.comm(), mesh.element_type(0), mesh.elements(0), mesh.points());
    std::shared_ptr<Mesh> surf = view;
    SharedBuffer<idx_t>   map;
    if (mesh.element_type(0) == TET4) {
        surf = skin(view);
        if (!surf) {
            return SMESH_FAILURE;
        }
        map = surf->node_mapping();
    }
    geom_t **const  p  = mesh.points()->data();
    const ptrdiff_t nn = mesh.n_nodes();
    const int       sd = mesh.spatial_dimension();
    mesh.clear_parametrizations();
    for (const auto &kv : stored) {
        if (!kv.second) {
            continue;
        }
        auto sp = std::dynamic_pointer_cast<SphereParametrization>(kv.second);
        if (!sp) {
            auto ns = kv.second->nodeset();
            int  ok = ns != nullptr;
            if (ok) {
                const idx_t *ids = ns->nodes()->data();
                for (ptrdiff_t i = 0; i < ns->size() && ok; ++i) {
                    ok = ids[i] >= 0 && (ptrdiff_t)ids[i] < nn;
                }
            }
            if (ok) {
                mesh.add_parametrization(kv.first, kv.second);
            }
            continue;
        }
        const geom_t lo = static_cast<geom_t>(0.85) * sp->radius();
        std::vector<idx_t> ids;
        auto consider = [&](idx_t v) {
            if (v < 0 || (ptrdiff_t)v >= nn) {
                return;
            }
            geom_t dx = p[0][v] - sp->cx();
            geom_t dy = p[1][v] - sp->cy();
            geom_t dz = sd >= 3 ? p[2][v] - sp->cz() : static_cast<geom_t>(0);
            const geom_t r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r >= lo) {
                ids.push_back(v);
            }
        };
        if (map) {
            const idx_t *mp = map->data();
            for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
                consider(mp[i]);
            }
        } else {
            for (ptrdiff_t i = 0; i < (mesh.element_type(0) == TET4 ? nn : surf->n_nodes()); ++i) {
                consider((idx_t)i);
            }
        }
        if (ids.empty()) {
            continue;
        }
        auto buf = create_host_buffer<idx_t>(ids.size());
        for (size_t i = 0; i < ids.size(); ++i) {
            buf->data()[i] = ids[i];
        }
        mesh.add_parametrization(kv.first, sp->with_nodeset(Nodeset::create(mesh.comm(), buf)));
    }
    return SMESH_SUCCESS;
}

}  // namespace

int smooth_enhance(Mesh &mesh, const AdaptRefineOptions &opt) {
    SMESH_TRACE_SCOPE("smooth_enhance");
    if (mesh.n_blocks() != 1) {
        fprintf(stderr, "smooth_enhance: single-block only\n");
        return SMESH_FAILURE;
    }
    if (mesh.comm() && mesh.comm()->size() > 1) {
        fprintf(stderr, "smooth_enhance: serial only\n");
        return SMESH_FAILURE;
    }
    const geom_t diag = bbox_diag(mesh);
    const geom_t band = static_cast<geom_t>(0.02) * diag;
    return constrained_smooth(mesh,
                              opt.smooth_iters,
                              opt.smooth_lambda,
                              opt.sharp_cos_threshold,
                              band,
                              band,
                              opt.use_parametrization,
                              nullptr,
                              0,
                              nullptr,
                              nullptr,
                              0,
                              opt.q_min);
}

std::shared_ptr<Mesh> adapt_refine(const std::shared_ptr<Mesh> &mesh,
                                   const AdaptRefineOptions    &opt) {
    SMESH_TRACE_SCOPE("adapt_refine");
    if (!mesh) {
        return nullptr;
    }
    if (mesh->n_blocks() != 1) {
        fprintf(stderr, "adapt_refine: single-block only\n");
        return nullptr;
    }
    if (mesh->comm() && mesh->comm()->size() > 1) {
        fprintf(stderr, "adapt_refine: serial only\n");
        return nullptr;
    }
    const enum ElemType et = mesh->element_type(0);
    if (!adapt_refine_type_supported(et)) {
        fprintf(stderr, "adapt_refine: unsupported type %s\n", type_to_string(et));
        return nullptr;
    }

    const int       sdim    = mesh->spatial_dimension();
    const ptrdiff_t n_nodes = mesh->n_nodes();
    const ptrdiff_t n_elem  = mesh->n_elements(0);
    geom_t **const  pts     = mesh->points()->data();
    idx_t **const   elems   = mesh->elements(0)->data();

    enum ElemType surf_et = et;
    idx_t       **surf_el = elems;
    ptrdiff_t     n_surf  = n_elem;
    std::shared_ptr<Buffer<idx_t *>> surf_buf;
    std::shared_ptr<Mesh>            surf_mesh = mesh;
    if (et == TET4) {
        auto ss = skin_sideset(mesh);
        if (!ss) {
            fprintf(stderr, "adapt_refine: TET4 skin failed\n");
            return nullptr;
        }
        auto created = create_surface_from_sideset(mesh, ss);
        surf_et      = created.first;
        surf_buf     = created.second;
        if (!surf_buf) {
            fprintf(stderr, "adapt_refine: TET4 surface extract failed\n");
            return nullptr;
        }
        surf_el = surf_buf->data();
        n_surf  = (ptrdiff_t)surf_buf->extent(1);
        surf_mesh = std::make_shared<Mesh>(mesh->comm(), surf_et, surf_buf, mesh->points());
    }

    geom_t  *kappa = (geom_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(geom_t));
    uint8_t *sharp = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    idx_t   *se0   = nullptr;
    idx_t   *se1   = nullptr;
    ptrdiff_t nsh  = 0;
    fill_lock_from_sharp(*surf_mesh, sharp, n_nodes, opt.sharp_cos_threshold, &se0, &se1, &nsh);
    uint8_t *sharp_node = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        sharp_node[i] = sharp[i] ? 1 : 0;
    }
    mesh_node_curvature<idx_t, geom_t>(
            surf_et, n_surf, surf_el, sdim, n_nodes, pts, nullptr, kappa);

    SurfaceSnap snap;
    if (surface_snap_copy(surf_et, n_surf, surf_el, sdim, n_nodes, pts, &snap) != SMESH_SUCCESS) {
        SMESH_FREE(kappa);
        SMESH_FREE(sharp);
        SMESH_FREE(sharp_node);
        SMESH_FREE(se0);
        SMESH_FREE(se1);
        return nullptr;
    }

    geom_t *h = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    const geom_t diag = bbox_diag(*mesh);
    mesh_size_from_curvature_error(n_nodes,
                                   kappa,
                                   opt.cells_per_radius,
                                   opt.geom_error,
                                   opt.h_min,
                                   opt.h_max,
                                   diag,
                                   h);
    const geom_t pi   = static_cast<geom_t>(3.14159265358979323846);
    const geom_t ncpr = opt.cells_per_radius > static_cast<geom_t>(0) ? opt.cells_per_radius
                                                                     : static_cast<geom_t>(8);
    geom_t       eps  = opt.geom_error;
    if (!(eps > static_cast<geom_t>(0))) {
        eps = (pi * pi) * diag / (static_cast<geom_t>(2) * ncpr * ncpr);
    }
    if (opt.use_parametrization) {
        tighten_h_from_param_error(*surf_mesh, *mesh, h, eps);
    }
    auto n2n = (et == TET4) ? mesh->node_to_node_graph() : surf_mesh->node_to_node_graph();
    mesh_grade_size_2to1<idx_t, count_t, geom_t>(
            n_nodes, n2n->rowptr()->data(), n2n->colidx()->data(), 8, h);

    ptrdiff_t  n_elem_o = 0, n_node_o = 0;
    idx_t    **elems_o  = nullptr;
    geom_t   **pts_o    = nullptr;
    ptrdiff_t *parent   = nullptr;
    count_t   *pptr     = nullptr;
    idx_t     *cid      = nullptr;
    idx_t     *na       = nullptr;
    idx_t     *nb       = nullptr;
    const int nxe       = elem_num_nodes(et);
    uint8_t *geom_node = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    if (geom_node) {
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            geom_node[i] = (kappa[i] * diag > static_cast<geom_t>(0.5)) ? 1 : 0;
        }
    }
    if (mesh_adapt_refine<idx_t, count_t, geom_t>(et,
                                                  n_elem,
                                                  elems,
                                                  sdim,
                                                  n_nodes,
                                                  pts,
                                                  h,
                                                  opt.element_mark,
                                                  opt.q_min,
                                                  geom_node,
                                                  opt.max_levels,
                                                  &n_elem_o,
                                                  &elems_o,
                                                  &n_node_o,
                                                  &pts_o,
                                                  &parent,
                                                  &pptr,
                                                  &cid,
                                                  &na,
                                                  &nb) != SMESH_SUCCESS) {
        SMESH_FREE(kappa);
        SMESH_FREE(sharp);
        SMESH_FREE(sharp_node);
        SMESH_FREE(h);
        SMESH_FREE(geom_node);
        SMESH_FREE(se0);
        SMESH_FREE(se1);
        surface_snap_free(snap);
        return nullptr;
    }

    auto fine = wrap_adapt_mesh(mesh, et, nxe, sdim, n_elem_o, elems_o, n_node_o, pts_o);
    remap_sets_through_adapt(mesh, fine, parent, na, nb);

    if (opt.use_parametrization) {
        for (const auto &kv : mesh->parametrizations()) {
            if (!kv.second) {
                continue;
            }
            auto src = kv.second->nodeset();
            auto exp = expand_nodeset(mesh, src, fine, na, nb);
            auto np  = kv.second->with_nodeset(exp ? exp : src);
            fine->add_parametrization(kv.first, np);
            np->apply(*fine);
        }
        pull_mids_of_inverted_tets(*fine, na, nb);
    }

    mesh_adapt_refine_free(nxe, sdim, elems_o, pts_o, parent, pptr, cid, na, nb);
    SMESH_FREE(kappa);
    SMESH_FREE(sharp);
    SMESH_FREE(sharp_node);
    SMESH_FREE(h);
    SMESH_FREE(geom_node);

    if (opt.smooth_iters > 0) {
        const geom_t band = static_cast<geom_t>(0.02) * diag;
        if (constrained_smooth(*fine,
                               opt.smooth_iters,
                               opt.smooth_lambda,
                               opt.sharp_cos_threshold,
                               band,
                               band,
                               opt.use_parametrization,
                               &snap,
                               nsh,
                               se0,
                               se1,
                               n_nodes,
                               opt.q_min) != SMESH_SUCCESS) {
            SMESH_FREE(se0);
            SMESH_FREE(se1);
            surface_snap_free(snap);
            return nullptr;
        }
    }
    SMESH_FREE(se0);
    SMESH_FREE(se1);
    surface_snap_free(snap);
    return fine;
}

static void install_improve_soA(Mesh       &mesh,
                         const int   nxe,
                         const int   sdim,
                         ptrdiff_t   n_elem,
                         idx_t     **elems,
                         ptrdiff_t   n_nodes,
                         geom_t    **pts) {
    auto ebuf = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_elem);
    auto pbuf = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_nodes);
    for (int d = 0; d < nxe; ++d) {
        std::memcpy(ebuf->data()[d], elems[d], (size_t)n_elem * sizeof(idx_t));
    }
    for (int d = 0; d < sdim; ++d) {
        std::memcpy(pbuf->data()[d], pts[d], (size_t)n_nodes * sizeof(geom_t));
    }
    mesh.set_points(pbuf);
    mesh.block(0)->set_elements(ebuf);
}

int improve(Mesh &mesh, const ImproveOptions &opt) {
    SMESH_TRACE_SCOPE("improve");
    if (mesh.n_blocks() != 1) {
        fprintf(stderr, "improve: single-block only\n");
        return SMESH_FAILURE;
    }
    if (mesh.comm() && mesh.comm()->size() > 1) {
        fprintf(stderr, "improve: serial only\n");
        return SMESH_FAILURE;
    }
    const enum ElemType et = mesh.element_type(0);
    if (!improve_type_supported(et)) {
        fprintf(stderr, "improve: unsupported type %s\n", type_to_string(et));
        return SMESH_FAILURE;
    }

    const int       sdim    = mesh.spatial_dimension();
    const int       nxe     = elem_num_nodes(et);
    const ptrdiff_t n_nodes = mesh.n_nodes();
    const ptrdiff_t n_elem  = mesh.n_elements(0);
    geom_t **const  pts     = mesh.points()->data();
    idx_t **const   elems   = mesh.elements(0)->data();

    geom_t **x0 = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    for (int d = 0; d < sdim; ++d) {
        x0[d] = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
        std::memcpy(x0[d], pts[d], (size_t)n_nodes * sizeof(geom_t));
    }
    uint8_t *lock = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    uint8_t *surf = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));

    std::shared_ptr<Mesh> surf_mesh = nullptr;
    std::shared_ptr<Buffer<idx_t *>> surf_buf;
    auto vol_view = std::make_shared<Mesh>(mesh.comm(), et, mesh.elements(0), mesh.points());
    if (et == TET4) {
        auto skin_ss = skin_sideset(vol_view);
        if (!skin_ss) {
            fprintf(stderr, "improve: TET4 skin failed\n");
            for (int d = 0; d < sdim; ++d) {
                SMESH_FREE(x0[d]);
            }
            SMESH_FREE(x0);
            SMESH_FREE(lock);
            SMESH_FREE(surf);
            return SMESH_FAILURE;
        }
        auto created = create_surface_from_sideset(vol_view, skin_ss);
        surf_buf     = created.second;
        if (!surf_buf) {
            fprintf(stderr, "improve: TET4 surface extract failed\n");
            for (int d = 0; d < sdim; ++d) {
                SMESH_FREE(x0[d]);
            }
            SMESH_FREE(x0);
            SMESH_FREE(lock);
            SMESH_FREE(surf);
            return SMESH_FAILURE;
        }
        const ptrdiff_t n_surf = (ptrdiff_t)surf_buf->extent(1);
        idx_t **const   sel    = surf_buf->data();
        const int       snxe   = elem_num_nodes(created.first);
        for (ptrdiff_t e = 0; e < n_surf; ++e) {
            for (int d = 0; d < snxe; ++d) {
                const idx_t v = sel[d][e];
                if (v >= 0 && (ptrdiff_t)v < n_nodes) {
                    surf[v] = 1;
                }
            }
        }
        surf_mesh = std::make_shared<Mesh>(mesh.comm(), created.first, surf_buf, mesh.points());
    } else {
        surf_mesh = vol_view;
        std::memset(surf, 1, (size_t)n_nodes);
    }

    idx_t    *e0  = nullptr;
    idx_t    *e1  = nullptr;
    ptrdiff_t nsh = 0;
    fill_lock_from_sharp(*surf_mesh, lock, n_nodes, opt.sharp_cos_threshold, &e0, &e1, &nsh);

    const geom_t diag = bbox_diag(mesh);
    geom_t       max_abs = opt.max_abs_dev;
    geom_t       max_n   = opt.max_normal_dev;
    if (!(max_abs < static_cast<geom_t>(0)) && max_abs == static_cast<geom_t>(0)) {
        max_abs = static_cast<geom_t>(0.02) * diag;
    }
    if (!(max_n < static_cast<geom_t>(0)) && max_n == static_cast<geom_t>(0)) {
        max_n = max_abs;
    }

    ptrdiff_t  n_elem_o = 0, n_node_o = 0, n_ops = 0;
    idx_t    **elems_o = nullptr;
    geom_t   **pts_o   = nullptr;
    uint8_t   *lock_o  = nullptr;
    uint8_t   *surf_o  = nullptr;
    geom_t   **x0_o    = nullptr;
    if (mesh_improve<idx_t, count_t, geom_t>(et,
                                             n_elem,
                                             elems,
                                             sdim,
                                             n_nodes,
                                             pts,
                                             lock,
                                             surf,
                                             x0,
                                             nsh,
                                             e0,
                                             e1,
                                             opt.q_min,
                                             max_abs,
                                             max_n,
                                             opt.max_passes,
                                             opt.allow_split ? 1 : 0,
                                             opt.allow_collapse ? 1 : 0,
                                             opt.allow_swap ? 1 : 0,
                                             &n_elem_o,
                                             &elems_o,
                                             &n_node_o,
                                             &pts_o,
                                             &lock_o,
                                             &surf_o,
                                             &x0_o,
                                             &n_ops) != SMESH_SUCCESS) {
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(x0[d]);
        }
        SMESH_FREE(x0);
        SMESH_FREE(lock);
        SMESH_FREE(surf);
        SMESH_FREE(e0);
        SMESH_FREE(e1);
        return SMESH_FAILURE;
    }

    install_improve_soA(mesh, nxe, sdim, n_elem_o, elems_o, n_node_o, pts_o);
    mesh.invalidate_derived_graphs();
    if (rebind_params_to_surface(mesh) != SMESH_SUCCESS) {
        mesh_improve_free(nxe, sdim, elems_o, pts_o, lock_o, surf_o, x0_o);
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(x0[d]);
        }
        SMESH_FREE(x0);
        SMESH_FREE(lock);
        SMESH_FREE(surf);
        SMESH_FREE(e0);
        SMESH_FREE(e1);
        return SMESH_FAILURE;
    }

    geom_t max_abs_s = max_abs;
    if (opt.use_parametrization && !(opt.max_abs_dev > static_cast<geom_t>(0))) {
        max_abs_s = static_cast<geom_t>(0.04) * diag;
    }

    if (opt.smooth_iters > 0) {
        if (constrained_smooth(mesh,
                               opt.smooth_iters,
                               opt.smooth_lambda,
                               opt.sharp_cos_threshold,
                               max_abs_s,
                               max_abs_s,
                               opt.use_parametrization,
                               nullptr,
                               0,
                               nullptr,
                               nullptr,
                               0,
                               opt.q_min) != SMESH_SUCCESS) {
            mesh_improve_free(nxe, sdim, elems_o, pts_o, lock_o, surf_o, x0_o);
            for (int d = 0; d < sdim; ++d) {
                SMESH_FREE(x0[d]);
            }
            SMESH_FREE(x0);
            SMESH_FREE(lock);
            SMESH_FREE(surf);
            SMESH_FREE(e0);
            SMESH_FREE(e1);
            return SMESH_FAILURE;
        }
        if (et == TET4 && opt.max_passes > 1 && lock_o && surf_o) {
            geom_t    qprev = static_cast<geom_t>(-1);
            ptrdiff_t nprev = -1;
            for (int extra = 0; extra < 6; ++extra) {
                const geom_t    qm = mesh_tet_qmin(mesh);
                const ptrdiff_t nb = mesh_tet_nbelow(mesh, opt.q_min);
                if (nb == 0) {
                    break;
                }
                if (extra > 0 && !(nb < nprev) && !(qm > qprev)) {
                    break;
                }
                qprev = qm;
                nprev = nb;
                const ptrdiff_t nn1     = mesh.n_nodes();
                const ptrdiff_t ne1     = mesh.n_elements(0);
                geom_t **const  pc      = mesh.points()->data();
                idx_t **const   ec      = mesh.elements(0)->data();
                idx_t         **esnap   = (idx_t **)SMESH_ALLOC((size_t)nxe * sizeof(idx_t *));
                geom_t        **psnap   = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
                int             snapok  = esnap && psnap;
                if (esnap) {
                    for (int d = 0; d < nxe; ++d) {
                        esnap[d] = nullptr;
                    }
                }
                if (psnap) {
                    for (int d = 0; d < sdim; ++d) {
                        psnap[d] = nullptr;
                    }
                }
                for (int d = 0; d < nxe && snapok; ++d) {
                    esnap[d] = (idx_t *)SMESH_ALLOC((size_t)ne1 * sizeof(idx_t));
                    if (!esnap[d]) {
                        snapok = 0;
                        break;
                    }
                    std::memcpy(esnap[d], ec[d], (size_t)ne1 * sizeof(idx_t));
                }
                for (int d = 0; d < sdim && snapok; ++d) {
                    psnap[d] = (geom_t *)SMESH_ALLOC((size_t)nn1 * sizeof(geom_t));
                    if (!psnap[d]) {
                        snapok = 0;
                        break;
                    }
                    std::memcpy(psnap[d], pc[d], (size_t)nn1 * sizeof(geom_t));
                }
                auto free_snap = [&]() {
                    if (esnap) {
                        for (int d = 0; d < nxe; ++d) {
                            SMESH_FREE(esnap[d]);
                        }
                        SMESH_FREE(esnap);
                    }
                    if (psnap) {
                        for (int d = 0; d < sdim; ++d) {
                            SMESH_FREE(psnap[d]);
                        }
                        SMESH_FREE(psnap);
                    }
                };
                if (!snapok) {
                    free_snap();
                    break;
                }
                ptrdiff_t n_elem2 = 0, n_node2 = 0, n_ops2 = 0;
                idx_t   **elems2 = nullptr;
                geom_t  **pts2   = nullptr;
                uint8_t  *lock2  = nullptr;
                uint8_t  *surf2  = nullptr;
                geom_t  **x0_2   = nullptr;
                const int ierr2  = mesh_improve<idx_t, count_t, geom_t>(et,
                                                                            ne1,
                                                                            mesh.elements(0)->data(),
                                                                            sdim,
                                                                            nn1,
                                                                            pc,
                                                                            lock_o,
                                                                            surf_o,
                                                                            psnap,
                                                                            0,
                                                                            nullptr,
                                                                            nullptr,
                                                                            opt.q_min,
                                                                            max_abs,
                                                                            max_n,
                                                                            4,
                                                                            opt.allow_split ? 1 : 0,
                                                                            opt.allow_collapse ? 1 : 0,
                                                                            opt.allow_swap ? 1 : 0,
                                                                            &n_elem2,
                                                                            &elems2,
                                                                            &n_node2,
                                                                            &pts2,
                                                                            &lock2,
                                                                            &surf2,
                                                                            &x0_2,
                                                                            &n_ops2);
                if (ierr2 != SMESH_SUCCESS) {
                    mesh_improve_free(nxe, sdim, elems2, pts2, lock2, surf2, x0_2);
                    free_snap();
                    break;
                }
                install_improve_soA(mesh, nxe, sdim, n_elem2, elems2, n_node2, pts2);
                mesh.invalidate_derived_graphs();
                const int r2 = rebind_params_to_surface(mesh);
                const int s2 = r2 == SMESH_SUCCESS
                                       ? constrained_smooth(mesh,
                                                            opt.smooth_iters,
                                                            opt.smooth_lambda,
                                                            opt.sharp_cos_threshold,
                                                            max_abs_s,
                                                            max_abs_s,
                                                            opt.use_parametrization,
                                                            nullptr,
                                                            0,
                                                            nullptr,
                                                            nullptr,
                                                            0,
                                                            opt.q_min)
                                       : SMESH_FAILURE;
                if (s2 != SMESH_SUCCESS) {
                    mesh_improve_free(nxe, sdim, elems2, pts2, lock2, surf2, x0_2);
                    free_snap();
                    mesh_improve_free(nxe, sdim, elems_o, pts_o, lock_o, surf_o, x0_o);
                    for (int d = 0; d < sdim; ++d) {
                        SMESH_FREE(x0[d]);
                    }
                    SMESH_FREE(x0);
                    SMESH_FREE(lock);
                    SMESH_FREE(surf);
                    SMESH_FREE(e0);
                    SMESH_FREE(e1);
                    return SMESH_FAILURE;
                }
                if (mesh_tet_nbelow(mesh, opt.q_min) > nprev ||
                    (mesh_tet_nbelow(mesh, opt.q_min) == nprev && mesh_tet_qmin(mesh) < qm)) {
                    install_improve_soA(mesh, nxe, sdim, ne1, esnap, nn1, psnap);
                    mesh.invalidate_derived_graphs();
                    mesh_improve_free(nxe, sdim, elems2, pts2, lock2, surf2, x0_2);
                    free_snap();
                    break;
                }
                free_snap();
                mesh_improve_free(nxe, sdim, elems_o, pts_o, lock_o, surf_o, x0_o);
                elems_o = elems2;
                pts_o   = pts2;
                lock_o  = lock2;
                surf_o  = surf2;
                x0_o    = x0_2;
                n_ops += n_ops2;
            }
        }
    } else if (opt.use_parametrization) {
        for (const auto &kv : mesh.parametrizations()) {
            if (kv.second && kv.second->apply(mesh) != SMESH_SUCCESS) {
                mesh_improve_free(nxe, sdim, elems_o, pts_o, lock_o, surf_o, x0_o);
                for (int d = 0; d < sdim; ++d) {
                    SMESH_FREE(x0[d]);
                }
                SMESH_FREE(x0);
                SMESH_FREE(lock);
                SMESH_FREE(surf);
                SMESH_FREE(e0);
                SMESH_FREE(e1);
                return SMESH_FAILURE;
            }
        }
    }

    mesh_improve_free(nxe, sdim, elems_o, pts_o, lock_o, surf_o, x0_o);
    for (int d = 0; d < sdim; ++d) {
        SMESH_FREE(x0[d]);
    }
    SMESH_FREE(x0);
    SMESH_FREE(lock);
    SMESH_FREE(surf);
    SMESH_FREE(e0);
    SMESH_FREE(e1);
    (void)n_ops;
    return SMESH_SUCCESS;
}

std::shared_ptr<Mesh> remesh(const std::shared_ptr<Mesh> &mesh, const ImproveOptions &opt) {
    SMESH_TRACE_SCOPE("remesh");
    if (!mesh) {
        return nullptr;
    }
    auto out = mesh->clone();
    if (improve(*out, opt) != SMESH_SUCCESS) {
        return nullptr;
    }
    return out;
}

}  // namespace smesh
