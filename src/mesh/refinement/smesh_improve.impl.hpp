#ifndef SMESH_IMPROVE_IMPL_HPP
#define SMESH_IMPROVE_IMPL_HPP

#include "smesh_improve.hpp"

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_graph.hpp"
#include "smesh_quality.hpp"
#include "smesh_search.hpp"

#include <cmath>
#include <string.h>

namespace smesh {

namespace {

static const int kQuad4Child[4][4] = {{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};

template <typename T>
T *irealloc(T *p, const ptrdiff_t oldn, const ptrdiff_t newn) {
    T *q = (T *)SMESH_ALLOC((size_t)newn * sizeof(T));
    if (!q) {
        return nullptr;
    }
    if (p && oldn > 0) {
        memcpy(q, p, (size_t)(oldn < newn ? oldn : newn) * sizeof(T));
    }
    if (newn > oldn) {
        memset(q + oldn, 0, (size_t)(newn - oldn) * sizeof(T));
    }
    SMESH_FREE(p);
    return q;
}

template <typename idx_t, typename count_t>
count_t find_n2n_slot(const count_t *const SMESH_RESTRICT rowptr,
                      const idx_t *const SMESH_RESTRICT   colidx,
                      idx_t                               a,
                      idx_t                               b) {
    if (a > b) {
        const idx_t t = a;
        a             = b;
        b             = t;
    }
    const count_t beg = rowptr[a];
    const count_t len = rowptr[a + 1] - beg;
    const count_t k   = binary_search(b, colidx + beg, len);
    if (k < len && colidx[beg + k] == b) {
        return beg + k;
    }
    return (count_t)(-1);
}

template <typename idx_t, typename geom_t>
geom_t edge_len2(const int                                                sdim,
                 const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p,
                 const idx_t                                              a,
                 const idx_t                                              b) {
    geom_t s = 0;
    for (int d = 0; d < sdim; ++d) {
        const geom_t t = p[d][b] - p[d][a];
        s += t * t;
    }
    return s;
}

template <typename idx_t, typename geom_t>
geom_t tet_orient(const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p,
                  const idx_t                                              a,
                  const idx_t                                              b,
                  const idx_t                                              c,
                  const idx_t                                              d) {
    const geom_t ux = p[0][b] - p[0][a], uy = p[1][b] - p[1][a], uz = p[2][b] - p[2][a];
    const geom_t vx = p[0][c] - p[0][a], vy = p[1][c] - p[1][a], vz = p[2][c] - p[2][a];
    const geom_t wx = p[0][d] - p[0][a], wy = p[1][d] - p[1][a], wz = p[2][d] - p[2][a];
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx);
}

template <typename idx_t, typename geom_t>
void fix_tet(const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p, idx_t *v) {
    if (tet_orient(p, v[0], v[1], v[2], v[3]) < static_cast<geom_t>(0)) {
        const idx_t t = v[2];
        v[2]          = v[3];
        v[3]          = t;
    }
}

template <typename idx_t, typename geom_t>
void fix_tri(const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p,
             const int                                                sdim,
             idx_t                                                   *v) {
    const geom_t ux = p[0][v[1]] - p[0][v[0]], uy = p[1][v[1]] - p[1][v[0]];
    const geom_t vx = p[0][v[2]] - p[0][v[0]], vy = p[1][v[2]] - p[1][v[0]];
    geom_t       o  = ux * vy - uy * vx;
    if (sdim >= 3) {
        const geom_t uz = p[2][v[1]] - p[2][v[0]], vz = p[2][v[2]] - p[2][v[0]];
        const geom_t nz = ux * vy - uy * vx;
        (void)uz;
        (void)vz;
        o = nz;
    }
    if (sdim < 3 && o < static_cast<geom_t>(0)) {
        const idx_t t = v[1];
        v[1]          = v[2];
        v[2]          = t;
    }
}

template <typename idx_t, typename geom_t>
geom_t q_from(const enum ElemType                                     et,
              const int                                               sdim,
              const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p,
              idx_t                                                   a,
              idx_t                                                   b,
              idx_t                                                   c,
              idx_t                                                   d) {
    idx_t *soa[4] = {&a, &b, &c, &d};
    return mesh_elem_mean_ratio<idx_t, geom_t>(et, sdim, soa, p, 0);
}

template <typename idx_t, typename geom_t>
struct ImproveW {
    enum ElemType et;
    int           sdim, nxe, n_lei, n_sides;
    int           is_tri, is_tet, is_quad;
    LocalEdgeTable let;
    LocalSideTable lst;
    ptrdiff_t      n_elem, n_nodes, cap_elem, cap_nodes;
    idx_t        **elems;
    geom_t       **pts;
    geom_t       **x0;
    uint8_t       *lock;
    uint8_t       *surf;
    uint8_t       *dead_e;
    geom_t         q_min, max_abs, max_n;
    ptrdiff_t      n_sharp;
    idx_t         *se0;
    idx_t         *se1;
};

template <typename idx_t, typename geom_t>
int grow_elem(ImproveW<idx_t, geom_t> &w, const ptrdiff_t need) {
    if (need <= w.cap_elem) {
        return SMESH_SUCCESS;
    }
    ptrdiff_t cap = w.cap_elem > 0 ? w.cap_elem : 8;
    while (cap < need) {
        cap *= 2;
    }
    for (int d = 0; d < w.nxe; ++d) {
        w.elems[d] = irealloc(w.elems[d], w.cap_elem, cap);
        if (!w.elems[d]) {
            return SMESH_FAILURE;
        }
    }
    w.dead_e = irealloc(w.dead_e, w.cap_elem, cap);
    if (!w.dead_e) {
        return SMESH_FAILURE;
    }
    w.cap_elem = cap;
    return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int grow_node(ImproveW<idx_t, geom_t> &w, const ptrdiff_t need) {
    if (need <= w.cap_nodes) {
        return SMESH_SUCCESS;
    }
    ptrdiff_t cap = w.cap_nodes > 0 ? w.cap_nodes : 8;
    while (cap < need) {
        cap *= 2;
    }
    for (int d = 0; d < w.sdim; ++d) {
        w.pts[d] = irealloc(w.pts[d], w.cap_nodes, cap);
        w.x0[d]  = irealloc(w.x0[d], w.cap_nodes, cap);
        if (!w.pts[d] || !w.x0[d]) {
            return SMESH_FAILURE;
        }
    }
    w.lock = irealloc(w.lock, w.cap_nodes, cap);
    w.surf = irealloc(w.surf, w.cap_nodes, cap);
    if (!w.lock || !w.surf) {
        return SMESH_FAILURE;
    }
    w.cap_nodes = cap;
    return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
void clamp_node(ImproveW<idx_t, geom_t> &w, const ptrdiff_t i) {
    if (w.lock[i] == 2) {
        for (int d = 0; d < w.sdim; ++d) {
            w.pts[d][i] = w.x0[d][i];
        }
        return;
    }
    if (!w.surf[i] && w.is_tet) {
        return;
    }
    if (!(w.max_abs > static_cast<geom_t>(0))) {
        return;
    }
    geom_t x = w.pts[0][i], y = w.pts[1][i], z = w.sdim >= 3 ? w.pts[2][i] : static_cast<geom_t>(0);
    mesh_clamp_to_band(w.sdim,
                       &x,
                       &y,
                       w.sdim >= 3 ? &z : nullptr,
                       w.x0[0][i],
                       w.x0[1][i],
                       w.sdim >= 3 ? w.x0[2][i] : static_cast<geom_t>(0),
                       static_cast<geom_t>(0),
                       static_cast<geom_t>(0),
                       static_cast<geom_t>(0),
                       w.max_abs,
                       static_cast<geom_t>(0));
    w.pts[0][i] = x;
    w.pts[1][i] = y;
    if (w.sdim >= 3) {
        w.pts[2][i] = z;
    }
}

template <typename idx_t, typename geom_t>
idx_t add_mid(ImproveW<idx_t, geom_t> &w, const idx_t a, const idx_t b, const int feature) {
    if (grow_node(w, w.n_nodes + 1) != SMESH_SUCCESS) {
        return static_cast<idx_t>(-1);
    }
    const idx_t m = (idx_t)w.n_nodes;
    for (int d = 0; d < w.sdim; ++d) {
        w.pts[d][m] = static_cast<geom_t>(0.5) * (w.pts[d][a] + w.pts[d][b]);
        w.x0[d][m]  = static_cast<geom_t>(0.5) * (w.x0[d][a] + w.x0[d][b]);
    }
    uint8_t lk = 0;
    if (feature) {
        lk = 1;
        if (w.lock[a] == 2 && w.lock[b] == 2) {
            lk = 1;
        }
    }
    w.lock[m] = lk;
    w.surf[m] = (uint8_t)((w.surf[a] && w.surf[b]) || (!w.is_tet));
    if (feature) {
        idx_t *n0 = irealloc(w.se0, w.n_sharp, w.n_sharp + 2);
        idx_t *n1 = irealloc(w.se1, w.n_sharp, w.n_sharp + 2);
        if (n0 && n1) {
            w.se0 = n0;
            w.se1 = n1;
            w.se0[w.n_sharp]     = a;
            w.se1[w.n_sharp]     = m;
            w.se0[w.n_sharp + 1] = m;
            w.se1[w.n_sharp + 1] = b;
            w.n_sharp += 2;
        }
    }
    if (w.lock[m] == 2) {
        for (int d = 0; d < w.sdim; ++d) {
            w.pts[d][m] = w.x0[d][m];
        }
    } else {
        clamp_node(w, (ptrdiff_t)m);
    }
    w.n_nodes += 1;
    return m;
}

template <typename idx_t, typename geom_t>
int write_elem(ImproveW<idx_t, geom_t> &w, const ptrdiff_t e, const idx_t *v) {
    if (e >= w.n_elem) {
        if (grow_elem(w, e + 1) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        w.n_elem = e + 1;
    }
    for (int d = 0; d < w.nxe; ++d) {
        w.elems[d][e] = v[d];
    }
    w.dead_e[e] = 0;
    return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int compact_w(ImproveW<idx_t, geom_t> &w) {
    ptrdiff_t we = 0;
    for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
        if (w.dead_e[e]) {
            continue;
        }
        if (we != e) {
            for (int d = 0; d < w.nxe; ++d) {
                w.elems[d][we] = w.elems[d][e];
            }
        }
        w.dead_e[we] = 0;
        ++we;
    }
    w.n_elem = we;
    uint8_t *used = (uint8_t *)SMESH_CALLOC((size_t)w.n_nodes, sizeof(uint8_t));
    if (!used) {
        return SMESH_FAILURE;
    }
    for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
        for (int d = 0; d < w.nxe; ++d) {
            const idx_t v = w.elems[d][e];
            if (v >= 0 && (ptrdiff_t)v < w.n_nodes) {
                used[v] = 1;
            }
        }
    }
    idx_t *map = (idx_t *)SMESH_ALLOC((size_t)w.n_nodes * sizeof(idx_t));
    if (!map) {
        SMESH_FREE(used);
        return SMESH_FAILURE;
    }
    ptrdiff_t wn = 0;
    for (ptrdiff_t i = 0; i < w.n_nodes; ++i) {
        if (!used[i]) {
            map[i] = static_cast<idx_t>(-1);
            continue;
        }
        map[i] = (idx_t)wn;
        if (wn != i) {
            for (int d = 0; d < w.sdim; ++d) {
                w.pts[d][wn] = w.pts[d][i];
                w.x0[d][wn]  = w.x0[d][i];
            }
            w.lock[wn] = w.lock[i];
            w.surf[wn] = w.surf[i];
        }
        ++wn;
    }
    for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
        for (int d = 0; d < w.nxe; ++d) {
            w.elems[d][e] = map[w.elems[d][e]];
        }
    }
    w.n_nodes = wn;
    ptrdiff_t ns = 0;
    for (ptrdiff_t s = 0; s < w.n_sharp; ++s) {
        const idx_t a = w.se0[s] >= 0 ? map[w.se0[s]] : static_cast<idx_t>(-1);
        const idx_t b = w.se1[s] >= 0 ? map[w.se1[s]] : static_cast<idx_t>(-1);
        if (a < 0 || b < 0) {
            continue;
        }
        w.se0[ns] = a;
        w.se1[ns] = b;
        ++ns;
    }
    w.n_sharp = ns;
    SMESH_FREE(used);
    SMESH_FREE(map);
    return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
void unique_edges(const ptrdiff_t n_nodes,
                  const count_t  *rowptr,
                  const idx_t    *colidx,
                  ptrdiff_t      *n_uedge,
                  ptrdiff_t     **uid_out,
                  idx_t         **eu_out,
                  idx_t         **ev_out) {
    ptrdiff_t nu = 0;
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
            if ((ptrdiff_t)colidx[k] > i) {
                ++nu;
            }
        }
    }
    ptrdiff_t *uid = (ptrdiff_t *)SMESH_ALLOC((size_t)rowptr[n_nodes] * sizeof(ptrdiff_t));
    idx_t     *eu  = (idx_t *)SMESH_ALLOC((size_t)nu * sizeof(idx_t));
    idx_t     *ev  = (idx_t *)SMESH_ALLOC((size_t)nu * sizeof(idx_t));
    ptrdiff_t  id  = 0;
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
            uid[k] = -1;
            if ((ptrdiff_t)colidx[k] > i) {
                uid[k] = id;
                eu[id] = (idx_t)i;
                ev[id] = colidx[k];
                ++id;
            }
        }
    }
    *n_uedge = nu;
    *uid_out = uid;
    *eu_out  = eu;
    *ev_out  = ev;
}

template <typename idx_t, typename count_t>
int edge_incidents(const count_t       *n2eptr,
                   const element_idx_t *elindex,
                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                   const int                                               nxe,
                   const idx_t                                             a,
                   const idx_t                                             b,
                   const uint8_t                                          *dead,
                   ptrdiff_t                                              *inc,
                   const int                                               max_inc) {
    int n = 0;
    for (count_t k = n2eptr[a]; k < n2eptr[a + 1]; ++k) {
        const ptrdiff_t e = (ptrdiff_t)elindex[k];
        if (e < 0 || (dead && dead[e])) {
            continue;
        }
        int             ha = 0, hb = 0;
        for (int d = 0; d < nxe; ++d) {
            ha |= elems[d][e] == a;
            hb |= elems[d][e] == b;
        }
        if (ha && hb) {
            if (n >= max_inc) {
                return max_inc + 1;
            }
            inc[n++] = e;
        }
    }
    return n;
}

template <typename idx_t, typename geom_t>
geom_t qe(ImproveW<idx_t, geom_t> &w, const ptrdiff_t e) {
    return mesh_elem_mean_ratio<idx_t, geom_t>(w.et, w.sdim, w.elems, w.pts, e);
}

template <typename idx_t, typename count_t, typename geom_t>
ptrdiff_t tri_flips(ImproveW<idx_t, geom_t> &w,
                    const count_t           *rowptr,
                    const idx_t             *colidx,
                    const ptrdiff_t         *uid,
                    const idx_t             *eu,
                    const idx_t             *ev,
                    const uint8_t           *feat,
                    const ptrdiff_t          n_uedge,
                    const count_t           *n2eptr,
                    const element_idx_t             *elindex,
                    uint8_t                 *busy) {
    ptrdiff_t nops = 0;
    for (ptrdiff_t id = 0; id < n_uedge; ++id) {
        if (feat[id]) {
            continue;
        }
        const idx_t a = eu[id], b = ev[id];
        if (busy[a] || busy[b] || w.lock[a] == 2 || w.lock[b] == 2) {
            continue;
        }
        ptrdiff_t inc[4];
        const int ni =
                edge_incidents<idx_t, count_t>(n2eptr, elindex, w.elems, w.nxe, a, b, w.dead_e, inc, 4);
        if (ni != 2) {
            continue;
        }
        idx_t c = static_cast<idx_t>(-1), d = static_cast<idx_t>(-1);
        for (int k = 0; k < 3; ++k) {
            const idx_t v = w.elems[k][inc[0]];
            if (v != a && v != b) {
                c = v;
            }
        }
        for (int k = 0; k < 3; ++k) {
            const idx_t v = w.elems[k][inc[1]];
            if (v != a && v != b) {
                d = v;
            }
        }
        if (c < 0 || d < 0 || c == d || busy[c] || busy[d]) {
            continue;
        }
        if (find_n2n_slot(rowptr, colidx, c, d) != (count_t)(-1)) {
            continue;
        }
        const geom_t qold = qe(w, inc[0]) < qe(w, inc[1]) ? qe(w, inc[0]) : qe(w, inc[1]);
        idx_t        t0[3] = {a, c, d};
        idx_t        t1[3] = {b, d, c};
        fix_tri(w.pts, w.sdim, t0);
        fix_tri(w.pts, w.sdim, t1);
        const geom_t q0 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, t0[0], t0[1], t0[2], 0);
        const geom_t q1 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, t1[0], t1[1], t1[2], 0);
        const geom_t qn = q0 < q1 ? q0 : q1;
        if (!(qn > qold)) {
            continue;
        }
        write_elem(w, inc[0], t0);
        write_elem(w, inc[1], t1);
        busy[a] = busy[b] = busy[c] = busy[d] = 1;
        ++nops;
    }
    (void)uid;
    return nops;
}

template <typename idx_t, typename count_t, typename geom_t>
ptrdiff_t edge_splits(ImproveW<idx_t, geom_t> &w,
                      const count_t           *rowptr,
                      const idx_t             *colidx,
                      const ptrdiff_t         *uid,
                      const idx_t             *eu,
                      const idx_t             *ev,
                      const uint8_t           *feat,
                      const ptrdiff_t          n_uedge,
                      const count_t           *n2eptr,
                      const element_idx_t             *elindex,
                      const geom_t            *q,
                      uint8_t                 *busy) {
    ptrdiff_t nops = 0;
    for (ptrdiff_t id = 0; id < n_uedge; ++id) {
        const idx_t a = eu[id], b = ev[id];
        if (busy[a] || busy[b]) {
            continue;
        }
        ptrdiff_t inc[16];
        const int ni = edge_incidents<idx_t, count_t>(
                n2eptr, elindex, w.elems, w.nxe, a, b, w.dead_e, inc, 16);
        if (ni < 1 || ni > 12) {
            continue;
        }
            int do_split = 0;
        const geom_t l2 = edge_len2(w.sdim, w.pts, a, b);
        for (int k = 0; k < ni; ++k) {
            const ptrdiff_t e = inc[k];
            if (!(q[e] < w.q_min)) {
                continue;
            }
            int longest = 1;
            for (int le = 0; le < w.n_lei; ++le) {
                const idx_t u = w.elems[w.let(le, 0)][e];
                const idx_t v = w.elems[w.let(le, 1)][e];
                if (edge_len2(w.sdim, w.pts, u, v) > l2 * static_cast<geom_t>(1.0000001)) {
                    longest = 0;
                }
            }
            do_split |= longest;
        }
        if (!do_split) {
            continue;
        }
        if (w.is_tet && ni == 1) {
            /* boundary edge with a single tet is allowed */
        }
        geom_t qold = q[inc[0]];
        for (int k = 1; k < ni; ++k) {
            if (q[inc[k]] < qold) {
                qold = q[inc[k]];
            }
        }
        const ptrdiff_t n_before  = w.n_nodes;
        const ptrdiff_t sh_before = w.n_sharp;
        const idx_t     m         = add_mid(w, a, b, feat[id]);
        if (m < 0) {
            break;
        }
        geom_t qnmin = static_cast<geom_t>(1);
        for (int k = 0; k < ni; ++k) {
            const ptrdiff_t e = inc[k];
            idx_t           v[8];
            for (int d = 0; d < w.nxe; ++d) {
                v[d] = w.elems[d][e];
            }
            if (w.is_tri) {
                idx_t other = 0;
                for (int d = 0; d < 3; ++d) {
                    if (v[d] != a && v[d] != b) {
                        other = v[d];
                    }
                }
                idx_t c0[3] = {a, m, other};
                idx_t c1[3] = {m, b, other};
                fix_tri(w.pts, w.sdim, c0);
                fix_tri(w.pts, w.sdim, c1);
                const geom_t q0 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c0[0], c0[1], c0[2], 0);
                const geom_t q1 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c1[0], c1[1], c1[2], 0);
                const geom_t qm = q0 < q1 ? q0 : q1;
                if (qm < qnmin) {
                    qnmin = qm;
                }
            } else {
                idx_t o0 = 0, o1 = 0, nv = 0;
                for (int d = 0; d < 4; ++d) {
                    if (v[d] != a && v[d] != b) {
                        if (nv == 0) {
                            o0 = v[d];
                        } else {
                            o1 = v[d];
                        }
                        ++nv;
                    }
                }
                idx_t c0[4] = {a, m, o0, o1};
                idx_t c1[4] = {m, b, o0, o1};
                fix_tet(w.pts, c0);
                fix_tet(w.pts, c1);
                const geom_t q0 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c0[0], c0[1], c0[2], c0[3]);
                const geom_t q1 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c1[0], c1[1], c1[2], c1[3]);
                const geom_t qm = q0 < q1 ? q0 : q1;
                if (qm < qnmin) {
                    qnmin = qm;
                }
            }
        }
        if (!(qnmin > qold)) {
            w.n_nodes = n_before;
            w.n_sharp = sh_before;
            continue;
        }
        for (int k = 0; k < ni; ++k) {
            const ptrdiff_t e = inc[k];
            idx_t           v[8];
            for (int d = 0; d < w.nxe; ++d) {
                v[d] = w.elems[d][e];
            }
            if (w.is_tri) {
                idx_t other = 0;
                for (int d = 0; d < 3; ++d) {
                    if (v[d] != a && v[d] != b) {
                        other = v[d];
                    }
                }
                idx_t c0[3] = {a, m, other};
                idx_t c1[3] = {m, b, other};
                fix_tri(w.pts, w.sdim, c0);
                fix_tri(w.pts, w.sdim, c1);
                write_elem(w, e, c0);
                write_elem(w, w.n_elem, c1);
            } else {
                idx_t o0 = 0, o1 = 0, nv = 0;
                for (int d = 0; d < 4; ++d) {
                    if (v[d] != a && v[d] != b) {
                        if (nv == 0) {
                            o0 = v[d];
                        } else {
                            o1 = v[d];
                        }
                        ++nv;
                    }
                }
                idx_t c0[4] = {a, m, o0, o1};
                idx_t c1[4] = {m, b, o0, o1};
                fix_tet(w.pts, c0);
                fix_tet(w.pts, c1);
                write_elem(w, e, c0);
                write_elem(w, w.n_elem, c1);
            }
        }
        busy[a] = busy[b] = 1;
        ++nops;
    }
    (void)rowptr;
    (void)colidx;
    (void)uid;
    return nops;
}

template <typename idx_t, typename count_t, typename geom_t>
ptrdiff_t edge_collapses(ImproveW<idx_t, geom_t> &w,
                         const count_t           *rowptr,
                         const idx_t             *colidx,
                         const idx_t             *eu,
                         const idx_t             *ev,
                         const uint8_t           *feat,
                         const ptrdiff_t          n_uedge,
                         const count_t           *n2eptr,
                         const element_idx_t             *elindex,
                         uint8_t                 *busy) {
    ptrdiff_t nops = 0;
    for (ptrdiff_t id = 0; id < n_uedge; ++id) {
        const idx_t a0 = eu[id], b0 = ev[id];
        if (busy[a0] || busy[b0]) {
            continue;
        }
        if (w.lock[a0] == 2 && w.lock[b0] == 2) {
            continue;
        }
        if (feat[id]) {
            if (!(w.lock[a0] == 1 && w.lock[b0] == 1)) {
                continue;
            }
        }
        idx_t survive = a0, killed = b0;
        if (w.lock[b0] > w.lock[a0]) {
            survive = b0;
            killed  = a0;
        } else if (w.lock[a0] == w.lock[b0] && b0 < a0) {
            survive = b0;
            killed  = a0;
        }
        if (w.lock[killed] == 2) {
            continue;
        }
        if (w.is_tet && w.surf[a0] != w.surf[b0]) {
            continue;
        }
        ptrdiff_t einc[64];
        const int ni = edge_incidents<idx_t, count_t>(
                n2eptr, elindex, w.elems, w.nxe, a0, b0, w.dead_e, einc, 64);
        if (ni < 1 || ni > 64) {
            continue;
        }
        if (w.is_tri) {
            int ncommon = 0;
            for (count_t k = rowptr[a0]; k < rowptr[a0 + 1]; ++k) {
                const idx_t nb = colidx[k];
                if (nb == b0) {
                    continue;
                }
                if (find_n2n_slot(rowptr, colidx, b0, nb) != (count_t)(-1)) {
                    ++ncommon;
                }
            }
            if (ni == 2 && ncommon != 2) {
                continue;
            }
            if (ni == 1 && ncommon != 1) {
                continue;
            }
        }
        if (w.is_tet) {
            idx_t ring[32];
            int   nr          = 0;
            int   n_surf_ring = 0;
            int   ring_full   = 0;
            for (int k = 0; k < ni; ++k) {
                for (int d = 0; d < 4; ++d) {
                    const idx_t v = w.elems[d][einc[k]];
                    if (v == a0 || v == b0) {
                        continue;
                    }
                    int seen = 0;
                    for (int t = 0; t < nr; ++t) {
                        seen |= ring[t] == v;
                    }
                    if (seen) {
                        continue;
                    }
                    if (nr >= 32) {
                        ring_full = 1;
                        break;
                    }
                    ring[nr++] = v;
                    n_surf_ring += w.surf[v] ? 1 : 0;
                }
            }
            if (ring_full) {
                continue;
            }
            idx_t star[64];
            int   ns = 0;
            int   star_full = 0;
            for (count_t k = n2eptr[a0]; k < n2eptr[a0 + 1] && !star_full; ++k) {
                const ptrdiff_t ea = (ptrdiff_t)elindex[k];
                if (ea < 0 || w.dead_e[ea]) {
                    continue;
                }
                for (int d = 0; d < 4; ++d) {
                    const idx_t v = w.elems[d][ea];
                    if (v == a0) {
                        continue;
                    }
                    int seen = 0;
                    for (int t = 0; t < ns; ++t) {
                        seen |= star[t] == v;
                    }
                    if (seen) {
                        continue;
                    }
                    if (ns >= 64) {
                        star_full = 1;
                        break;
                    }
                    star[ns++] = v;
                }
            }
            if (star_full) {
                continue;
            }
            int ncom = 0;
            for (int t = 0; t < ns; ++t) {
                const idx_t v = star[t];
                if (v == b0) {
                    continue;
                }
                int in_b = 0;
                for (count_t k2 = n2eptr[b0]; k2 < n2eptr[b0 + 1] && !in_b; ++k2) {
                    const ptrdiff_t eb = (ptrdiff_t)elindex[k2];
                    if (eb < 0 || w.dead_e[eb]) {
                        continue;
                    }
                    for (int d2 = 0; d2 < 4; ++d2) {
                        in_b |= w.elems[d2][eb] == v;
                    }
                }
                ncom += in_b ? 1 : 0;
            }
            if (ncom != nr) {
                continue;
            }
            if (w.surf[a0] && w.surf[b0] && n_surf_ring != 2) {
                continue;
            }
        }
        geom_t qold = qe(w, einc[0]);
        for (int k = 1; k < ni; ++k) {
            const geom_t qq = qe(w, einc[k]);
            if (qq < qold) {
                qold = qq;
            }
        }
        ptrdiff_t cav[128];
        int       ncav = 0;
        for (count_t k = n2eptr[killed]; k < n2eptr[killed + 1]; ++k) {
            const ptrdiff_t e = (ptrdiff_t)elindex[k];
            if (w.dead_e[e]) {
                continue;
            }
            int skip = 0;
            for (int j = 0; j < ni; ++j) {
                if (einc[j] == e) {
                    skip = 1;
                }
            }
            if (skip) {
                continue;
            }
            if (ncav < 128) {
                cav[ncav++] = e;
            }
        }
        for (int k = 0; k < ncav; ++k) {
            const geom_t qq = qe(w, cav[k]);
            if (qq < qold) {
                qold = qq;
            }
        }
        int ok = 1;
        geom_t qnmin = static_cast<geom_t>(1);
        for (int k = 0; k < ncav; ++k) {
            idx_t v[8];
            int   dup = 0;
            for (int d = 0; d < w.nxe; ++d) {
                v[d] = w.elems[d][cav[k]] == killed ? survive : w.elems[d][cav[k]];
            }
            for (int d = 0; d < w.nxe; ++d) {
                for (int j = d + 1; j < w.nxe; ++j) {
                    if (v[d] == v[j]) {
                        dup = 1;
                    }
                }
            }
            if (dup) {
                ok = 0;
                break;
            }
            if (w.is_tet) {
                fix_tet(w.pts, v);
            } else if (w.is_tri) {
                fix_tri(w.pts, w.sdim, v);
            }
            const geom_t qq =
                    q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, v[0], v[1], v[2], w.is_tet ? v[3] : 0);
            if (qq < qnmin) {
                qnmin = qq;
            }
            if (!(qq > static_cast<geom_t>(0))) {
                ok = 0;
                break;
            }
        }
        if (!ok || (ncav > 0 && !(qnmin > qold))) {
            continue;
        }
        if (ncav == 0) {
            continue;
        }
        for (int k = 0; k < ncav; ++k) {
            idx_t v[8];
            for (int d = 0; d < w.nxe; ++d) {
                v[d] = w.elems[d][cav[k]] == killed ? survive : w.elems[d][cav[k]];
            }
            if (w.is_tet) {
                fix_tet(w.pts, v);
            } else if (w.is_tri) {
                fix_tri(w.pts, w.sdim, v);
            }
            write_elem(w, cav[k], v);
        }
        for (int k = 0; k < ni; ++k) {
            w.dead_e[einc[k]] = 1;
        }
        busy[survive] = 1;
        busy[killed]  = 1;
        ++nops;
    }
    return nops;
}

template <typename idx_t, typename count_t, typename geom_t>
ptrdiff_t tet_23(ImproveW<idx_t, geom_t> &w,
                 const element_idx_t     *adj,
                 const count_t           *rowptr,
                 const idx_t             *colidx,
                 const ptrdiff_t         *uid,
                 const uint8_t           *feat,
                 uint8_t                 *busy) {
    ptrdiff_t nops = 0;
    const int ns   = w.n_sides;
    const ptrdiff_t n0 = w.n_elem;
    for (ptrdiff_t e = 0; e < n0; ++e) {
        if (w.dead_e[e]) {
            continue;
        }
        for (int f = 0; f < ns; ++f) {
            const element_idx_t nb = adj[e * ns + f];
            if (nb == invalid_idx<element_idx_t>() || (ptrdiff_t)nb <= e || (ptrdiff_t)nb >= n0 ||
                w.dead_e[(ptrdiff_t)nb]) {
                continue;
            }
            idx_t face[3];
            for (int k = 0; k < 3; ++k) {
                face[k] = w.elems[w.lst(f, k)][e];
            }
            if (feat && rowptr && colidx && uid) {
                int on_feat = 0;
                for (int k = 0; k < 3; ++k) {
                    const count_t slot = find_n2n_slot(rowptr, colidx, face[k], face[(k + 1) % 3]);
                    if (slot != (count_t)(-1) && uid[slot] >= 0 && feat[uid[slot]]) {
                        on_feat = 1;
                        break;
                    }
                }
                if (on_feat) {
                    continue;
                }
            }
            idx_t te = static_cast<idx_t>(-1), tn = static_cast<idx_t>(-1);
            for (int d = 0; d < 4; ++d) {
                const idx_t v = w.elems[d][e];
                if (v != face[0] && v != face[1] && v != face[2]) {
                    te = v;
                }
            }
            for (int d = 0; d < 4; ++d) {
                const idx_t v = w.elems[d][(ptrdiff_t)nb];
                if (v != face[0] && v != face[1] && v != face[2]) {
                    tn = v;
                }
            }
            if (te < 0 || tn < 0) {
                continue;
            }
            if (rowptr && colidx && find_n2n_slot(rowptr, colidx, te, tn) != (count_t)(-1)) {
                continue;
            }
            const int nsurf = (w.surf[face[0]] ? 1 : 0) + (w.surf[face[1]] ? 1 : 0) +
                              (w.surf[face[2]] ? 1 : 0) + (w.surf[te] ? 1 : 0) + (w.surf[tn] ? 1 : 0);
            if (nsurf >= 4) {
                continue;
            }
            const idx_t vs[5] = {face[0], face[1], face[2], te, tn};
            int         hit   = 0;
            for (int k = 0; k < 5; ++k) {
                hit |= busy[vs[k]];
            }
            if (hit) {
                continue;
            }
            const geom_t qold =
                    qe(w, e) < qe(w, (ptrdiff_t)nb) ? qe(w, e) : qe(w, (ptrdiff_t)nb);
            idx_t c0[4] = {face[0], face[1], te, tn};
            idx_t c1[4] = {face[1], face[2], te, tn};
            idx_t c2[4] = {face[2], face[0], te, tn};
            fix_tet(w.pts, c0);
            fix_tet(w.pts, c1);
            fix_tet(w.pts, c2);
            const geom_t q0 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c0[0], c0[1], c0[2], c0[3]);
            const geom_t q1 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c1[0], c1[1], c1[2], c1[3]);
            const geom_t q2 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c2[0], c2[1], c2[2], c2[3]);
            geom_t       qn = q0 < q1 ? q0 : q1;
            if (q2 < qn) {
                qn = q2;
            }
            if (!(qn > qold)) {
                continue;
            }
            write_elem(w, e, c0);
            write_elem(w, (ptrdiff_t)nb, c1);
            write_elem(w, w.n_elem, c2);
            for (int k = 0; k < 5; ++k) {
                busy[vs[k]] = 1;
            }
            ++nops;
        }
    }
    return nops;
}

template <typename idx_t, typename count_t, typename geom_t>
ptrdiff_t tet_32(ImproveW<idx_t, geom_t> &w,
                 const idx_t             *eu,
                 const idx_t             *ev,
                 const uint8_t           *feat,
                 const ptrdiff_t          n_uedge,
                 const count_t           *n2eptr,
                         const element_idx_t             *elindex,
                 uint8_t                 *busy) {
    ptrdiff_t nops = 0;
    for (ptrdiff_t id = 0; id < n_uedge; ++id) {
        if (feat[id]) {
            continue;
        }
        const idx_t a = eu[id], b = ev[id];
        if (busy[a] || busy[b] || w.lock[a] == 2 || w.lock[b] == 2) {
            continue;
        }
        if (w.surf[a] || w.surf[b]) {
            continue;
        }
        ptrdiff_t inc[8];
        const int ni =
                edge_incidents<idx_t, count_t>(n2eptr, elindex, w.elems, w.nxe, a, b, w.dead_e, inc, 8);
        if (ni != 3) {
            continue;
        }
        idx_t eq[3];
        int   neq = 0;
        for (int k = 0; k < 3; ++k) {
            for (int d = 0; d < 4; ++d) {
                const idx_t v = w.elems[d][inc[k]];
                if (v == a || v == b) {
                    continue;
                }
                int seen = 0;
                for (int t = 0; t < neq; ++t) {
                    seen |= eq[t] == v;
                }
                if (!seen && neq < 3) {
                    eq[neq++] = v;
                }
            }
        }
        if (neq != 3) {
            continue;
        }
        if (busy[eq[0]] || busy[eq[1]] || busy[eq[2]]) {
            continue;
        }
        geom_t qold = qe(w, inc[0]);
        for (int k = 1; k < 3; ++k) {
            const geom_t qq = qe(w, inc[k]);
            if (qq < qold) {
                qold = qq;
            }
        }
        idx_t c0[4] = {eq[0], eq[1], eq[2], a};
        idx_t c1[4] = {eq[0], eq[1], eq[2], b};
        fix_tet(w.pts, c0);
        fix_tet(w.pts, c1);
        const geom_t q0 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c0[0], c0[1], c0[2], c0[3]);
        const geom_t q1 = q_from<idx_t, geom_t>(w.et, w.sdim, w.pts, c1[0], c1[1], c1[2], c1[3]);
        const geom_t qn = q0 < q1 ? q0 : q1;
        if (!(qn > qold)) {
            continue;
        }
        write_elem(w, inc[0], c0);
        write_elem(w, inc[1], c1);
        w.dead_e[inc[2]] = 1;
        busy[a] = busy[b] = busy[eq[0]] = busy[eq[1]] = busy[eq[2]] = 1;
        ++nops;
    }
    return nops;
}

template <typename idx_t, typename count_t, typename geom_t>
ptrdiff_t quad_split(ImproveW<idx_t, geom_t> &w,
                     const count_t           *rowptr,
                     const idx_t             *colidx,
                     const ptrdiff_t         *uid,
                     const ptrdiff_t          n_uedge,
                     const geom_t            *q) {
    uint8_t *emark = (uint8_t *)SMESH_CALLOC((size_t)w.n_elem, sizeof(uint8_t));
    uint8_t *eflag = (uint8_t *)SMESH_CALLOC((size_t)n_uedge, sizeof(uint8_t));
    if (!emark || !eflag) {
        SMESH_FREE(emark);
        SMESH_FREE(eflag);
        return 0;
    }
    for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
        emark[e] = q[e] < w.q_min ? 1 : 0;
    }
    int changed = 1;
    while (changed) {
        changed = 0;
        for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
            if (!emark[e]) {
                continue;
            }
            for (int le = 0; le < 4; ++le) {
                const count_t slot = find_n2n_slot(rowptr, colidx, w.elems[w.let(le, 0)][e],
                                                   w.elems[w.let(le, 1)][e]);
                if (slot != (count_t)(-1) && uid[slot] >= 0) {
                    eflag[uid[slot]] = 1;
                }
            }
        }
        for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
            if (emark[e]) {
                continue;
            }
            for (int le = 0; le < 4; ++le) {
                const count_t slot = find_n2n_slot(rowptr, colidx, w.elems[w.let(le, 0)][e],
                                                   w.elems[w.let(le, 1)][e]);
                if (slot != (count_t)(-1) && uid[slot] >= 0 && eflag[uid[slot]]) {
                    emark[e] = 1;
                    changed  = 1;
                    break;
                }
            }
        }
    }
    ptrdiff_t *mid_of = (ptrdiff_t *)SMESH_ALLOC((size_t)n_uedge * sizeof(ptrdiff_t));
    for (ptrdiff_t i = 0; i < n_uedge; ++i) {
        mid_of[i] = -1;
    }
    ptrdiff_t n_split = 0;
    for (ptrdiff_t e = 0; e < w.n_elem; ++e) {
        n_split += emark[e];
    }
    if (n_split == 0) {
        SMESH_FREE(emark);
        SMESH_FREE(eflag);
        SMESH_FREE(mid_of);
        return 0;
    }
    for (ptrdiff_t i = 0; i < n_uedge; ++i) {
        if (!eflag[i]) {
            continue;
        }
        /* endpoints recovered when applying */
    }
    for (ptrdiff_t i = 0; i < w.n_nodes; ++i) {
        for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
            if ((ptrdiff_t)colidx[k] <= i) {
                continue;
            }
            const ptrdiff_t id = uid[k];
            if (id < 0 || !eflag[id]) {
                continue;
            }
            const idx_t a = (idx_t)i, b = colidx[k];
            const int   feature = w.lock[a] >= 1 && w.lock[b] >= 1;
            const idx_t m       = add_mid(w, a, b, feature);
            mid_of[id]          = (ptrdiff_t)m;
        }
    }
    const ptrdiff_t n_old = w.n_elem;
    for (ptrdiff_t e = 0; e < n_old; ++e) {
        if (!emark[e]) {
            continue;
        }
        idx_t macro[9];
        for (int k = 0; k < 4; ++k) {
            macro[k] = w.elems[k][e];
        }
        for (int le = 0; le < 4; ++le) {
            const count_t slot =
                    find_n2n_slot(rowptr, colidx, w.elems[w.let(le, 0)][e], w.elems[w.let(le, 1)][e]);
            macro[4 + le] = (idx_t)mid_of[uid[slot]];
        }
        if (grow_node(w, w.n_nodes + 1) != SMESH_SUCCESS) {
            break;
        }
        const idx_t fc = (idx_t)w.n_nodes;
        for (int d = 0; d < w.sdim; ++d) {
            w.pts[d][fc] = static_cast<geom_t>(0.25) *
                           (w.pts[d][macro[0]] + w.pts[d][macro[1]] + w.pts[d][macro[2]] +
                            w.pts[d][macro[3]]);
            w.x0[d][fc] = static_cast<geom_t>(0.25) *
                          (w.x0[d][macro[0]] + w.x0[d][macro[1]] + w.x0[d][macro[2]] + w.x0[d][macro[3]]);
        }
        w.lock[fc] = 0;
        w.surf[fc] = 1;
        clamp_node(w, (ptrdiff_t)fc);
        w.n_nodes += 1;
        macro[8] = fc;
        idx_t c0[4];
        for (int d = 0; d < 4; ++d) {
            c0[d] = macro[kQuad4Child[0][d]];
        }
        write_elem(w, e, c0);
        for (int c = 1; c < 4; ++c) {
            idx_t ch[4];
            for (int d = 0; d < 4; ++d) {
                ch[d] = macro[kQuad4Child[c][d]];
            }
            write_elem(w, w.n_elem, ch);
        }
    }
    SMESH_FREE(emark);
    SMESH_FREE(eflag);
    SMESH_FREE(mid_of);
    return n_split;
}

template <typename idx_t, typename count_t, typename geom_t>
int rebuild_improve_graphs(ImproveW<idx_t, geom_t> &w,
                           count_t                 **n2eptr,
                           element_idx_t           **elindex,
                           count_t                 **rowptr,
                           idx_t                   **colidx,
                           ptrdiff_t                *n_uedge,
                           ptrdiff_t               **uid,
                           idx_t                   **eu,
                           idx_t                   **ev,
                           uint8_t                 **feat) {
    SMESH_FREE(*n2eptr);
    SMESH_FREE(*elindex);
    SMESH_FREE(*rowptr);
    SMESH_FREE(*colidx);
    SMESH_FREE(*uid);
    SMESH_FREE(*eu);
    SMESH_FREE(*ev);
    SMESH_FREE(*feat);
    *n2eptr = nullptr;
    *elindex = nullptr;
    *rowptr = nullptr;
    *colidx = nullptr;
    *uid = nullptr;
    *eu = nullptr;
    *ev = nullptr;
    *feat = nullptr;
    *n_uedge = 0;
    if (create_n2e<idx_t, count_t, element_idx_t>(
                w.n_elem, w.n_nodes, w.nxe, w.elems, n2eptr, elindex) != SMESH_SUCCESS ||
        create_edge_graph_for_element_from_n2e<idx_t, count_t>(
                w.et, w.n_elem, w.n_nodes, w.elems, *n2eptr, *elindex, rowptr, colidx) !=
                SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    unique_edges<idx_t, count_t>(w.n_nodes, *rowptr, *colidx, n_uedge, uid, eu, ev);
    *feat = (uint8_t *)SMESH_CALLOC((size_t)(*n_uedge > 0 ? *n_uedge : 1), sizeof(uint8_t));
    for (ptrdiff_t s = 0; s < w.n_sharp; ++s) {
        const count_t slot = find_n2n_slot(*rowptr, *colidx, w.se0[s], w.se1[s]);
        if (slot != (count_t)(-1) && *uid && (*uid)[slot] >= 0) {
            (*feat)[(*uid)[slot]] = 1;
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace

template <typename idx_t, typename geom_t>
void improve_free_t(const int nxe, const int sdim, idx_t **elements, geom_t **points, uint8_t *lock,
                    uint8_t *surface, geom_t **x0) {
    if (elements) {
        for (int d = 0; d < nxe; ++d) {
            SMESH_FREE(elements[d]);
        }
        SMESH_FREE(elements);
    }
    if (points) {
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(points[d]);
        }
        SMESH_FREE(points);
    }
    SMESH_FREE(lock);
    SMESH_FREE(surface);
    if (x0) {
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(x0[d]);
        }
        SMESH_FREE(x0);
    }
}

void mesh_improve_free(const int nxe, const int sdim, idx_t **elements, geom_t **points, uint8_t *lock,
                       uint8_t *surface, geom_t **x0) {
    improve_free_t<idx_t, geom_t>(nxe, sdim, elements, points, lock, surface, x0);
}

template <typename idx_t, typename count_t, typename geom_t>
int mesh_improve(const enum ElemType                                      element_type,
                 const ptrdiff_t                                          n_elements,
                 const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                 const int                                                sdim,
                 const ptrdiff_t                                          n_nodes,
                 const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                 const uint8_t *const SMESH_RESTRICT                      lock,
                 const uint8_t *const SMESH_RESTRICT                      surface,
                 const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT x0,
                 const ptrdiff_t                                          n_sharp,
                 const idx_t *const SMESH_RESTRICT                        se0,
                 const idx_t *const SMESH_RESTRICT                        se1,
                 const geom_t                                             q_min,
                 const geom_t                                             max_abs_dev,
                 const geom_t                                             max_normal_dev,
                 const int                                                max_passes,
                 const int                                                allow_split,
                 const int                                                allow_collapse,
                 const int                                                allow_swap,
                 ptrdiff_t                                               *n_elements_out,
                 idx_t                                                 ***elements_out,
                 ptrdiff_t                                               *n_nodes_out,
                 geom_t                                                ***points_out,
                 uint8_t                                                **lock_out,
                 uint8_t                                                **surface_out,
                 geom_t                                                ***x0_out,
                 ptrdiff_t                                               *n_ops_out) {
    if (!improve_type_supported(element_type) || !elements || !points || !lock || !x0 || !n_elements_out ||
        !elements_out || !n_nodes_out || !points_out) {
        return SMESH_FAILURE;
    }
    ImproveW<idx_t, geom_t> w;
    memset(&w, 0, sizeof(w));
    w.et      = element_type;
    w.sdim    = sdim;
    w.nxe     = elem_num_nodes(element_type);
    w.n_lei   = elem_num_edges(element_type);
    w.n_sides = elem_num_sides(element_type);
    w.is_tri  = element_type == TRI3 || element_type == TRISHELL3;
    w.is_tet  = element_type == TET4;
    w.is_quad = element_type == QUAD4 || element_type == QUADSHELL4;
    w.q_min   = q_min;
    w.max_abs = max_abs_dev;
    w.max_n   = max_normal_dev;
    if (w.let.fill(element_type) != SMESH_SUCCESS || w.lst.fill(element_type) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    w.cap_elem  = n_elements > 0 ? n_elements : 1;
    w.cap_nodes = n_nodes > 0 ? n_nodes : 1;
    w.n_elem    = n_elements;
    w.n_nodes   = n_nodes;
    w.elems     = (idx_t **)SMESH_ALLOC((size_t)w.nxe * sizeof(idx_t *));
    w.pts       = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    w.x0        = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    w.dead_e    = (uint8_t *)SMESH_CALLOC((size_t)w.cap_elem, sizeof(uint8_t));
    w.lock      = (uint8_t *)SMESH_ALLOC((size_t)w.cap_nodes * sizeof(uint8_t));
    w.surf      = (uint8_t *)SMESH_ALLOC((size_t)w.cap_nodes * sizeof(uint8_t));
    for (int d = 0; d < w.nxe; ++d) {
        w.elems[d] = (idx_t *)SMESH_ALLOC((size_t)w.cap_elem * sizeof(idx_t));
        memcpy(w.elems[d], elements[d], (size_t)n_elements * sizeof(idx_t));
    }
    for (int d = 0; d < sdim; ++d) {
        w.pts[d] = (geom_t *)SMESH_ALLOC((size_t)w.cap_nodes * sizeof(geom_t));
        w.x0[d]  = (geom_t *)SMESH_ALLOC((size_t)w.cap_nodes * sizeof(geom_t));
        memcpy(w.pts[d], points[d], (size_t)n_nodes * sizeof(geom_t));
        memcpy(w.x0[d], x0[d], (size_t)n_nodes * sizeof(geom_t));
    }
    memcpy(w.lock, lock, (size_t)n_nodes);
    if (surface) {
        memcpy(w.surf, surface, (size_t)n_nodes);
    } else {
        memset(w.surf, 1, (size_t)n_nodes);
    }
    w.n_sharp = n_sharp > 0 ? n_sharp : 0;
    w.se0     = w.n_sharp > 0 ? (idx_t *)SMESH_ALLOC((size_t)w.n_sharp * sizeof(idx_t)) : nullptr;
    w.se1     = w.n_sharp > 0 ? (idx_t *)SMESH_ALLOC((size_t)w.n_sharp * sizeof(idx_t)) : nullptr;
    if (w.n_sharp > 0) {
        memcpy(w.se0, se0, (size_t)w.n_sharp * sizeof(idx_t));
        memcpy(w.se1, se1, (size_t)w.n_sharp * sizeof(idx_t));
    }

    ptrdiff_t total_ops = 0;
    const int npass     = max_passes > 0 ? max_passes : 1;
    for (int pass = 0; pass < npass; ++pass) {
        if (compact_w(w) != SMESH_SUCCESS) {
            SMESH_FREE(w.se0);
            SMESH_FREE(w.se1);
            improve_free_t<idx_t, geom_t>(w.nxe, w.sdim, w.elems, w.pts, w.lock, w.surf, w.x0);
            return SMESH_FAILURE;
        }
        count_t       *n2eptr  = nullptr;
        element_idx_t *elindex = nullptr;
        count_t       *rowptr  = nullptr;
        idx_t         *colidx  = nullptr;
        if (create_n2e<idx_t, count_t, element_idx_t>(
                    w.n_elem, w.n_nodes, w.nxe, w.elems, &n2eptr, &elindex) != SMESH_SUCCESS ||
            create_edge_graph_for_element_from_n2e<idx_t, count_t>(
                    w.et, w.n_elem, w.n_nodes, w.elems, n2eptr, elindex, &rowptr, &colidx) !=
                    SMESH_SUCCESS) {
            SMESH_FREE(n2eptr);
            SMESH_FREE(elindex);
            SMESH_FREE(w.se0);
            SMESH_FREE(w.se1);
            improve_free_t<idx_t, geom_t>(w.nxe, w.sdim, w.elems, w.pts, w.lock, w.surf, w.x0);
            return SMESH_FAILURE;
        }
        ptrdiff_t  n_uedge = 0;
        ptrdiff_t *uid     = nullptr;
        idx_t     *eu = nullptr, *ev = nullptr;
        unique_edges<idx_t, count_t>(w.n_nodes, rowptr, colidx, &n_uedge, &uid, &eu, &ev);
        uint8_t *feat = (uint8_t *)SMESH_CALLOC((size_t)(n_uedge > 0 ? n_uedge : 1), sizeof(uint8_t));
        for (ptrdiff_t s = 0; s < w.n_sharp; ++s) {
            const count_t slot = find_n2n_slot(rowptr, colidx, w.se0[s], w.se1[s]);
            if (slot != (count_t)(-1) && uid[slot] >= 0) {
                feat[uid[slot]] = 1;
            }
        }
        geom_t *q = (geom_t *)SMESH_ALLOC((size_t)w.n_elem * sizeof(geom_t));
        mesh_element_quality<idx_t, geom_t>(w.et, w.n_elem, w.elems, w.sdim, w.pts, q);
        uint8_t        *busy   = (uint8_t *)SMESH_CALLOC((size_t)w.n_nodes, sizeof(uint8_t));
        ptrdiff_t       n_busy = w.n_nodes;
        ptrdiff_t       nops   = 0;
        auto            reset_busy = [&]() -> int {
            if (w.n_nodes > n_busy) {
                uint8_t *nb = (uint8_t *)SMESH_CALLOC((size_t)w.n_nodes, sizeof(uint8_t));
                if (!nb) {
                    return SMESH_FAILURE;
                }
                SMESH_FREE(busy);
                busy   = nb;
                n_busy = w.n_nodes;
            } else {
                memset(busy, 0, (size_t)n_busy);
            }
            return SMESH_SUCCESS;
        };
        auto refresh_q = [&]() {
            SMESH_FREE(q);
            q = (geom_t *)SMESH_ALLOC((size_t)w.n_elem * sizeof(geom_t));
            if (q) {
                mesh_element_quality<idx_t, geom_t>(w.et, w.n_elem, w.elems, w.sdim, w.pts, q);
            }
        };
        auto fail_pass = [&]() -> int {
            SMESH_FREE(rowptr);
            SMESH_FREE(colidx);
            SMESH_FREE(uid);
            SMESH_FREE(eu);
            SMESH_FREE(ev);
            SMESH_FREE(feat);
            SMESH_FREE(n2eptr);
            SMESH_FREE(elindex);
            SMESH_FREE(q);
            SMESH_FREE(busy);
            SMESH_FREE(w.se0);
            SMESH_FREE(w.se1);
            improve_free_t<idx_t, geom_t>(w.nxe, w.sdim, w.elems, w.pts, w.lock, w.surf, w.x0);
            return SMESH_FAILURE;
        };
        if (w.is_quad) {
            if (allow_split) {
                nops += quad_split<idx_t, count_t, geom_t>(w, rowptr, colidx, uid, n_uedge, q);
            }
        } else {
            if (allow_swap && w.is_tri) {
                nops += tri_flips<idx_t, count_t, geom_t>(
                        w, rowptr, colidx, uid, eu, ev, feat, n_uedge, n2eptr, elindex, busy);
            }
            if (allow_swap && w.is_tet) {
                element_idx_t *adj = nullptr;
                create_element_adj_table<idx_t, count_t, element_idx_t>(
                        w.n_elem, w.n_nodes, w.et, w.elems, &adj);
                nops += tet_23<idx_t, count_t, geom_t>(w, adj, rowptr, colidx, uid, feat, busy);
                SMESH_FREE(adj);
                if (compact_w(w) != SMESH_SUCCESS ||
                    rebuild_improve_graphs<idx_t, count_t, geom_t>(w,
                                                                   &n2eptr,
                                                                   &elindex,
                                                                   &rowptr,
                                                                   &colidx,
                                                                   &n_uedge,
                                                                   &uid,
                                                                   &eu,
                                                                   &ev,
                                                                   &feat) != SMESH_SUCCESS) {
                    return fail_pass();
                }
                if (reset_busy() != SMESH_SUCCESS) {
                    return fail_pass();
                }
                nops += tet_32<idx_t, count_t, geom_t>(
                        w, eu, ev, feat, n_uedge, n2eptr, elindex, busy);
            }
            if (compact_w(w) != SMESH_SUCCESS ||
                rebuild_improve_graphs<idx_t, count_t, geom_t>(w,
                                                               &n2eptr,
                                                               &elindex,
                                                               &rowptr,
                                                               &colidx,
                                                               &n_uedge,
                                                               &uid,
                                                               &eu,
                                                               &ev,
                                                               &feat) != SMESH_SUCCESS) {
                return fail_pass();
            }
            refresh_q();
            if (!q || reset_busy() != SMESH_SUCCESS) {
                return fail_pass();
            }
            if (allow_split) {
                nops += edge_splits<idx_t, count_t, geom_t>(
                        w, rowptr, colidx, uid, eu, ev, feat, n_uedge, n2eptr, elindex, q, busy);
            }
            if (allow_collapse) {
                if (compact_w(w) != SMESH_SUCCESS ||
                    rebuild_improve_graphs<idx_t, count_t, geom_t>(w,
                                                                   &n2eptr,
                                                                   &elindex,
                                                                   &rowptr,
                                                                   &colidx,
                                                                   &n_uedge,
                                                                   &uid,
                                                                   &eu,
                                                                   &ev,
                                                                   &feat) != SMESH_SUCCESS) {
                    return fail_pass();
                }
                if (reset_busy() != SMESH_SUCCESS) {
                    return fail_pass();
                }
                nops += edge_collapses<idx_t, count_t, geom_t>(
                        w, rowptr, colidx, eu, ev, feat, n_uedge, n2eptr, elindex, busy);
            }
        }
        total_ops += nops;
        SMESH_FREE(rowptr);
        SMESH_FREE(colidx);
        SMESH_FREE(uid);
        SMESH_FREE(eu);
        SMESH_FREE(ev);
        SMESH_FREE(feat);
        SMESH_FREE(n2eptr);
        SMESH_FREE(elindex);
        SMESH_FREE(q);
        SMESH_FREE(busy);
        if (nops == 0) {
            break;
        }
    }
    compact_w(w);

    idx_t **eo = (idx_t **)SMESH_ALLOC((size_t)w.nxe * sizeof(idx_t *));
    for (int d = 0; d < w.nxe; ++d) {
        eo[d] = (idx_t *)SMESH_ALLOC((size_t)w.n_elem * sizeof(idx_t));
        memcpy(eo[d], w.elems[d], (size_t)w.n_elem * sizeof(idx_t));
    }
    geom_t **po = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    geom_t **xo = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    for (int d = 0; d < sdim; ++d) {
        po[d] = (geom_t *)SMESH_ALLOC((size_t)w.n_nodes * sizeof(geom_t));
        xo[d] = (geom_t *)SMESH_ALLOC((size_t)w.n_nodes * sizeof(geom_t));
        memcpy(po[d], w.pts[d], (size_t)w.n_nodes * sizeof(geom_t));
        memcpy(xo[d], w.x0[d], (size_t)w.n_nodes * sizeof(geom_t));
    }
    uint8_t *lo = (uint8_t *)SMESH_ALLOC((size_t)w.n_nodes);
    uint8_t *so = (uint8_t *)SMESH_ALLOC((size_t)w.n_nodes);
    memcpy(lo, w.lock, (size_t)w.n_nodes);
    memcpy(so, w.surf, (size_t)w.n_nodes);
    *n_elements_out = w.n_elem;
    *elements_out   = eo;
    *n_nodes_out    = w.n_nodes;
    *points_out     = po;
    if (lock_out) {
        *lock_out = lo;
    } else {
        SMESH_FREE(lo);
    }
    if (surface_out) {
        *surface_out = so;
    } else {
        SMESH_FREE(so);
    }
    if (x0_out) {
        *x0_out = xo;
    } else {
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(xo[d]);
        }
        SMESH_FREE(xo);
    }
    if (n_ops_out) {
        *n_ops_out = total_ops;
    }
    SMESH_FREE(w.se0);
    SMESH_FREE(w.se1);
    improve_free_t<idx_t, geom_t>(w.nxe, w.sdim, w.elems, w.pts, w.lock, w.surf, w.x0);
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
