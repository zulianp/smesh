#ifndef SMESH_REFINE_ADAPTIVE_IMPL_HPP
#define SMESH_REFINE_ADAPTIVE_IMPL_HPP

#include "smesh_alloc.hpp"
#include "smesh_refine_adaptive.hpp"

#include <algorithm>
#include <cstring>
#include <stdint.h>

namespace smesh {

static ptrdiff_t adapt_next_cap(ptrdiff_t old_cap, ptrdiff_t need) {
    if (need <= old_cap) {
        return old_cap;
    }
    ptrdiff_t ncap = (old_cap > 0) ? old_cap : (ptrdiff_t)1;
    while (ncap < need) {
        if (ncap > PTRDIFF_MAX / 2) {
            return need;
        }
        ncap *= 2;
    }
    return ncap;
}

template <typename T>
static int adapt_grow(T **ptr, ptrdiff_t old_cap, ptrdiff_t new_cap) {
    if (new_cap <= old_cap && *ptr) {
        return SMESH_SUCCESS;
    }
    T *np = (T *)SMESH_REALLOC(*ptr, (size_t)new_cap * sizeof(T));
    if (!np) {
        return SMESH_FAILURE;
    }
    *ptr = np;
    return SMESH_SUCCESS;
}

template <typename idx_t>
static void adapt_order_edge(idx_t *lo, idx_t *hi) {
    if (*lo > *hi) {
        const idx_t tmp = *lo;
        *lo             = *hi;
        *hi             = tmp;
    }
}

template <typename idx_t>
static ptrdiff_t adapt_find_mid(const idx_t *lo, const idx_t *hi, ptrdiff_t n, idx_t a, idx_t b) {
    adapt_order_edge(&a, &b);
    ptrdiff_t L = 0;
    ptrdiff_t R = n;
    while (L < R) {
        const ptrdiff_t m = L + (R - L) / 2;
        if (lo[m] < a || (lo[m] == a && hi[m] < b)) {
            L = m + 1;
        } else {
            R = m;
        }
    }
    if (L < n && lo[L] == a && hi[L] == b) {
        return L;
    }
    return (ptrdiff_t)-1;
}

template <typename idx_t>
struct AdaptEdgeRec {
    idx_t lo;
    idx_t hi;
    idx_t leaf;
};

template <typename idx_t>
struct AdaptFaceRec {
    idx_t n0;
    idx_t n1;
    idx_t n2;
    idx_t leaf;
};

template <typename idx_t>
static void adapt_order_face(idx_t *a, idx_t *b, idx_t *c) {
    if (*a > *b) {
        const idx_t t = *a;
        *a            = *b;
        *b            = t;
    }
    if (*b > *c) {
        const idx_t t = *b;
        *b            = *c;
        *c            = t;
    }
    if (*a > *b) {
        const idx_t t = *a;
        *a            = *b;
        *b            = t;
    }
}

template <typename idx_t>
struct AdaptPair {
    idx_t lo;
    idx_t hi;
};

template <typename idx_t>
static ptrdiff_t adapt_lower_edge(const AdaptEdgeRec<idx_t> *e, ptrdiff_t n, idx_t lo, idx_t hi) {
    ptrdiff_t L = 0;
    ptrdiff_t R = n;
    while (L < R) {
        const ptrdiff_t m = L + (R - L) / 2;
        if (e[m].lo < lo || (e[m].lo == lo && e[m].hi < hi)) {
            L = m + 1;
        } else {
            R = m;
        }
    }
    return L;
}

/// Inward face corners, so the tet centroid is on the positive side.
static constexpr int tet4_inward_face[4][3] = {{0, 3, 1}, {1, 3, 2}, {0, 2, 3}, {0, 1, 2}};
/// Tet-edge index of each side of `tet4_inward_face`, in corner order.
static constexpr int tet4_inward_face_edge[4][3] = {{3, 4, 0}, {4, 5, 1}, {2, 5, 3}, {0, 1, 2}};
/// The two local faces incident to each tet edge. The lower index owns the edgeset entry.
static constexpr int tet4_edge_faces[6][2] = {{0, 3}, {1, 3}, {2, 3}, {0, 2}, {0, 1}, {1, 2}};

static int adapt_slot_parent_face(u8 sl, int face) {
    if (sl >= TET4_ADAPT_SLOT_FAN) {
        if (face != 3) {
            return -1;
        }
        return (int)sl - TET4_ADAPT_SLOT_FAN;
    }
    if (sl < TET4_ADAPT_SLOT_GREEN_EDGE) {
        return tet4_child_face_parent[sl][face];
    }
    if (sl < TET4_ADAPT_SLOT_GREEN_FACE) {
        const int code = (int)sl - TET4_ADAPT_SLOT_GREEN_EDGE;
        return tet4_green_edge_face[code >> 1][code & 1][face];
    }
    const int code = (int)sl - TET4_ADAPT_SLOT_GREEN_FACE;
    return tet4_green_face_face[code >> 2][code & 3][face];
}

static int adapt_slot_parent_edge(u8 sl, int lei) {
    if (sl < TET4_ADAPT_SLOT_GREEN_EDGE) {
        return tet4_child_edge_parent[sl][lei];
    }
    if (sl < TET4_ADAPT_SLOT_GREEN_FACE) {
        const int code = (int)sl - TET4_ADAPT_SLOT_GREEN_EDGE;
        return tet4_green_edge_edge[code >> 1][code & 1][lei];
    }
    const int code = (int)sl - TET4_ADAPT_SLOT_GREEN_FACE;
    return tet4_green_face_edge[code >> 2][code & 3][lei];
}

template <typename idx_t>
static int adapt_map_face(const idx_t *parent, const u8 *slot, const idx_t *root, idx_t leaf, int face, idx_t *root_e,
                          int *root_face) {
    const idx_t none = (idx_t)-1;
    idx_t       cur  = leaf;
    int         f    = face;
    for (int step = 0; parent[cur] != none; ++step) {
        if (step > 64) {
            return 0;
        }
        const int pf = adapt_slot_parent_face(slot[cur], f);
        if (pf < 0) {
            return 0;
        }
        f   = pf;
        cur = parent[cur];
    }
    *root_e    = root[cur];
    *root_face = f;
    return 1;
}

template <typename idx_t>
static int adapt_point_on_edge(idx_t n, idx_t a, idx_t b, const idx_t *as_lo, const idx_t *as_hi, ptrdiff_t n_pts) {
    if (n == a || n == b) {
        return 1;
    }
    const idx_t none = (idx_t)-1;
    if ((ptrdiff_t)n < 0 || (ptrdiff_t)n >= n_pts || as_lo[n] == none) {
        return 0;
    }
    idx_t lo = as_lo[n];
    idx_t hi = as_hi[n];
    adapt_order_edge(&lo, &hi);
    adapt_order_edge(&a, &b);
    return lo == a && hi == b;
}

template <typename idx_t>
static int adapt_parent_edge_of(const idx_t *const *vx, idx_t tet, idx_t p, idx_t q, const idx_t *as_lo,
                                const idx_t *as_hi, ptrdiff_t n_pts) {
    for (int e = 0; e < 6; ++e) {
        const idx_t a = vx[tet4_refine_edges[e][0]][tet];
        const idx_t b = vx[tet4_refine_edges[e][1]][tet];
        if (adapt_point_on_edge(p, a, b, as_lo, as_hi, n_pts) && adapt_point_on_edge(q, a, b, as_lo, as_hi, n_pts)) {
            return e;
        }
    }
    return -1;
}

/// Triangles of one face, oriented so the tet centroid is on the positive side.
/// Missing edge midpoints are `none`. At most 4 triangles.
template <typename idx_t>
static int adapt_face_tris(idx_t c0, idx_t c1, idx_t c2, idx_t m01, idx_t m12, idx_t m20, idx_t (*out)[3]) {
    const idx_t none = (idx_t)-1;
    const int   h01  = m01 != none;
    const int   h12  = m12 != none;
    const int   h20  = m20 != none;
    int         w    = 0;
    auto        push = [&](idx_t a, idx_t b, idx_t c) {
        out[w][0] = a;
        out[w][1] = b;
        out[w][2] = c;
        ++w;
    };
    if (h01 + h12 + h20 == 0) {
        push(c0, c1, c2);
    } else if (h01 + h12 + h20 == 3) {
        push(c0, m01, m20);
        push(m01, c1, m12);
        push(m20, m12, c2);
        push(m01, m12, m20);
    } else if (h01 && !h12 && !h20) {
        push(c0, m01, c2);
        push(m01, c1, c2);
    } else if (!h01 && h12 && !h20) {
        push(c0, c1, m12);
        push(c0, m12, c2);
    } else if (!h01 && !h12 && h20) {
        push(c0, c1, m20);
        push(c1, c2, m20);
    } else if (h01 && h12) {
        push(m01, c1, m12);
        push(c0, m01, m12);
        push(c0, m12, c2);
    } else if (h12 && h20) {
        push(m12, c2, m20);
        push(c0, c1, m12);
        push(c0, m12, m20);
    } else {
        push(c0, m01, m20);
        push(m01, c1, c2);
        push(m01, c2, m20);
    }
    return w;
}

template <typename idx_t>
static int adapt_map_edge(const idx_t *parent, const u8 *slot, const idx_t *root, const idx_t *const *vx,
                          const idx_t *as_lo, const idx_t *as_hi, ptrdiff_t n_pts, idx_t leaf, int lei, idx_t *root_e,
                          int *root_lei) {
    const idx_t none = (idx_t)-1;
    idx_t       cur  = leaf;
    int         e    = lei;
    for (int step = 0; parent[cur] != none; ++step) {
        if (step > 64) {
            return 0;
        }
        if (slot[cur] >= TET4_ADAPT_SLOT_FAN) {
            if (e >= 3) {
                return 0;
            }
            const idx_t par = parent[cur];
            const idx_t p   = vx[tet4_refine_edges[e][0]][cur];
            const idx_t q   = vx[tet4_refine_edges[e][1]][cur];
            const int   pe  = adapt_parent_edge_of(vx, par, p, q, as_lo, as_hi, n_pts);
            if (pe < 0) {
                return 0;
            }
            const int face  = (int)slot[cur] - TET4_ADAPT_SLOT_FAN;
            const int low   = tet4_edge_faces[pe][0] < tet4_edge_faces[pe][1] ? tet4_edge_faces[pe][0]
                                                                              : tet4_edge_faces[pe][1];
            if (face != low) {
                return 0;
            }
            e   = pe;
            cur = par;
            continue;
        }
        const int pe = adapt_slot_parent_edge(slot[cur], e);
        if (pe < 0) {
            return 0;
        }
        e   = pe;
        cur = parent[cur];
    }
    *root_e   = root[cur];
    *root_lei = e;
    return 1;
}

template <typename idx_t>
static int adapt_leaf_mask(const idx_t *const *vx, const idx_t *mlo, const idx_t *mhi, ptrdiff_t n_mids, idx_t leaf) {
    int mask = 0;
    for (int e = 0; e < 6; ++e) {
        const idx_t a = vx[tet4_refine_edges[e][0]][leaf];
        const idx_t b = vx[tet4_refine_edges[e][1]][leaf];
        if (adapt_find_mid(mlo, mhi, n_mids, a, b) >= 0) {
            mask |= (1 << e);
        }
    }
    return mask;
}

template <typename idx_t, typename FaceFn, typename EdgeFn>
static int adapt_visit_boundary(const idx_t *parent, const u8 *slot, const idx_t *root, const idx_t *const *vx,
                                const idx_t *as_lo, const idx_t *as_hi, ptrdiff_t n_pts, idx_t leaf, int kind, int which,
                                FaceFn face_fn, EdgeFn edge_fn) {
    int nc = 1;
    if (kind == TET4_ADAPT_EDGE) {
        nc = 2;
    } else if (kind == TET4_ADAPT_FACE) {
        nc = 4;
    } else if (kind != TET4_ADAPT_KEEP) {
        return SMESH_FAILURE;
    }

    for (int g = 0; g < nc; ++g) {
        for (int lf = 0; lf < 4; ++lf) {
            int leaf_face = lf;
            if (kind == TET4_ADAPT_EDGE) {
                leaf_face = tet4_green_edge_face[which][g][lf];
            } else if (kind == TET4_ADAPT_FACE) {
                leaf_face = tet4_green_face_face[which][g][lf];
            }
            if (leaf_face < 0) {
                continue;
            }
            idx_t root_e = (idx_t)-1;
            int   root_f = -1;
            if (!adapt_map_face(parent, slot, root, leaf, leaf_face, &root_e, &root_f)) {
                continue;
            }
            face_fn(root_e, root_f, lf);
        }
        for (int lei = 0; lei < 6; ++lei) {
            int leaf_edge = lei;
            if (kind == TET4_ADAPT_EDGE) {
                leaf_edge = tet4_green_edge_edge[which][g][lei];
            } else if (kind == TET4_ADAPT_FACE) {
                leaf_edge = tet4_green_face_edge[which][g][lei];
            }
            if (leaf_edge < 0) {
                continue;
            }
            idx_t root_e = (idx_t)-1;
            int   root_l = -1;
            if (!adapt_map_edge(parent, slot, root, vx, as_lo, as_hi, n_pts, leaf, leaf_edge, &root_e, &root_l)) {
                continue;
            }
            edge_fn(root_e, root_l, lei);
        }
    }
    return SMESH_SUCCESS;
}

template <typename idx_t>
static int adapt_idx_room(ptrdiff_t n) {
    if (sizeof(idx_t) >= sizeof(ptrdiff_t)) {
        return n >= 0;
    }
    const idx_t idx_max = (idx_t)(((uintmax_t)1 << (8 * sizeof(idx_t) - 1)) - 1);
    return n >= 0 && n <= (ptrdiff_t)idx_max;
}

template <typename count_t>
static int adapt_count_room(ptrdiff_t n) {
    if (n < 0) {
        return 0;
    }
    if (sizeof(count_t) >= sizeof(ptrdiff_t)) {
        return 1;
    }
    const count_t count_max = (count_t)(((uintmax_t)1 << (8 * sizeof(count_t) - 1)) - 1);
    return n <= (ptrdiff_t)count_max;
}

template <typename idx_t, typename count_t, typename geom_t>
int mesh_refine_adaptive(
        const ptrdiff_t n_elements, const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
        const int spatial_dim, const ptrdiff_t n_nodes, const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
        const idx_t *const SMESH_RESTRICT marked, const ptrdiff_t n_marked, const int level, idx_t ***out_elements,
        ptrdiff_t *n_out_elements, geom_t ***out_points, ptrdiff_t *n_out_nodes, count_t **face_ptr, idx_t **face_elem,
        idx_t **face_lfi, count_t **edge_ptr, idx_t **edge_elem, idx_t **edge_lei, idx_t **mid_lo, idx_t **mid_hi,
        idx_t **mid_id, ptrdiff_t *n_mids_out, idx_t **out_root, const int *mark_level) {
    if (out_elements) {
        *out_elements = nullptr;
    }
    if (out_points) {
        *out_points = nullptr;
    }
    if (n_out_elements) {
        *n_out_elements = 0;
    }
    if (n_out_nodes) {
        *n_out_nodes = 0;
    }
    if (face_ptr) {
        *face_ptr = nullptr;
    }
    if (face_elem) {
        *face_elem = nullptr;
    }
    if (face_lfi) {
        *face_lfi = nullptr;
    }
    if (edge_ptr) {
        *edge_ptr = nullptr;
    }
    if (edge_elem) {
        *edge_elem = nullptr;
    }
    if (edge_lei) {
        *edge_lei = nullptr;
    }
    if (mid_lo) {
        *mid_lo = nullptr;
    }
    if (mid_hi) {
        *mid_hi = nullptr;
    }
    if (mid_id) {
        *mid_id = nullptr;
    }
    if (n_mids_out) {
        *n_mids_out = 0;
    }
    if (out_root) {
        *out_root = nullptr;
    }

    if (!out_elements || !n_out_elements || !out_points || !n_out_nodes || !face_ptr || !face_elem || !face_lfi ||
        !edge_ptr || !edge_elem || !edge_lei || !mid_lo || !mid_hi || !mid_id || !n_mids_out || level < 1 ||
        n_elements < 1 || n_nodes < 4 || spatial_dim < 1 || spatial_dim > 3 || !elements || !points ||
        (n_marked > 0 && !marked)) {
        fprintf(stderr, "refine_adaptive: invalid kernel arguments\n");
        return SMESH_FAILURE;
    }

    const idx_t none      = (idx_t)-1;
    int         depth_cap = level;

    idx_t *vx[4]     = {nullptr, nullptr, nullptr, nullptr};
    int   *tet_level = nullptr;
    int   *tet_target = nullptr;
    idx_t *child0    = nullptr;
    idx_t *parent    = nullptr;
    u8    *slot      = nullptr;
    idx_t *root      = nullptr;
    ptrdiff_t n_forest = 0;
    ptrdiff_t cap_forest = 0;

    geom_t   *pt[3]  = {nullptr, nullptr, nullptr};
    ptrdiff_t n_pts  = 0;
    ptrdiff_t cap_pts = 0;
    idx_t    *as_lo  = nullptr;
    idx_t    *as_hi  = nullptr;

    idx_t    *mlo    = nullptr;
    idx_t    *mhi    = nullptr;
    idx_t    *mid    = nullptr;
    ptrdiff_t n_mids = 0;

    idx_t    *clo      = nullptr;
    idx_t    *chi      = nullptr;
    idx_t    *cid      = nullptr;
    ptrdiff_t n_created = 0;
    ptrdiff_t cap_created = 0;

    idx_t  **out_e   = nullptr;
    geom_t **out_p   = nullptr;
    count_t *f_ptr   = nullptr;
    idx_t   *f_elem  = nullptr;
    idx_t   *f_lfi   = nullptr;
    count_t *e_ptr   = nullptr;
    idx_t   *e_elem  = nullptr;
    idx_t   *e_lei   = nullptr;
    idx_t   *out_root_e = nullptr;

    cap_forest = adapt_next_cap(0, n_elements);
    if (adapt_grow(&tet_level, (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&tet_target, (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&child0, (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&parent, (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&slot, (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&root, (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&vx[0], (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&vx[1], (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&vx[2], (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS ||
        adapt_grow(&vx[3], (ptrdiff_t)0, cap_forest) != SMESH_SUCCESS) {
        fprintf(stderr, "refine_adaptive: forest allocation failed\n");
        goto adapt_fail;
    }
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
        for (int k = 0; k < 4; ++k) {
            vx[k][e] = elements[k][e];
        }
        tet_level[e]  = 0;
        tet_target[e] = 0;
        child0[e]     = none;
        parent[e]     = none;
        slot[e]       = 0;
        root[e]       = (idx_t)e;
    }
    n_forest = n_elements;

    for (ptrdiff_t i = 0; i < n_marked; ++i) {
        const idx_t m = marked[i];
        if (m < 0 || (ptrdiff_t)m >= n_elements) {
            fprintf(stderr, "refine_adaptive: marked element out of range\n");
            goto adapt_fail;
        }
        int depth = level;
        if (mark_level) {
            depth = mark_level[i];
            if (depth < 0) {
                fprintf(stderr, "refine_adaptive: level must be non-negative\n");
                goto adapt_fail;
            }
            if (depth > depth_cap) {
                depth_cap = depth;
            }
        }
        if (tet_target[m] < depth) {
            tet_target[m] = depth;
        }
    }

    // Target L forces every tet that shares a face up to at least L-1.
    // Each pass extends that floor by one face, so a depth-L mark grades a
    // halo L-1 elements wide and does not promote the neighbor to L.
    if (depth_cap > 1) {
        const ptrdiff_t      n_fr = n_elements * 4;
        AdaptFaceRec<idx_t> *fr =
                (AdaptFaceRec<idx_t> *)SMESH_ALLOC((size_t)(n_fr > 0 ? n_fr : 1) * sizeof(AdaptFaceRec<idx_t>));
        if (!fr) {
            fprintf(stderr, "refine_adaptive: face allocation failed\n");
            goto adapt_fail;
        }
        ptrdiff_t nfr = 0;
        for (ptrdiff_t e = 0; e < n_elements; ++e) {
            for (int f = 0; f < 4; ++f) {
                idx_t n0 = vx[tet4_inward_face[f][0]][e];
                idx_t n1 = vx[tet4_inward_face[f][1]][e];
                idx_t n2 = vx[tet4_inward_face[f][2]][e];
                adapt_order_face(&n0, &n1, &n2);
                fr[nfr].n0   = n0;
                fr[nfr].n1   = n1;
                fr[nfr].n2   = n2;
                fr[nfr].leaf = (idx_t)e;
                ++nfr;
            }
        }
        std::sort(fr, fr + nfr, [](const AdaptFaceRec<idx_t> &a, const AdaptFaceRec<idx_t> &b) {
            if (a.n0 != b.n0) {
                return a.n0 < b.n0;
            }
            if (a.n1 != b.n1) {
                return a.n1 < b.n1;
            }
            if (a.n2 != b.n2) {
                return a.n2 < b.n2;
            }
            return a.leaf < b.leaf;
        });
        for (int hop = 1; hop < depth_cap; ++hop) {
            int       changed = 0;
            ptrdiff_t k       = 0;
            while (k < nfr) {
                ptrdiff_t k1 = k + 1;
                while (k1 < nfr && fr[k1].n0 == fr[k].n0 && fr[k1].n1 == fr[k].n1 && fr[k1].n2 == fr[k].n2) {
                    ++k1;
                }
                int hi_target = tet_target[fr[k].leaf];
                for (ptrdiff_t t = k + 1; t < k1; ++t) {
                    const int target = tet_target[fr[t].leaf];
                    if (target > hi_target) {
                        hi_target = target;
                    }
                }
                const int floor_target = hi_target - 1;
                if (floor_target > 0) {
                    for (ptrdiff_t t = k; t < k1; ++t) {
                        int *const slot = &tet_target[fr[t].leaf];
                        if (*slot < floor_target) {
                            *slot   = floor_target;
                            changed = 1;
                        }
                    }
                }
                k = k1;
            }
            if (!changed) {
                break;
            }
        }
        SMESH_FREE(fr);
    }

    cap_pts = adapt_next_cap(0, n_nodes);
    if (adapt_grow(&as_lo, (ptrdiff_t)0, cap_pts) != SMESH_SUCCESS ||
        adapt_grow(&as_hi, (ptrdiff_t)0, cap_pts) != SMESH_SUCCESS) {
        fprintf(stderr, "refine_adaptive: point allocation failed\n");
        goto adapt_fail;
    }
    for (int d = 0; d < spatial_dim; ++d) {
        if (adapt_grow(&pt[d], (ptrdiff_t)0, cap_pts) != SMESH_SUCCESS) {
            fprintf(stderr, "refine_adaptive: point allocation failed\n");
            goto adapt_fail;
        }
        memcpy(pt[d], points[d], (size_t)n_nodes * sizeof(geom_t));
    }
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        as_lo[i] = none;
        as_hi[i] = none;
    }
    n_pts = n_nodes;

    {
        const ptrdiff_t pass_limit = n_elements * (ptrdiff_t)(depth_cap + 3) + 8;
        for (ptrdiff_t pass = 0; pass < pass_limit; ++pass) {
            ptrdiff_t n_leaves = 0;
            for (ptrdiff_t i = 0; i < n_forest; ++i) {
                if (child0[i] == none) {
                    ++n_leaves;
                }
            }
            idx_t *leaves = (idx_t *)SMESH_ALLOC((size_t)n_leaves * sizeof(idx_t));
            if (!leaves) {
                fprintf(stderr, "refine_adaptive: leaf allocation failed\n");
                goto adapt_fail;
            }
            ptrdiff_t nl = 0;
            for (ptrdiff_t i = 0; i < n_forest; ++i) {
                if (child0[i] == none) {
                    leaves[nl++] = (idx_t)i;
                }
            }

            const ptrdiff_t     n_er = n_leaves * 6;
            AdaptEdgeRec<idx_t> *er  = (AdaptEdgeRec<idx_t> *)SMESH_ALLOC((size_t)(n_er > 0 ? n_er : 1) * sizeof(AdaptEdgeRec<idx_t>));
            if (!er) {
                SMESH_FREE(leaves);
                fprintf(stderr, "refine_adaptive: edge allocation failed\n");
                goto adapt_fail;
            }
            ptrdiff_t ner = 0;
            for (ptrdiff_t i = 0; i < n_leaves; ++i) {
                const idx_t f = leaves[i];
                for (int e = 0; e < 6; ++e) {
                    idx_t a = vx[tet4_refine_edges[e][0]][f];
                    idx_t b = vx[tet4_refine_edges[e][1]][f];
                    adapt_order_edge(&a, &b);
                    er[ner].lo   = a;
                    er[ner].hi   = b;
                    er[ner].leaf = f;
                    ++ner;
                }
            }
            std::sort(er, er + ner, [](const AdaptEdgeRec<idx_t> &a, const AdaptEdgeRec<idx_t> &b) {
                if (a.lo != b.lo) {
                    return a.lo < b.lo;
                }
                if (a.hi != b.hi) {
                    return a.hi < b.hi;
                }
                return a.leaf < b.leaf;
            });

            u8 *blocked = (u8 *)SMESH_CALLOC((size_t)(n_leaves > 0 ? n_leaves : 1), sizeof(u8));
            if (!blocked) {
                SMESH_FREE(er);
                SMESH_FREE(leaves);
                fprintf(stderr, "refine_adaptive: allocation failed\n");
                goto adapt_fail;
            }

            // A coarser neighbor blocks this red split until the shared edge is a
            // real edge of a transition child. Do not raise that neighbor to a full
            // 8-split: on the 24-tet hex that bisects every spoke and the closure
            // walks the mesh. Green splits below match the hanging node only.
            for (ptrdiff_t i = 0; i < n_leaves; ++i) {
                const idx_t f = leaves[i];
                if (tet_target[f] <= tet_level[f]) {
                    continue;
                }
                const int my_level = tet_level[f];
                for (int e = 0; e < 6; ++e) {
                    idx_t a = vx[tet4_refine_edges[e][0]][f];
                    idx_t b = vx[tet4_refine_edges[e][1]][f];
                    adapt_order_edge(&a, &b);
                    const ptrdiff_t at = adapt_lower_edge(er, ner, a, b);
                    for (ptrdiff_t k = at; k < ner && er[k].lo == a && er[k].hi == b; ++k) {
                        const idx_t other = er[k].leaf;
                        if (other == f) {
                            continue;
                        }
                        if (tet_level[other] < my_level) {
                            blocked[i] = 1;
                        }
                    }
                    const idx_t ends[2] = {vx[tet4_refine_edges[e][0]][f], vx[tet4_refine_edges[e][1]][f]};
                    for (int h = 0; h < 2; ++h) {
                        const idx_t node = ends[h];
                        if ((ptrdiff_t)node >= n_pts || as_lo[node] == none) {
                            continue;
                        }
                        idx_t       plo = as_lo[node];
                        idx_t       phi = as_hi[node];
                        adapt_order_edge(&plo, &phi);
                        const ptrdiff_t pt_at = adapt_lower_edge(er, ner, plo, phi);
                        for (ptrdiff_t k = pt_at; k < ner && er[k].lo == plo && er[k].hi == phi; ++k) {
                            const idx_t other = er[k].leaf;
                            if (tet_level[other] < my_level) {
                                blocked[i] = 1;
                            }
                        }
                    }
                }
            }

            ptrdiff_t n_sel = 0;
            for (ptrdiff_t i = 0; i < n_leaves; ++i) {
                const idx_t f = leaves[i];
                if (!blocked[i] && tet_target[f] > tet_level[f]) {
                    ++n_sel;
                }
            }
            idx_t *sel = nullptr;
            if (n_sel > 0) {
                sel = (idx_t *)SMESH_ALLOC((size_t)n_sel * sizeof(idx_t));
                if (!sel) {
                    SMESH_FREE(blocked);
                    SMESH_FREE(er);
                    SMESH_FREE(leaves);
                    fprintf(stderr, "refine_adaptive: allocation failed\n");
                    goto adapt_fail;
                }
                ptrdiff_t ns = 0;
                for (ptrdiff_t i = 0; i < n_leaves; ++i) {
                    const idx_t f = leaves[i];
                    if (!blocked[i] && tet_target[f] > tet_level[f]) {
                        sel[ns++] = f;
                    }
                }
            }

            int refined = 0;
            if (n_sel > 0) {
                const ptrdiff_t     n_pairs_cap = n_sel * 6;
                AdaptPair<idx_t>   *pairs =
                        (AdaptPair<idx_t> *)SMESH_ALLOC((size_t)n_pairs_cap * sizeof(AdaptPair<idx_t>));
                if (!pairs) {
                    SMESH_FREE(sel);
                    SMESH_FREE(blocked);
                    SMESH_FREE(er);
                    SMESH_FREE(leaves);
                    fprintf(stderr, "refine_adaptive: allocation failed\n");
                    goto adapt_fail;
                }
                ptrdiff_t npairs = 0;
                for (ptrdiff_t s = 0; s < n_sel; ++s) {
                    const idx_t f = sel[s];
                    for (int e = 0; e < 6; ++e) {
                        idx_t a = vx[tet4_refine_edges[e][0]][f];
                        idx_t b = vx[tet4_refine_edges[e][1]][f];
                        adapt_order_edge(&a, &b);
                        pairs[npairs].lo = a;
                        pairs[npairs].hi = b;
                        ++npairs;
                    }
                }
                std::sort(pairs, pairs + npairs, [](const AdaptPair<idx_t> &a, const AdaptPair<idx_t> &b) {
                    if (a.lo != b.lo) {
                        return a.lo < b.lo;
                    }
                    return a.hi < b.hi;
                });
                ptrdiff_t nu = 0;
                for (ptrdiff_t p = 0; p < npairs; ++p) {
                    if (nu > 0 && pairs[nu - 1].lo == pairs[p].lo && pairs[nu - 1].hi == pairs[p].hi) {
                        continue;
                    }
                    pairs[nu++] = pairs[p];
                }

                ptrdiff_t n_new = 0;
                for (ptrdiff_t p = 0; p < nu; ++p) {
                    if (adapt_find_mid(mlo, mhi, n_mids, pairs[p].lo, pairs[p].hi) < 0) {
                        ++n_new;
                    }
                }

                if (n_pts > PTRDIFF_MAX - n_new || !adapt_idx_room<idx_t>(n_pts + n_new)) {
                    SMESH_FREE(pairs);
                    SMESH_FREE(sel);
                    SMESH_FREE(blocked);
                    SMESH_FREE(er);
                    SMESH_FREE(leaves);
                    fprintf(stderr, "refine_adaptive: node index overflow\n");
                    goto adapt_fail;
                }

                const ptrdiff_t old_pts = cap_pts;
                const ptrdiff_t old_created = cap_created;
                cap_pts                     = adapt_next_cap(old_pts, n_pts + n_new);
                cap_created                 = adapt_next_cap(old_created, n_created + n_new);
                if (adapt_grow(&as_lo, old_pts, cap_pts) != SMESH_SUCCESS ||
                    adapt_grow(&as_hi, old_pts, cap_pts) != SMESH_SUCCESS ||
                    adapt_grow(&clo, old_created, cap_created) != SMESH_SUCCESS ||
                    adapt_grow(&chi, old_created, cap_created) != SMESH_SUCCESS ||
                    adapt_grow(&cid, old_created, cap_created) != SMESH_SUCCESS ||
                    adapt_grow(&pt[0], old_pts, cap_pts) != SMESH_SUCCESS ||
                    (spatial_dim > 1 && adapt_grow(&pt[1], old_pts, cap_pts) != SMESH_SUCCESS) ||
                    (spatial_dim > 2 && adapt_grow(&pt[2], old_pts, cap_pts) != SMESH_SUCCESS)) {
                    SMESH_FREE(pairs);
                    SMESH_FREE(sel);
                    SMESH_FREE(blocked);
                    SMESH_FREE(er);
                    SMESH_FREE(leaves);
                    fprintf(stderr, "refine_adaptive: midpoint allocation failed\n");
                    goto adapt_fail;
                }
                for (ptrdiff_t i = n_pts; i < cap_pts; ++i) {
                    as_lo[i] = none;
                    as_hi[i] = none;
                }

                idx_t *new_lo = nullptr;
                idx_t *new_hi = nullptr;
                idx_t *new_id = nullptr;
                if (n_new > 0) {
                    new_lo = (idx_t *)SMESH_ALLOC((size_t)n_new * sizeof(idx_t));
                    new_hi = (idx_t *)SMESH_ALLOC((size_t)n_new * sizeof(idx_t));
                    new_id = (idx_t *)SMESH_ALLOC((size_t)n_new * sizeof(idx_t));
                    if (!new_lo || !new_hi || !new_id) {
                        SMESH_FREE(new_lo);
                        SMESH_FREE(new_hi);
                        SMESH_FREE(new_id);
                        SMESH_FREE(pairs);
                        SMESH_FREE(sel);
                        SMESH_FREE(blocked);
                        SMESH_FREE(er);
                        SMESH_FREE(leaves);
                        fprintf(stderr, "refine_adaptive: midpoint allocation failed\n");
                        goto adapt_fail;
                    }
                }

                ptrdiff_t n_fill = 0;
                idx_t     next   = (idx_t)n_pts;
                for (ptrdiff_t p = 0; p < nu; ++p) {
                    if (adapt_find_mid(mlo, mhi, n_mids, pairs[p].lo, pairs[p].hi) >= 0) {
                        continue;
                    }
                    const idx_t lo = pairs[p].lo;
                    const idx_t hi = pairs[p].hi;
                    new_lo[n_fill] = lo;
                    new_hi[n_fill] = hi;
                    new_id[n_fill] = next;
                    for (int d = 0; d < spatial_dim; ++d) {
                        pt[d][next] = (pt[d][lo] + pt[d][hi]) / (geom_t)2;
                    }
                    as_lo[next]       = lo;
                    as_hi[next]       = hi;
                    clo[n_created]    = lo;
                    chi[n_created]    = hi;
                    cid[n_created]    = next;
                    ++n_created;
                    ++n_fill;
                    next = (idx_t)((ptrdiff_t)next + 1);
                }
                n_pts = (ptrdiff_t)next;
                SMESH_FREE(pairs);

                if (n_new > 0) {
                    const ptrdiff_t n_merged = n_mids + n_new;
                    idx_t          *jlo      = (idx_t *)SMESH_ALLOC((size_t)n_merged * sizeof(idx_t));
                    idx_t          *jhi      = (idx_t *)SMESH_ALLOC((size_t)n_merged * sizeof(idx_t));
                    idx_t          *jid      = (idx_t *)SMESH_ALLOC((size_t)n_merged * sizeof(idx_t));
                    if (!jlo || !jhi || !jid) {
                        SMESH_FREE(jlo);
                        SMESH_FREE(jhi);
                        SMESH_FREE(jid);
                        SMESH_FREE(new_lo);
                        SMESH_FREE(new_hi);
                        SMESH_FREE(new_id);
                        SMESH_FREE(sel);
                        SMESH_FREE(blocked);
                        SMESH_FREE(er);
                        SMESH_FREE(leaves);
                        fprintf(stderr, "refine_adaptive: midpoint allocation failed\n");
                        goto adapt_fail;
                    }
                    ptrdiff_t ia = 0;
                    ptrdiff_t ib = 0;
                    ptrdiff_t iw = 0;
                    while (ia < n_mids && ib < n_new) {
                        const int old_first = mlo[ia] < new_lo[ib] || (mlo[ia] == new_lo[ib] && mhi[ia] < new_hi[ib]);
                        if (old_first) {
                            jlo[iw] = mlo[ia];
                            jhi[iw] = mhi[ia];
                            jid[iw] = mid[ia];
                            ++ia;
                        } else {
                            jlo[iw] = new_lo[ib];
                            jhi[iw] = new_hi[ib];
                            jid[iw] = new_id[ib];
                            ++ib;
                        }
                        ++iw;
                    }
                    while (ia < n_mids) {
                        jlo[iw] = mlo[ia];
                        jhi[iw] = mhi[ia];
                        jid[iw] = mid[ia];
                        ++ia;
                        ++iw;
                    }
                    while (ib < n_new) {
                        jlo[iw] = new_lo[ib];
                        jhi[iw] = new_hi[ib];
                        jid[iw] = new_id[ib];
                        ++ib;
                        ++iw;
                    }
                    SMESH_FREE(mlo);
                    SMESH_FREE(mhi);
                    SMESH_FREE(mid);
                    mlo    = jlo;
                    mhi    = jhi;
                    mid    = jid;
                    n_mids = n_merged;
                }
                SMESH_FREE(new_lo);
                SMESH_FREE(new_hi);
                SMESH_FREE(new_id);

                const ptrdiff_t n_add = n_sel * 8;
                if (n_forest > PTRDIFF_MAX - n_add || !adapt_idx_room<idx_t>(n_forest + n_add)) {
                    SMESH_FREE(sel);
                    SMESH_FREE(blocked);
                    SMESH_FREE(er);
                    SMESH_FREE(leaves);
                    fprintf(stderr, "refine_adaptive: element index overflow\n");
                    goto adapt_fail;
                }
                const ptrdiff_t old_forest = cap_forest;
                cap_forest                 = adapt_next_cap(old_forest, n_forest + n_add);
                if (adapt_grow(&vx[0], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[1], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[2], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[3], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&tet_level, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&tet_target, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&child0, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&parent, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&slot, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&root, old_forest, cap_forest) != SMESH_SUCCESS) {
                    SMESH_FREE(sel);
                    SMESH_FREE(blocked);
                    SMESH_FREE(er);
                    SMESH_FREE(leaves);
                    fprintf(stderr, "refine_adaptive: forest allocation failed\n");
                    goto adapt_fail;
                }

                const ptrdiff_t base = n_forest;
                for (ptrdiff_t s = 0; s < n_sel; ++s) {
                    const idx_t f = sel[s];
                    idx_t       macro[10];
                    for (int k = 0; k < 4; ++k) {
                        macro[k] = vx[k][f];
                    }
                    for (int e = 0; e < 6; ++e) {
                        const ptrdiff_t mi = adapt_find_mid(mlo, mhi, n_mids, macro[tet4_refine_edges[e][0]],
                                                            macro[tet4_refine_edges[e][1]]);
                        if (mi < 0) {
                            SMESH_FREE(sel);
                            SMESH_FREE(blocked);
                            SMESH_FREE(er);
                            SMESH_FREE(leaves);
                            fprintf(stderr, "refine_adaptive: missing edge midpoint\n");
                            goto adapt_fail;
                        }
                        macro[4 + e] = mid[mi];
                    }
                    const ptrdiff_t c0 = base + s * 8;
                    child0[f]          = (idx_t)c0;
                    for (int c = 0; c < 8; ++c) {
                        const ptrdiff_t dst = c0 + c;
                        for (int k = 0; k < 4; ++k) {
                            vx[k][dst] = macro[tet4_refine_pattern[c][k]];
                        }
                        tet_level[dst]  = tet_level[f] + 1;
                        tet_target[dst] = tet_target[f];
                        child0[dst]     = none;
                        parent[dst]     = f;
                        slot[dst]       = (u8)c;
                        root[dst]       = root[f];
                    }
                }
                n_forest = base + n_add;
                refined  = 1;
            }

            SMESH_FREE(sel);
            SMESH_FREE(blocked);
            SMESH_FREE(er);
            SMESH_FREE(leaves);

            ptrdiff_t n_leaves2 = 0;
            for (ptrdiff_t i = 0; i < n_forest; ++i) {
                if (child0[i] == none) {
                    ++n_leaves2;
                }
            }
            idx_t *leaves2 = (idx_t *)SMESH_ALLOC((size_t)(n_leaves2 > 0 ? n_leaves2 : 1) * sizeof(idx_t));
            if (!leaves2) {
                fprintf(stderr, "refine_adaptive: leaf allocation failed\n");
                goto adapt_fail;
            }
            ptrdiff_t nl2 = 0;
            for (ptrdiff_t i = 0; i < n_forest; ++i) {
                if (child0[i] == none) {
                    leaves2[nl2++] = (idx_t)i;
                }
            }

            idx_t *gleaf  = (idx_t *)SMESH_ALLOC((size_t)(n_leaves2 > 0 ? n_leaves2 : 1) * sizeof(idx_t));
            u8    *gkind  = (u8 *)SMESH_ALLOC((size_t)(n_leaves2 > 0 ? n_leaves2 : 1) * sizeof(u8));
            u8    *gwhich = (u8 *)SMESH_ALLOC((size_t)(n_leaves2 > 0 ? n_leaves2 : 1) * sizeof(u8));
            idx_t *fleaf  = (idx_t *)SMESH_ALLOC((size_t)(n_leaves2 > 0 ? n_leaves2 : 1) * sizeof(idx_t));
            if (!gleaf || !gkind || !gwhich || !fleaf) {
                SMESH_FREE(gleaf);
                SMESH_FREE(gkind);
                SMESH_FREE(gwhich);
                SMESH_FREE(fleaf);
                SMESH_FREE(leaves2);
                fprintf(stderr, "refine_adaptive: leaf allocation failed\n");
                goto adapt_fail;
            }

            ptrdiff_t n_green = 0;
            ptrdiff_t n_fan   = 0;
            for (ptrdiff_t i = 0; i < n_leaves2; ++i) {
                const idx_t           f    = leaves2[i];
                const int             mask = adapt_leaf_mask<idx_t>(vx, mlo, mhi, n_mids, f);
                const Tet4AdaptAction act  = tet4_adapt_action[mask];
                if (tet_target[f] > tet_level[f]) {
                    continue;
                }
                if (act.kind == TET4_ADAPT_PROMOTE) {
                    fleaf[n_fan++] = f;
                } else if (act.kind == TET4_ADAPT_EDGE || act.kind == TET4_ADAPT_FACE) {
                    gleaf[n_green]  = f;
                    gkind[n_green]  = act.kind;
                    gwhich[n_green] = act.which;
                    ++n_green;
                }
            }
            SMESH_FREE(leaves2);

            int fanned = 0;
            if (n_fan > 0) {
                ptrdiff_t n_add = 0;
                for (ptrdiff_t i = 0; i < n_fan; ++i) {
                    const idx_t f = fleaf[i];
                    idx_t       corner[4];
                    idx_t       mid_of[6];
                    for (int k = 0; k < 4; ++k) {
                        corner[k] = vx[k][f];
                    }
                    for (int e = 0; e < 6; ++e) {
                        const ptrdiff_t mi = adapt_find_mid(mlo, mhi, n_mids, corner[tet4_refine_edges[e][0]],
                                                            corner[tet4_refine_edges[e][1]]);
                        mid_of[e]           = (mi >= 0) ? mid[mi] : none;
                    }
                    for (int face = 0; face < 4; ++face) {
                        const idx_t c0  = corner[tet4_inward_face[face][0]];
                        const idx_t c1  = corner[tet4_inward_face[face][1]];
                        const idx_t c2  = corner[tet4_inward_face[face][2]];
                        const idx_t m01 = mid_of[tet4_inward_face_edge[face][0]];
                        const idx_t m12 = mid_of[tet4_inward_face_edge[face][1]];
                        const idx_t m20 = mid_of[tet4_inward_face_edge[face][2]];
                        idx_t       tmp[4][3];
                        n_add += adapt_face_tris(c0, c1, c2, m01, m12, m20, tmp);
                    }
                }
                if (n_forest > PTRDIFF_MAX - n_add || n_pts > PTRDIFF_MAX - n_fan ||
                    !adapt_idx_room<idx_t>(n_forest + n_add) || !adapt_idx_room<idx_t>(n_pts + n_fan)) {
                    SMESH_FREE(gleaf);
                    SMESH_FREE(gkind);
                    SMESH_FREE(gwhich);
                    SMESH_FREE(fleaf);
                    fprintf(stderr, "refine_adaptive: element index overflow\n");
                    goto adapt_fail;
                }
                const ptrdiff_t old_forest = cap_forest;
                const ptrdiff_t old_pts    = cap_pts;
                cap_forest                 = adapt_next_cap(old_forest, n_forest + n_add);
                cap_pts                    = adapt_next_cap(old_pts, n_pts + n_fan);
                if (adapt_grow(&vx[0], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[1], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[2], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[3], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&tet_level, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&tet_target, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&child0, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&parent, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&slot, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&root, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&as_lo, old_pts, cap_pts) != SMESH_SUCCESS ||
                    adapt_grow(&as_hi, old_pts, cap_pts) != SMESH_SUCCESS ||
                    adapt_grow(&pt[0], old_pts, cap_pts) != SMESH_SUCCESS ||
                    (spatial_dim > 1 && adapt_grow(&pt[1], old_pts, cap_pts) != SMESH_SUCCESS) ||
                    (spatial_dim > 2 && adapt_grow(&pt[2], old_pts, cap_pts) != SMESH_SUCCESS)) {
                    SMESH_FREE(gleaf);
                    SMESH_FREE(gkind);
                    SMESH_FREE(gwhich);
                    SMESH_FREE(fleaf);
                    fprintf(stderr, "refine_adaptive: forest allocation failed\n");
                    goto adapt_fail;
                }
                for (ptrdiff_t i = n_pts; i < cap_pts; ++i) {
                    as_lo[i] = none;
                    as_hi[i] = none;
                }

                ptrdiff_t cursor = n_forest;
                idx_t     next   = (idx_t)n_pts;
                for (ptrdiff_t i = 0; i < n_fan; ++i) {
                    const idx_t f = fleaf[i];
                    idx_t       corner[4];
                    idx_t       mid_of[6];
                    for (int k = 0; k < 4; ++k) {
                        corner[k] = vx[k][f];
                    }
                    for (int e = 0; e < 6; ++e) {
                        const ptrdiff_t mi = adapt_find_mid(mlo, mhi, n_mids, corner[tet4_refine_edges[e][0]],
                                                            corner[tet4_refine_edges[e][1]]);
                        mid_of[e]           = (mi >= 0) ? mid[mi] : none;
                    }
                    const idx_t cid = next;
                    next            = (idx_t)((ptrdiff_t)next + 1);
                    for (int d = 0; d < spatial_dim; ++d) {
                        pt[d][cid] = (pt[d][corner[0]] + pt[d][corner[1]] + pt[d][corner[2]] + pt[d][corner[3]]) /
                                     (geom_t)4;
                    }
                    child0[f] = (idx_t)cursor;
                    for (int face = 0; face < 4; ++face) {
                        const idx_t c0  = corner[tet4_inward_face[face][0]];
                        const idx_t c1  = corner[tet4_inward_face[face][1]];
                        const idx_t c2  = corner[tet4_inward_face[face][2]];
                        const idx_t m01 = mid_of[tet4_inward_face_edge[face][0]];
                        const idx_t m12 = mid_of[tet4_inward_face_edge[face][1]];
                        const idx_t m20 = mid_of[tet4_inward_face_edge[face][2]];
                        idx_t       tris[4][3];
                        const int   nt = adapt_face_tris(c0, c1, c2, m01, m12, m20, tris);
                        for (int t = 0; t < nt; ++t) {
                            vx[0][cursor]      = tris[t][0];
                            vx[1][cursor]      = tris[t][1];
                            vx[2][cursor]      = tris[t][2];
                            vx[3][cursor]      = cid;
                            tet_level[cursor]  = tet_level[f] + 1;
                            tet_target[cursor] = tet_level[cursor];
                            child0[cursor]     = none;
                            parent[cursor]     = f;
                            slot[cursor]       = (u8)(TET4_ADAPT_SLOT_FAN + face);
                            root[cursor]       = root[f];
                            ++cursor;
                        }
                    }
                }
                n_pts    = (ptrdiff_t)next;
                n_forest = cursor;
                fanned   = 1;
            }
            SMESH_FREE(fleaf);

            int greened = 0;
            if (n_green > 0) {
                ptrdiff_t n_add = 0;
                for (ptrdiff_t i = 0; i < n_green; ++i) {
                    n_add += (gkind[i] == TET4_ADAPT_EDGE) ? (ptrdiff_t)2 : (ptrdiff_t)4;
                }
                if (n_forest > PTRDIFF_MAX - n_add || !adapt_idx_room<idx_t>(n_forest + n_add)) {
                    SMESH_FREE(gleaf);
                    SMESH_FREE(gkind);
                    SMESH_FREE(gwhich);
                    fprintf(stderr, "refine_adaptive: element index overflow\n");
                    goto adapt_fail;
                }
                const ptrdiff_t old_forest = cap_forest;
                cap_forest                 = adapt_next_cap(old_forest, n_forest + n_add);
                if (adapt_grow(&vx[0], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[1], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[2], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&vx[3], old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&tet_level, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&tet_target, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&child0, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&parent, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&slot, old_forest, cap_forest) != SMESH_SUCCESS ||
                    adapt_grow(&root, old_forest, cap_forest) != SMESH_SUCCESS) {
                    SMESH_FREE(gleaf);
                    SMESH_FREE(gkind);
                    SMESH_FREE(gwhich);
                    fprintf(stderr, "refine_adaptive: forest allocation failed\n");
                    goto adapt_fail;
                }

                ptrdiff_t cursor = n_forest;
                for (ptrdiff_t i = 0; i < n_green; ++i) {
                    const idx_t f  = gleaf[i];
                    const int   nk = gkind[i];
                    const int   wh = gwhich[i];
                    const int   nc = (nk == TET4_ADAPT_EDGE) ? 2 : 4;
                    idx_t       macro[10];
                    for (int k = 0; k < 4; ++k) {
                        macro[k] = vx[k][f];
                    }
                    for (int e = 0; e < 6; ++e) {
                        macro[4 + e] = none;
                        const ptrdiff_t mi = adapt_find_mid(mlo, mhi, n_mids, macro[tet4_refine_edges[e][0]],
                                                            macro[tet4_refine_edges[e][1]]);
                        if (mi >= 0) {
                            macro[4 + e] = mid[mi];
                        }
                    }
                    child0[f] = (idx_t)cursor;
                    for (int g = 0; g < nc; ++g) {
                        const int *pat = (nk == TET4_ADAPT_EDGE) ? tet4_green_edge[wh][g] : tet4_green_face[wh][g];
                        const u8   sl  = (nk == TET4_ADAPT_EDGE)
                                               ? (u8)(TET4_ADAPT_SLOT_GREEN_EDGE + wh * 2 + g)
                                               : (u8)(TET4_ADAPT_SLOT_GREEN_FACE + wh * 4 + g);
                        for (int k = 0; k < 4; ++k) {
                            const int src = pat[k];
                            if (src < 0 || src > 9 || macro[src] == none) {
                                SMESH_FREE(gleaf);
                                SMESH_FREE(gkind);
                                SMESH_FREE(gwhich);
                                fprintf(stderr, "refine_adaptive: green pattern referenced a missing node\n");
                                goto adapt_fail;
                            }
                            vx[k][cursor] = macro[src];
                        }
                        tet_level[cursor]  = tet_level[f] + 1;
                        tet_target[cursor] = tet_level[cursor];
                        child0[cursor]     = none;
                        parent[cursor]     = f;
                        slot[cursor]       = sl;
                        root[cursor]       = root[f];
                        ++cursor;
                    }
                }
                n_forest = cursor;
                greened  = 1;
            }
            SMESH_FREE(gleaf);
            SMESH_FREE(gkind);
            SMESH_FREE(gwhich);

            if (!refined && !fanned && !greened) {
                break;
            }
            if (pass + 1 == pass_limit) {
                fprintf(stderr, "refine_adaptive: closure did not finish\n");
                goto adapt_fail;
            }
        }
    }

    {
        ptrdiff_t n_leaves = 0;
        for (ptrdiff_t i = 0; i < n_forest; ++i) {
            if (child0[i] == none) {
                ++n_leaves;
            }
        }
        idx_t *leaves = (idx_t *)SMESH_ALLOC((size_t)(n_leaves > 0 ? n_leaves : 1) * sizeof(idx_t));
        u8    *kind   = (u8 *)SMESH_ALLOC((size_t)(n_leaves > 0 ? n_leaves : 1) * sizeof(u8));
        u8    *which  = (u8 *)SMESH_ALLOC((size_t)(n_leaves > 0 ? n_leaves : 1) * sizeof(u8));
        if (!leaves || !kind || !which) {
            SMESH_FREE(leaves);
            SMESH_FREE(kind);
            SMESH_FREE(which);
            fprintf(stderr, "refine_adaptive: leaf allocation failed\n");
            goto adapt_fail;
        }
        ptrdiff_t nl      = 0;
        ptrdiff_t n_out   = 0;
        for (ptrdiff_t i = 0; i < n_forest; ++i) {
            if (child0[i] != none) {
                continue;
            }
            const int             mask = adapt_leaf_mask<idx_t>(vx, mlo, mhi, n_mids, (idx_t)i);
            const Tet4AdaptAction act  = tet4_adapt_action[mask];
            if (act.kind == TET4_ADAPT_PROMOTE || tet_target[i] > tet_level[i]) {
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: non-conforming leaf remained\n");
                goto adapt_fail;
            }
            leaves[nl] = (idx_t)i;
            kind[nl]   = act.kind;
            which[nl]  = act.which;
            ++nl;
            n_out += (act.kind == TET4_ADAPT_EDGE) ? 2 : (act.kind == TET4_ADAPT_FACE) ? 4 : 1;
        }

        const ptrdiff_t n_face_rows = n_elements * 4;
        const ptrdiff_t n_edge_rows = n_elements * 6;
        ptrdiff_t      *face_count =
                (ptrdiff_t *)SMESH_CALLOC((size_t)(n_face_rows > 0 ? n_face_rows : 1), sizeof(ptrdiff_t));
        ptrdiff_t *edge_count =
                (ptrdiff_t *)SMESH_CALLOC((size_t)(n_edge_rows > 0 ? n_edge_rows : 1), sizeof(ptrdiff_t));
        if (!face_count || !edge_count) {
            SMESH_FREE(face_count);
            SMESH_FREE(edge_count);
            SMESH_FREE(leaves);
            SMESH_FREE(kind);
            SMESH_FREE(which);
            fprintf(stderr, "refine_adaptive: allocation failed\n");
            goto adapt_fail;
        }

        for (ptrdiff_t i = 0; i < n_leaves; ++i) {
            const int rc = adapt_visit_boundary<idx_t>(
                    parent, slot, root, vx, as_lo, as_hi, n_pts, leaves[i], kind[i], which[i],
                    [&](idx_t root_e, int root_f, int) {
                        face_count[(ptrdiff_t)root_e * 4 + root_f] += 1;
                    },
                    [&](idx_t root_e, int root_l, int) {
                        edge_count[(ptrdiff_t)root_e * 6 + root_l] += 1;
                    });
            if (rc != SMESH_SUCCESS) {
                SMESH_FREE(face_count);
                SMESH_FREE(edge_count);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: bad green mask at emit\n");
                goto adapt_fail;
            }
        }

        for (ptrdiff_t r = 0; r < n_face_rows; ++r) {
            if (face_count[r] < 1) {
                SMESH_FREE(face_count);
                SMESH_FREE(edge_count);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: coarse face has no child\n");
                goto adapt_fail;
            }
        }
        for (ptrdiff_t r = 0; r < n_edge_rows; ++r) {
            if (edge_count[r] < 1) {
                SMESH_FREE(face_count);
                SMESH_FREE(edge_count);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: coarse edge has no child\n");
                goto adapt_fail;
            }
        }

        ptrdiff_t n_face_nnz = 0;
        ptrdiff_t n_edge_nnz = 0;
        f_ptr                = (count_t *)SMESH_CALLOC((size_t)n_face_rows + 1, sizeof(count_t));
        e_ptr                = (count_t *)SMESH_CALLOC((size_t)n_edge_rows + 1, sizeof(count_t));
        if (!f_ptr || !e_ptr) {
            SMESH_FREE(face_count);
            SMESH_FREE(edge_count);
            SMESH_FREE(leaves);
            SMESH_FREE(kind);
            SMESH_FREE(which);
            fprintf(stderr, "refine_adaptive: allocation failed\n");
            goto adapt_fail;
        }
        for (ptrdiff_t r = 0; r < n_face_rows; ++r) {
            if (!adapt_count_room<count_t>(n_face_nnz)) {
                SMESH_FREE(face_count);
                SMESH_FREE(edge_count);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: face CRS overflow\n");
                goto adapt_fail;
            }
            f_ptr[r] = (count_t)n_face_nnz;
            n_face_nnz += face_count[r];
        }
        f_ptr[n_face_rows] = (count_t)n_face_nnz;
        for (ptrdiff_t r = 0; r < n_edge_rows; ++r) {
            if (!adapt_count_room<count_t>(n_edge_nnz)) {
                SMESH_FREE(face_count);
                SMESH_FREE(edge_count);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: edge CRS overflow\n");
                goto adapt_fail;
            }
            e_ptr[r] = (count_t)n_edge_nnz;
            n_edge_nnz += edge_count[r];
        }
        e_ptr[n_edge_rows] = (count_t)n_edge_nnz;
        SMESH_FREE(face_count);
        SMESH_FREE(edge_count);
        face_count = nullptr;
        edge_count = nullptr;

        f_elem = (idx_t *)SMESH_ALLOC((size_t)(n_face_nnz > 0 ? n_face_nnz : 1) * sizeof(idx_t));
        f_lfi  = (idx_t *)SMESH_ALLOC((size_t)(n_face_nnz > 0 ? n_face_nnz : 1) * sizeof(idx_t));
        e_elem = (idx_t *)SMESH_ALLOC((size_t)(n_edge_nnz > 0 ? n_edge_nnz : 1) * sizeof(idx_t));
        e_lei  = (idx_t *)SMESH_ALLOC((size_t)(n_edge_nnz > 0 ? n_edge_nnz : 1) * sizeof(idx_t));
        out_e  = (idx_t **)SMESH_ALLOC(4 * sizeof(idx_t *));
        out_p  = (geom_t **)SMESH_ALLOC((size_t)spatial_dim * sizeof(geom_t *));
        if (!f_elem || !f_lfi || !e_elem || !e_lei || !out_e || !out_p) {
            SMESH_FREE(leaves);
            SMESH_FREE(kind);
            SMESH_FREE(which);
            fprintf(stderr, "refine_adaptive: output allocation failed\n");
            goto adapt_fail;
        }
        for (int k = 0; k < 4; ++k) {
            out_e[k] = nullptr;
        }
        for (int d = 0; d < spatial_dim; ++d) {
            out_p[d] = nullptr;
        }
        for (int k = 0; k < 4; ++k) {
            out_e[k] = (idx_t *)SMESH_ALLOC((size_t)(n_out > 0 ? n_out : 1) * sizeof(idx_t));
            if (!out_e[k]) {
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: output allocation failed\n");
                goto adapt_fail;
            }
        }
        if (out_root) {
            out_root_e = (idx_t *)SMESH_ALLOC((size_t)(n_out > 0 ? n_out : 1) * sizeof(idx_t));
            if (!out_root_e) {
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: output allocation failed\n");
                goto adapt_fail;
            }
        }
        for (int d = 0; d < spatial_dim; ++d) {
            out_p[d] = (geom_t *)SMESH_ALLOC((size_t)n_pts * sizeof(geom_t));
            if (!out_p[d]) {
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: output allocation failed\n");
                goto adapt_fail;
            }
            memcpy(out_p[d], pt[d], (size_t)n_pts * sizeof(geom_t));
        }

        ptrdiff_t *face_cursor = (ptrdiff_t *)SMESH_ALLOC((size_t)(n_face_rows > 0 ? n_face_rows : 1) * sizeof(ptrdiff_t));
        ptrdiff_t *edge_cursor = (ptrdiff_t *)SMESH_ALLOC((size_t)(n_edge_rows > 0 ? n_edge_rows : 1) * sizeof(ptrdiff_t));
        if (!face_cursor || !edge_cursor) {
            SMESH_FREE(face_cursor);
            SMESH_FREE(edge_cursor);
            SMESH_FREE(leaves);
            SMESH_FREE(kind);
            SMESH_FREE(which);
            fprintf(stderr, "refine_adaptive: allocation failed\n");
            goto adapt_fail;
        }
        for (ptrdiff_t r = 0; r < n_face_rows; ++r) {
            face_cursor[r] = (ptrdiff_t)f_ptr[r];
        }
        for (ptrdiff_t r = 0; r < n_edge_rows; ++r) {
            edge_cursor[r] = (ptrdiff_t)e_ptr[r];
        }

        ptrdiff_t out_at = 0;
        for (ptrdiff_t i = 0; i < n_leaves; ++i) {
            const idx_t f    = leaves[i];
            const int   nk  = kind[i];
            const int   wh  = which[i];
            int         nc  = 1;
            if (nk == TET4_ADAPT_EDGE) {
                nc = 2;
            } else if (nk == TET4_ADAPT_FACE) {
                nc = 4;
            }
            idx_t macro[10];
            for (int k = 0; k < 4; ++k) {
                macro[k] = vx[k][f];
            }
            for (int e = 0; e < 6; ++e) {
                macro[4 + e] = none;
            }
            if (nk != TET4_ADAPT_KEEP) {
                for (int e = 0; e < 6; ++e) {
                    const ptrdiff_t mi = adapt_find_mid(mlo, mhi, n_mids, macro[tet4_refine_edges[e][0]],
                                                        macro[tet4_refine_edges[e][1]]);
                    if (mi >= 0) {
                        macro[4 + e] = mid[mi];
                    }
                }
            }

            const ptrdiff_t group = out_at;
            for (int g = 0; g < nc; ++g) {
                const int *pat = nullptr;
                int        idpat[4];
                if (nk == TET4_ADAPT_KEEP) {
                    idpat[0] = 0;
                    idpat[1] = 1;
                    idpat[2] = 2;
                    idpat[3] = 3;
                    pat      = idpat;
                } else if (nk == TET4_ADAPT_EDGE) {
                    pat = tet4_green_edge[wh][g];
                } else {
                    pat = tet4_green_face[wh][g];
                }
                for (int k = 0; k < 4; ++k) {
                    const int src = pat[k];
                    if (src < 0 || macro[src] == none) {
                        SMESH_FREE(face_cursor);
                        SMESH_FREE(edge_cursor);
                        SMESH_FREE(leaves);
                        SMESH_FREE(kind);
                        SMESH_FREE(which);
                        fprintf(stderr, "refine_adaptive: green pattern referenced a missing node\n");
                        goto adapt_fail;
                    }
                    out_e[k][group + g] = macro[src];
                }
            }
            if (out_root_e) {
                const idx_t coarse = root[f];
                for (int g = 0; g < nc; ++g) {
                    out_root_e[group + g] = coarse;
                }
            }

            for (int g = 0; g < nc; ++g) {
                for (int lf = 0; lf < 4; ++lf) {
                    int leaf_face = lf;
                    if (nk == TET4_ADAPT_EDGE) {
                        leaf_face = tet4_green_edge_face[wh][g][lf];
                    } else if (nk == TET4_ADAPT_FACE) {
                        leaf_face = tet4_green_face_face[wh][g][lf];
                    }
                    if (leaf_face < 0) {
                        continue;
                    }
                    idx_t root_e = none;
                    int   root_f = -1;
                    if (!adapt_map_face(parent, slot, root, f, leaf_face, &root_e, &root_f)) {
                        continue;
                    }
                    const ptrdiff_t row = (ptrdiff_t)root_e * 4 + root_f;
                    if (row < 0 || row >= n_face_rows) {
                        SMESH_FREE(face_cursor);
                        SMESH_FREE(edge_cursor);
                        SMESH_FREE(leaves);
                        SMESH_FREE(kind);
                        SMESH_FREE(which);
                        fprintf(stderr, "refine_adaptive: face child row out of range\n");
                        goto adapt_fail;
                    }
                    const ptrdiff_t at = face_cursor[row]++;
                    if (at < (ptrdiff_t)f_ptr[row] || at >= (ptrdiff_t)f_ptr[row + 1]) {
                        SMESH_FREE(face_cursor);
                        SMESH_FREE(edge_cursor);
                        SMESH_FREE(leaves);
                        SMESH_FREE(kind);
                        SMESH_FREE(which);
                        fprintf(stderr, "refine_adaptive: face child CRS overrun\n");
                        goto adapt_fail;
                    }
                    f_elem[at] = (idx_t)(group + g);
                    f_lfi[at]  = (idx_t)lf;
                }
                for (int lei = 0; lei < 6; ++lei) {
                    int leaf_edge = lei;
                    if (nk == TET4_ADAPT_EDGE) {
                        leaf_edge = tet4_green_edge_edge[wh][g][lei];
                    } else if (nk == TET4_ADAPT_FACE) {
                        leaf_edge = tet4_green_face_edge[wh][g][lei];
                    }
                    if (leaf_edge < 0) {
                        continue;
                    }
                    idx_t root_e = none;
                    int   root_l = -1;
                    if (!adapt_map_edge(parent, slot, root, vx, as_lo, as_hi, n_pts, f, leaf_edge, &root_e, &root_l)) {
                        continue;
                    }
                    const ptrdiff_t row = (ptrdiff_t)root_e * 6 + root_l;
                    if (row < 0 || row >= n_edge_rows) {
                        SMESH_FREE(face_cursor);
                        SMESH_FREE(edge_cursor);
                        SMESH_FREE(leaves);
                        SMESH_FREE(kind);
                        SMESH_FREE(which);
                        fprintf(stderr, "refine_adaptive: edge child row out of range\n");
                        goto adapt_fail;
                    }
                    const ptrdiff_t at = edge_cursor[row]++;
                    if (at < (ptrdiff_t)e_ptr[row] || at >= (ptrdiff_t)e_ptr[row + 1]) {
                        SMESH_FREE(face_cursor);
                        SMESH_FREE(edge_cursor);
                        SMESH_FREE(leaves);
                        SMESH_FREE(kind);
                        SMESH_FREE(which);
                        fprintf(stderr, "refine_adaptive: edge child CRS overrun\n");
                        goto adapt_fail;
                    }
                    e_elem[at] = (idx_t)(group + g);
                    e_lei[at]  = (idx_t)lei;
                }
            }
            out_at = group + nc;
        }
        for (ptrdiff_t r = 0; r < n_face_rows; ++r) {
            if (face_cursor[r] != (ptrdiff_t)f_ptr[r + 1]) {
                SMESH_FREE(face_cursor);
                SMESH_FREE(edge_cursor);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: face child count mismatch\n");
                goto adapt_fail;
            }
        }
        for (ptrdiff_t r = 0; r < n_edge_rows; ++r) {
            if (edge_cursor[r] != (ptrdiff_t)e_ptr[r + 1]) {
                SMESH_FREE(face_cursor);
                SMESH_FREE(edge_cursor);
                SMESH_FREE(leaves);
                SMESH_FREE(kind);
                SMESH_FREE(which);
                fprintf(stderr, "refine_adaptive: edge child count mismatch\n");
                goto adapt_fail;
            }
        }

        SMESH_FREE(face_cursor);
        SMESH_FREE(edge_cursor);
        SMESH_FREE(leaves);
        SMESH_FREE(kind);
        SMESH_FREE(which);

        *out_elements   = out_e;
        *n_out_elements = n_out;
        *out_points     = out_p;
        *n_out_nodes    = n_pts;
        *face_ptr       = f_ptr;
        *face_elem      = f_elem;
        *face_lfi       = f_lfi;
        *edge_ptr       = e_ptr;
        *edge_elem      = e_elem;
        *edge_lei       = e_lei;
        if (out_root) {
            *out_root = out_root_e;
            out_root_e = nullptr;
        }
        if (n_created > 0) {
            *mid_lo = clo;
            *mid_hi = chi;
            *mid_id = cid;
            clo     = nullptr;
            chi     = nullptr;
            cid     = nullptr;
        }
        *n_mids_out = n_created;
        out_e       = nullptr;
        out_p       = nullptr;
        f_ptr       = nullptr;
        f_elem      = nullptr;
        f_lfi       = nullptr;
        e_ptr       = nullptr;
        e_elem      = nullptr;
        e_lei       = nullptr;
    }

    for (int k = 0; k < 4; ++k) {
        SMESH_FREE(vx[k]);
    }
    SMESH_FREE(tet_level);
    SMESH_FREE(tet_target);
    SMESH_FREE(child0);
    SMESH_FREE(parent);
    SMESH_FREE(slot);
    SMESH_FREE(root);
    for (int d = 0; d < spatial_dim; ++d) {
        SMESH_FREE(pt[d]);
    }
    SMESH_FREE(as_lo);
    SMESH_FREE(as_hi);
    SMESH_FREE(mlo);
    SMESH_FREE(mhi);
    SMESH_FREE(mid);
    SMESH_FREE(clo);
    SMESH_FREE(chi);
    SMESH_FREE(cid);
    return SMESH_SUCCESS;

adapt_fail:
    for (int k = 0; k < 4; ++k) {
        SMESH_FREE(vx[k]);
    }
    SMESH_FREE(tet_level);
    SMESH_FREE(tet_target);
    SMESH_FREE(child0);
    SMESH_FREE(parent);
    SMESH_FREE(slot);
    SMESH_FREE(root);
    for (int d = 0; d < 3; ++d) {
        SMESH_FREE(pt[d]);
    }
    SMESH_FREE(as_lo);
    SMESH_FREE(as_hi);
    SMESH_FREE(mlo);
    SMESH_FREE(mhi);
    SMESH_FREE(mid);
    SMESH_FREE(clo);
    SMESH_FREE(chi);
    SMESH_FREE(cid);
    if (out_e) {
        for (int k = 0; k < 4; ++k) {
            SMESH_FREE(out_e[k]);
        }
        SMESH_FREE(out_e);
    }
    if (out_p) {
        for (int d = 0; d < spatial_dim; ++d) {
            SMESH_FREE(out_p[d]);
        }
        SMESH_FREE(out_p);
    }
    SMESH_FREE(f_ptr);
    SMESH_FREE(f_elem);
    SMESH_FREE(f_lfi);
    SMESH_FREE(e_ptr);
    SMESH_FREE(e_elem);
    SMESH_FREE(e_lei);
    SMESH_FREE(out_root_e);
    return SMESH_FAILURE;
}

}  // namespace smesh

#endif  // SMESH_REFINE_ADAPTIVE_IMPL_HPP
