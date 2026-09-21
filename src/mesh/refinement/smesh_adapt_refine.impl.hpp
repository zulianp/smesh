#ifndef SMESH_ADAPT_REFINE_IMPL_HPP
#define SMESH_ADAPT_REFINE_IMPL_HPP

#include "smesh_adapt_refine.hpp"

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_graph.hpp"
#include "smesh_quality.hpp"
#include "smesh_search.hpp"

#include <cmath>
#include <cstdlib>
#include <string.h>

namespace smesh {

namespace {

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

/// Size field pairs that participate in geometry marking. Flat/interior nodes
/// keep a sentinel h ≫ curved h; mixed pairs must not mark (that is ℓ-based
/// volume fill). 2:1 neighbors may differ by at most 2, so 32× is unconstrained.
template <typename geom_t>
int h_both_constrained(const geom_t ha, const geom_t hb) {
    const geom_t lo = ha < hb ? ha : hb;
    const geom_t hi = ha < hb ? hb : ha;
    return lo > static_cast<geom_t>(0) && hi < lo * static_cast<geom_t>(32);
}

template <typename idx_t>
int local_edge_of(const LocalEdgeTable &let, const int n_lei, const idx_t *const v, const idx_t a, const idx_t b) {
    for (int le = 0; le < n_lei; ++le) {
        const idx_t x = v[let(le, 0)];
        const idx_t y = v[let(le, 1)];
        if ((x == a && y == b) || (x == b && y == a)) {
            return le;
        }
    }
    return 0;
}

template <typename idx_t, typename geom_t>
geom_t tri_orient3(const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p,
                   const int                                                sdim,
                   const idx_t                                              a,
                   const idx_t                                              b,
                   const idx_t                                              c) {
    const geom_t ux = p[0][b] - p[0][a], uy = p[1][b] - p[1][a];
    const geom_t vx = p[0][c] - p[0][a], vy = p[1][c] - p[1][a];
    if (sdim < 3) {
        return ux * vy - uy * vx;
    }
    const geom_t uz = p[2][b] - p[2][a], vz = p[2][c] - p[2][a];
    const geom_t nx = uy * vz - uz * vy;
    const geom_t ny = uz * vx - ux * vz;
    const geom_t nz = ux * vy - uy * vx;
    return nx * nx + ny * ny + nz * nz > static_cast<geom_t>(0) ? nz : ux * vy - uy * vx;
}

template <typename idx_t, typename geom_t>
void fix_tri_orient(const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT p,
                    const int                                                sdim,
                    idx_t                                                   *c,
                    const geom_t                                             oref) {
    const geom_t o = tri_orient3(p, sdim, c[0], c[1], c[2]);
    if (o * oref < static_cast<geom_t>(0)) {
        const idx_t t = c[1];
        c[1]          = c[2];
        c[2]          = t;
    }
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

static const int kQuad4Child[4][4] = {{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};

struct UniqueEdgePri {
    double    l2;
    ptrdiff_t id;
};

int unique_edge_pri_desc(const void *a, const void *b) {
    const UniqueEdgePri *x = static_cast<const UniqueEdgePri *>(a);
    const UniqueEdgePri *y = static_cast<const UniqueEdgePri *>(b);
    if (x->l2 < y->l2) {
        return 1;
    }
    if (x->l2 > y->l2) {
        return -1;
    }
    return (x->id > y->id) - (x->id < y->id);
}

}  // namespace

template <typename idx_t, typename geom_t>
void adapt_refine_free_t(const int         nxe,
                         const int         sdim,
                         idx_t           **elements,
                         geom_t          **points,
                         ptrdiff_t        *parent_elem,
                         void             *parent_ptr,
                         idx_t            *child_id,
                         idx_t            *node_a,
                         idx_t            *node_b) {
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
    SMESH_FREE(parent_elem);
    SMESH_FREE(parent_ptr);
    SMESH_FREE(child_id);
    SMESH_FREE(node_a);
    SMESH_FREE(node_b);
}

void mesh_adapt_refine_free(const int         nxe,
                            const int         sdim,
                            idx_t           **elements,
                            geom_t          **points,
                            ptrdiff_t        *parent_elem,
                            void             *parent_ptr,
                            idx_t            *child_id,
                            idx_t            *node_a,
                            idx_t            *node_b) {
    adapt_refine_free_t<idx_t, geom_t>(nxe, sdim, elements, points, parent_elem, parent_ptr, child_id, node_a, node_b);
}

template <typename idx_t, typename count_t, typename geom_t>
int mesh_adapt_refine(const enum ElemType                                     element_type,
                      const ptrdiff_t                                         n_elements_in,
                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements_in,
                      const int                                               sdim,
                      const ptrdiff_t                                         n_nodes_in,
                      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points_in,
                      const geom_t *const SMESH_RESTRICT                      h,
                      const uint8_t *const SMESH_RESTRICT                     element_mark,
                      const geom_t                                            q_min,
                      const uint8_t *const SMESH_RESTRICT                     geom_node,
                      const int                                               max_levels,
                      ptrdiff_t                                              *n_elements_out,
                      idx_t                                                ***elements_out,
                      ptrdiff_t                                              *n_nodes_out,
                      geom_t                                               ***points_out,
                      ptrdiff_t                                             **parent_elem_out,
                      count_t                                               **parent_ptr_out,
                      idx_t                                                 **child_id_out,
                      idx_t                                                 **node_a_out,
                      idx_t                                                 **node_b_out) {
    if (!adapt_refine_type_supported(element_type) || !elements_in || !points_in || !h ||
        !n_elements_out || !elements_out || !n_nodes_out || !points_out || !parent_elem_out ||
        !parent_ptr_out || !child_id_out || !node_a_out || !node_b_out || sdim < 2) {
        return SMESH_FAILURE;
    }
    const int nxe   = elem_num_nodes(element_type);
    const int n_lei = elem_num_edges(element_type);
    const bool is_quad =
            element_type == QUAD4 || element_type == QUADSHELL4;
    const bool is_tet = element_type == TET4;
    LocalEdgeTable let;
    if (let.fill(element_type) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }

    ptrdiff_t n_elem  = n_elements_in;
    ptrdiff_t n_nodes = n_nodes_in;
    idx_t   **elems   = (idx_t **)SMESH_ALLOC((size_t)nxe * sizeof(idx_t *));
    geom_t  **pts     = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
    if (!elems || !pts) {
        SMESH_FREE(elems);
        SMESH_FREE(pts);
        return SMESH_FAILURE;
    }
    for (int d = 0; d < nxe; ++d) {
        elems[d] = (idx_t *)SMESH_ALLOC((size_t)n_elem * sizeof(idx_t));
        if (!elems[d]) {
            adapt_refine_free_t<idx_t, geom_t>(d, 0, elems, (geom_t **)nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
            SMESH_FREE(pts);
            return SMESH_FAILURE;
        }
        memcpy(elems[d], elements_in[d], (size_t)n_elem * sizeof(idx_t));
    }
    for (int d = 0; d < sdim; ++d) {
        pts[d] = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
        if (!pts[d]) {
            adapt_refine_free_t<idx_t, geom_t>(nxe, d, elems, pts, nullptr, nullptr, nullptr, nullptr, nullptr);
            return SMESH_FAILURE;
        }
        memcpy(pts[d], points_in[d], (size_t)n_nodes * sizeof(geom_t));
    }
    ptrdiff_t *parent = (ptrdiff_t *)SMESH_ALLOC((size_t)n_elem * sizeof(ptrdiff_t));
    idx_t     *node_a = (idx_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(idx_t));
    idx_t     *node_b = (idx_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(idx_t));
    uint8_t   *force  = (uint8_t *)SMESH_CALLOC((size_t)n_elem, sizeof(uint8_t));
    uint8_t   *lei    = (uint8_t *)SMESH_ALLOC((size_t)n_elem * sizeof(uint8_t));
    if (!parent || !node_a || !node_b || !force || !lei) {
        adapt_refine_free_t<idx_t, geom_t>(nxe, sdim, elems, pts, parent, nullptr, nullptr, node_a, node_b);
        SMESH_FREE(force);
        SMESH_FREE(lei);
        return SMESH_FAILURE;
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t e = 0; e < n_elem; ++e) {
        parent[e] = e;
        force[e]  = (element_mark && element_mark[e]) ? 1 : 0;
        geom_t best = -1;
        int    ble  = 0;
        for (int le = 0; le < n_lei; ++le) {
            const idx_t  a  = elems[let(le, 0)][e];
            const idx_t  b  = elems[let(le, 1)][e];
            const geom_t l2 = edge_len2(sdim, pts, a, b);
            if (l2 > best) {
                best = l2;
                ble  = le;
            }
        }
        lei[e] = (uint8_t)ble;
    }
#pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        node_a[i] = (idx_t)i;
        node_b[i] = (idx_t)i;
    }

    const int levels = max_levels > 0 ? max_levels : 1;
    geom_t   *h_work = (geom_t *)SMESH_ALLOC((size_t)n_nodes * sizeof(geom_t));
    if (!h_work) {
        adapt_refine_free_t<idx_t, geom_t>(nxe, sdim, elems, pts, parent, nullptr, nullptr, node_a, node_b);
        SMESH_FREE(force);
        SMESH_FREE(lei);
        return SMESH_FAILURE;
    }
    memcpy(h_work, h, (size_t)n_nodes * sizeof(geom_t));
    uint8_t *g_work = (uint8_t *)SMESH_CALLOC((size_t)n_nodes, sizeof(uint8_t));
    if (!g_work) {
        adapt_refine_free_t<idx_t, geom_t>(nxe, sdim, elems, pts, parent, nullptr, nullptr, node_a, node_b);
        SMESH_FREE(force);
        SMESH_FREE(lei);
        SMESH_FREE(h_work);
        return SMESH_FAILURE;
    }
    if (geom_node) {
        memcpy(g_work, geom_node, (size_t)n_nodes);
    }

    for (int wave = 0; wave < levels; ++wave) {
        count_t       *n2eptr  = nullptr;
        element_idx_t *elindex = nullptr;
        count_t       *rowptr  = nullptr;
        idx_t         *colidx  = nullptr;
        if (create_n2e<idx_t, count_t, element_idx_t>(
                    n_elem, n_nodes, nxe, elems, &n2eptr, &elindex) != SMESH_SUCCESS ||
            create_edge_graph_for_element_from_n2e<idx_t, count_t>(
                    element_type, n_elem, n_nodes, elems, n2eptr, elindex, &rowptr, &colidx) !=
                    SMESH_SUCCESS) {
            SMESH_FREE(n2eptr);
            SMESH_FREE(elindex);
            SMESH_FREE(h_work);
            SMESH_FREE(g_work);
            SMESH_FREE(force);
            SMESH_FREE(lei);
            adapt_refine_free_t<idx_t, geom_t>(nxe, sdim, elems, pts, parent, nullptr, nullptr, node_a, node_b);
            return SMESH_FAILURE;
        }
        SMESH_FREE(n2eptr);
        SMESH_FREE(elindex);
        const count_t nnz = rowptr[n_nodes];
        ptrdiff_t    *uid = (ptrdiff_t *)SMESH_ALLOC((size_t)nnz * sizeof(ptrdiff_t));
#pragma omp parallel for schedule(static)
        for (count_t k = 0; k < nnz; ++k) {
            uid[k] = -1;
        }
        ptrdiff_t n_uedge = 0;
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
                if ((ptrdiff_t)colidx[k] > i) {
                    uid[k] = n_uedge++;
                }
            }
        }

        uint8_t *eflag = (uint8_t *)SMESH_CALLOC((size_t)n_uedge, sizeof(uint8_t));
        uint8_t *emark = (uint8_t *)SMESH_CALLOC((size_t)n_elem, sizeof(uint8_t));
        int     *longest = (int *)SMESH_ALLOC((size_t)n_elem * sizeof(int));
        ptrdiff_t *long_uid = (ptrdiff_t *)SMESH_ALLOC((size_t)n_elem * sizeof(ptrdiff_t));
        if (!uid || !eflag || !emark || !longest || !long_uid) {
            SMESH_FREE(rowptr);
            SMESH_FREE(colidx);
            SMESH_FREE(uid);
            SMESH_FREE(eflag);
            SMESH_FREE(emark);
            SMESH_FREE(longest);
            SMESH_FREE(long_uid);
            SMESH_FREE(h_work);
            SMESH_FREE(g_work);
            SMESH_FREE(force);
            SMESH_FREE(lei);
            adapt_refine_free_t<idx_t, geom_t>(nxe, sdim, elems, pts, parent, nullptr, nullptr, node_a, node_b);
            return SMESH_FAILURE;
        }

#pragma omp parallel for schedule(static)
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
            geom_t best = -1;
            int    ble  = 0;
            for (int le = 0; le < n_lei; ++le) {
                const idx_t  a  = elems[let(le, 0)][e];
                const idx_t  b  = elems[let(le, 1)][e];
                const geom_t l2 = edge_len2(sdim, pts, a, b);
                if (l2 > best) {
                    best = l2;
                    ble  = le;
                }
            }
            lei[e] = (uint8_t)ble;
        }

#pragma omp parallel for schedule(static)
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
            const int le_ref = (int)lei[e];
            int       need   = force[e];
            ptrdiff_t buid   = -1;
            for (int le = 0; le < n_lei; ++le) {
                const idx_t a = elems[let(le, 0)][e];
                const idx_t b = elems[let(le, 1)][e];
                const geom_t l2 = edge_len2(sdim, pts, a, b);
                const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                const ptrdiff_t id = (slot == (count_t)(-1)) ? (ptrdiff_t)-1 : uid[slot];
                if (le == le_ref) {
                    buid = id;
                }
                const geom_t ha = h_work[a];
                const geom_t hb = h_work[b];
                const geom_t hi = ha < hb ? ha : hb;
                if (h_both_constrained(ha, hb) && l2 > hi * hi) {
                    need = 1;
                }
            }
            longest[e]  = le_ref;
            long_uid[e] = buid;
            emark[e]    = (uint8_t)need;
        }

        if (is_quad) {
            int changed = 1;
            while (changed) {
                changed = 0;
                for (ptrdiff_t e = 0; e < n_elem; ++e) {
                    if (!emark[e]) {
                        continue;
                    }
                    for (int le = 0; le < n_lei; ++le) {
                        const idx_t a = elems[let(le, 0)][e];
                        const idx_t b = elems[let(le, 1)][e];
                        const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                        if (slot != (count_t)(-1) && uid[slot] >= 0) {
                            eflag[uid[slot]] = 1;
                        }
                    }
                }
                for (ptrdiff_t e = 0; e < n_elem; ++e) {
                    if (emark[e]) {
                        continue;
                    }
                    for (int le = 0; le < n_lei; ++le) {
                        const idx_t a = elems[let(le, 0)][e];
                        const idx_t b = elems[let(le, 1)][e];
                        const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                        if (slot != (count_t)(-1) && uid[slot] >= 0 && eflag[uid[slot]]) {
                            emark[e] = 1;
                            changed  = 1;
                            break;
                        }
                    }
                }
            }
        } else {
            geom_t *q = (geom_t *)SMESH_ALLOC((size_t)n_elem * sizeof(geom_t));
            if (q) {
                mesh_element_quality<idx_t, geom_t>(element_type, n_elem, elems, sdim, pts, q);
            }
            for (ptrdiff_t i = 0; i < n_uedge; ++i) {
                eflag[i] = 0;
            }
            for (ptrdiff_t e = 0; e < n_elem; ++e) {
                int has_size = force[e];
                int has_geom = 0;
                for (int le = 0; le < n_lei; ++le) {
                    const idx_t a = elems[let(le, 0)][e];
                    const idx_t b = elems[let(le, 1)][e];
                    const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                    if (slot == (count_t)(-1) || uid[slot] < 0) {
                        continue;
                    }
                    const geom_t l2 = edge_len2(sdim, pts, a, b);
                    const geom_t ha = h_work[a];
                    const geom_t hb = h_work[b];
                    const geom_t hi = ha < hb ? ha : hb;
                    if (h_both_constrained(ha, hb) && l2 > hi * hi) {
                        has_size = 1;
                        if (g_work[a] && g_work[b]) {
                            eflag[uid[slot]] = 1;
                            has_geom         = 1;
                        }
                    }
                }
                const int poor = q && q_min > static_cast<geom_t>(0) && q[e] < q_min;
                if (force[e] || poor || (has_size && !has_geom)) {
                    const ptrdiff_t id = long_uid[e];
                    if (id >= 0) {
                        eflag[id] = 1;
                    }
                }
            }
            SMESH_FREE(q);
        }

        ptrdiff_t *mid_of = (ptrdiff_t *)SMESH_ALLOC((size_t)n_uedge * sizeof(ptrdiff_t));
#pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < n_uedge; ++i) {
            mid_of[i] = -1;
        }

        ptrdiff_t n_split_edges = 0;
        if (is_quad) {
            for (ptrdiff_t i = 0; i < n_uedge; ++i) {
                if (eflag[i]) {
                    mid_of[i] = n_nodes + n_split_edges;
                    ++n_split_edges;
                }
            }
        } else {
            UniqueEdgePri *ord = (UniqueEdgePri *)SMESH_ALLOC((size_t)n_uedge * sizeof(UniqueEdgePri));
            ptrdiff_t      n_marked = 0;
            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
                    if ((ptrdiff_t)colidx[k] <= i) {
                        continue;
                    }
                    const ptrdiff_t id = uid[k];
                    if (id < 0 || !eflag[id]) {
                        continue;
                    }
                    ord[n_marked].l2 = (double)edge_len2(sdim, pts, (idx_t)i, colidx[k]);
                    ord[n_marked].id = id;
                    ++n_marked;
                }
            }
            if (n_marked > 1) {
                qsort(ord, (size_t)n_marked, sizeof(UniqueEdgePri), unique_edge_pri_desc);
            }

            count_t *e2tptr = (count_t *)SMESH_CALLOC((size_t)(n_uedge + 1), sizeof(count_t));
            for (ptrdiff_t e = 0; e < n_elem; ++e) {
                for (int le = 0; le < n_lei; ++le) {
                    const idx_t a = elems[let(le, 0)][e];
                    const idx_t b = elems[let(le, 1)][e];
                    const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                    if (slot != (count_t)(-1) && uid[slot] >= 0) {
                        e2tptr[uid[slot] + 1] += 1;
                    }
                }
            }
            for (ptrdiff_t i = 0; i < n_uedge; ++i) {
                e2tptr[i + 1] += e2tptr[i];
            }
            ptrdiff_t *e2t  = (ptrdiff_t *)SMESH_ALLOC((size_t)e2tptr[n_uedge] * sizeof(ptrdiff_t));
            count_t   *fill = (count_t *)SMESH_ALLOC((size_t)(n_uedge + 1) * sizeof(count_t));
            memcpy(fill, e2tptr, (size_t)(n_uedge + 1) * sizeof(count_t));
            for (ptrdiff_t e = 0; e < n_elem; ++e) {
                for (int le = 0; le < n_lei; ++le) {
                    const idx_t a = elems[let(le, 0)][e];
                    const idx_t b = elems[let(le, 1)][e];
                    const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                    if (slot != (count_t)(-1) && uid[slot] >= 0) {
                        e2t[fill[uid[slot]]++] = e;
                    }
                }
            }
            SMESH_FREE(fill);

            ptrdiff_t *assign = (ptrdiff_t *)SMESH_ALLOC((size_t)n_elem * sizeof(ptrdiff_t));
            for (ptrdiff_t e = 0; e < n_elem; ++e) {
                assign[e] = -1;
            }
            for (ptrdiff_t k = 0; k < n_marked; ++k) {
                const ptrdiff_t id = ord[k].id;
                int             busy = 0;
                for (count_t t = e2tptr[id]; t < e2tptr[id + 1]; ++t) {
                    if (assign[e2t[t]] >= 0) {
                        busy = 1;
                        break;
                    }
                }
                if (busy) {
                    continue;
                }
                mid_of[id] = n_nodes + n_split_edges;
                ++n_split_edges;
                for (count_t t = e2tptr[id]; t < e2tptr[id + 1]; ++t) {
                    assign[e2t[t]] = id;
                }
            }

            for (ptrdiff_t e = 0; e < n_elem; ++e) {
                longest[e]  = 0;
                long_uid[e] = -1;
                if (assign[e] < 0) {
                    continue;
                }
                long_uid[e] = assign[e];
                for (int le = 0; le < n_lei; ++le) {
                    const idx_t a = elems[let(le, 0)][e];
                    const idx_t b = elems[let(le, 1)][e];
                    const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                    if (slot != (count_t)(-1) && uid[slot] == assign[e]) {
                        longest[e] = le;
                        break;
                    }
                }
            }

            SMESH_FREE(ord);
            SMESH_FREE(e2tptr);
            SMESH_FREE(e2t);
            SMESH_FREE(assign);
        }

        ptrdiff_t n_split_elem = 0;
        uint8_t *split = (uint8_t *)SMESH_CALLOC((size_t)n_elem, sizeof(uint8_t));
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
            if (is_quad) {
                if (emark[e]) {
                    split[e] = 1;
                    ++n_split_elem;
                }
            } else if (long_uid[e] >= 0 && mid_of[long_uid[e]] >= 0) {
                split[e] = 1;
                ++n_split_elem;
            }
        }

        if (n_split_elem == 0) {
            SMESH_FREE(rowptr);
            SMESH_FREE(colidx);
            SMESH_FREE(uid);
            SMESH_FREE(eflag);
            SMESH_FREE(emark);
            SMESH_FREE(longest);
            SMESH_FREE(long_uid);
            SMESH_FREE(mid_of);
            SMESH_FREE(split);
            break;
        }

        const ptrdiff_t n_fc = is_quad ? n_split_elem : 0;
        const ptrdiff_t n_nodes_new =
                n_nodes + n_split_edges + n_fc;
        const ptrdiff_t n_elem_new =
                is_quad ? (n_elem - n_split_elem + 4 * n_split_elem)
                        : (n_elem - n_split_elem + 2 * n_split_elem);

        idx_t **nelems = (idx_t **)SMESH_ALLOC((size_t)nxe * sizeof(idx_t *));
        for (int d = 0; d < nxe; ++d) {
            nelems[d] = (idx_t *)SMESH_ALLOC((size_t)n_elem_new * sizeof(idx_t));
        }
        geom_t **npts = (geom_t **)SMESH_ALLOC((size_t)sdim * sizeof(geom_t *));
        for (int d = 0; d < sdim; ++d) {
            npts[d] = (geom_t *)SMESH_ALLOC((size_t)n_nodes_new * sizeof(geom_t));
            memcpy(npts[d], pts[d], (size_t)n_nodes * sizeof(geom_t));
        }
        ptrdiff_t *nparent = (ptrdiff_t *)SMESH_ALLOC((size_t)n_elem_new * sizeof(ptrdiff_t));
        uint8_t   *nforce  = (uint8_t *)SMESH_CALLOC((size_t)n_elem_new, sizeof(uint8_t));
        uint8_t   *nlei    = (uint8_t *)SMESH_ALLOC((size_t)n_elem_new * sizeof(uint8_t));
        idx_t     *nna     = (idx_t *)SMESH_ALLOC((size_t)n_nodes_new * sizeof(idx_t));
        idx_t     *nnb     = (idx_t *)SMESH_ALLOC((size_t)n_nodes_new * sizeof(idx_t));
        memcpy(nna, node_a, (size_t)n_nodes * sizeof(idx_t));
        memcpy(nnb, node_b, (size_t)n_nodes * sizeof(idx_t));

        // Edge mids in unique-edge order.
        {
            ptrdiff_t written = 0;
            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
                    if ((ptrdiff_t)colidx[k] <= i) {
                        continue;
                    }
                    const ptrdiff_t id = uid[k];
                    if (id < 0 || mid_of[id] < 0) {
                        continue;
                    }
                    const idx_t j = colidx[k];
                    const idx_t mid = (idx_t)mid_of[id];
                    for (int d = 0; d < sdim; ++d) {
                        npts[d][mid] = static_cast<geom_t>(0.5) * (pts[d][i] + pts[d][j]);
                    }
                    nna[mid] = (idx_t)i;
                    nnb[mid] = j;
                    ++written;
                }
            }
            (void)written;
        }

        geom_t *h_new = (geom_t *)SMESH_ALLOC((size_t)n_nodes_new * sizeof(geom_t));
        uint8_t *g_new = (uint8_t *)SMESH_CALLOC((size_t)n_nodes_new, sizeof(uint8_t));
        memcpy(h_new, h_work, (size_t)n_nodes * sizeof(geom_t));
        memcpy(g_new, g_work, (size_t)n_nodes);
        for (ptrdiff_t i = 0; i < n_uedge; ++i) {
            if (mid_of[i] < 0) {
                continue;
            }
        }
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            for (count_t k = rowptr[i]; k < rowptr[i + 1]; ++k) {
                if ((ptrdiff_t)colidx[k] <= i) {
                    continue;
                }
                const ptrdiff_t id = uid[k];
                if (id < 0 || mid_of[id] < 0) {
                    continue;
                }
                const idx_t j   = colidx[k];
                const idx_t mid = (idx_t)mid_of[id];
                h_new[mid]      = h_work[i] < h_work[j] ? h_work[i] : h_work[j];
                g_new[mid]      = (uint8_t)(g_work[i] && g_work[j]);
            }
        }

        ptrdiff_t w  = 0;
        ptrdiff_t fc = 0;
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
            if (!split[e]) {
                for (int d = 0; d < nxe; ++d) {
                    nelems[d][w] = elems[d][e];
                }
                nparent[w] = parent[e];
                nforce[w]  = force[e];
                nlei[w]    = lei[e];
                ++w;
                continue;
            }
            if (is_quad) {
                idx_t macro[9];
                for (int k = 0; k < 4; ++k) {
                    macro[k] = elems[k][e];
                }
                for (int le = 0; le < 4; ++le) {
                    const idx_t a = elems[let(le, 0)][e];
                    const idx_t b = elems[let(le, 1)][e];
                    const count_t slot = find_n2n_slot(rowptr, colidx, a, b);
                    macro[4 + le]      = (idx_t)mid_of[uid[slot]];
                }
                const idx_t fcid = (idx_t)(n_nodes + n_split_edges + fc);
                macro[8]         = fcid;
                for (int d = 0; d < sdim; ++d) {
                    npts[d][fcid] = static_cast<geom_t>(0.25) *
                                    (pts[d][macro[0]] + pts[d][macro[1]] + pts[d][macro[2]] +
                                     pts[d][macro[3]]);
                }
                nna[fcid]    = macro[0];
                nnb[fcid]    = macro[2];
                h_new[fcid]  = h_work[macro[0]];
                g_new[fcid]  = g_work[macro[0]];
                for (int k = 1; k < 4; ++k) {
                    if (h_work[macro[k]] < h_new[fcid]) {
                        h_new[fcid] = h_work[macro[k]];
                    }
                    g_new[fcid] = (uint8_t)(g_new[fcid] && g_work[macro[k]]);
                }
                ++fc;
                for (int c = 0; c < 4; ++c) {
                    for (int d = 0; d < 4; ++d) {
                        nelems[d][w] = macro[kQuad4Child[c][d]];
                    }
                    nparent[w] = parent[e];
                    nforce[w]  = force[e];
                    nlei[w]    = 0;
                    ++w;
                }
            } else {
                const int   le = longest[e];
                const idx_t a  = elems[let(le, 0)][e];
                const idx_t b  = elems[let(le, 1)][e];
                const idx_t m  = (idx_t)mid_of[long_uid[e]];
                if (is_tet) {
                    idx_t v[4];
                    int   nv = 0;
                    for (int d = 0; d < 4; ++d) {
                        const idx_t q = elems[d][e];
                        if (q != a && q != b) {
                            v[nv++] = q;
                        }
                    }
                    if (nv != 2) {
                        SMESH_ERROR("adapt_refine: tet bisection expected 2 opposite nodes, got %d\n", nv);
                        return SMESH_FAILURE;
                    }
                    idx_t c0[4] = {a, m, v[0], v[1]};
                    idx_t c1[4] = {m, b, v[0], v[1]};
                    if (sdim >= 3) {
                        if (tet_orient(npts, c0[0], c0[1], c0[2], c0[3]) < static_cast<geom_t>(0)) {
                            const idx_t t = c0[2];
                            c0[2]         = c0[3];
                            c0[3]         = t;
                        }
                        if (tet_orient(npts, c1[0], c1[1], c1[2], c1[3]) < static_cast<geom_t>(0)) {
                            const idx_t t = c1[2];
                            c1[2]         = c1[3];
                            c1[3]         = t;
                        }
                    }
                    const int lei0 = local_edge_of(let, n_lei, c0, v[0], v[1]);
                    const int lei1 = local_edge_of(let, n_lei, c1, v[0], v[1]);
                    for (int d = 0; d < 4; ++d) {
                        nelems[d][w] = c0[d];
                    }
                    nparent[w] = parent[e];
                    nforce[w]  = force[e];
                    nlei[w]    = (uint8_t)lei0;
                    ++w;
                    for (int d = 0; d < 4; ++d) {
                        nelems[d][w] = c1[d];
                    }
                    nparent[w] = parent[e];
                    nforce[w]  = force[e];
                    nlei[w]    = (uint8_t)lei1;
                    ++w;
                } else {
                    idx_t other = 0;
                    for (int d = 0; d < 3; ++d) {
                        const idx_t q = elems[d][e];
                        if (q != a && q != b) {
                            other = q;
                        }
                    }
                    const geom_t oref = tri_orient3(pts, sdim, elems[0][e], elems[1][e], elems[2][e]);
                    idx_t c0[3] = {a, m, other};
                    idx_t c1[3] = {m, b, other};
                    fix_tri_orient(npts, sdim, c0, oref);
                    fix_tri_orient(npts, sdim, c1, oref);
                    for (int d = 0; d < 3; ++d) {
                        nelems[d][w] = c0[d];
                    }
                    nparent[w] = parent[e];
                    nforce[w]  = force[e];
                    nlei[w]    = (uint8_t)local_edge_of(let, n_lei, c0, a, other);
                    ++w;
                    for (int d = 0; d < 3; ++d) {
                        nelems[d][w] = c1[d];
                    }
                    nparent[w] = parent[e];
                    nforce[w]  = force[e];
                    nlei[w]    = (uint8_t)local_edge_of(let, n_lei, c1, b, other);
                    ++w;
                }
            }
        }

        for (int d = 0; d < nxe; ++d) {
            SMESH_FREE(elems[d]);
        }
        SMESH_FREE(elems);
        for (int d = 0; d < sdim; ++d) {
            SMESH_FREE(pts[d]);
        }
        SMESH_FREE(pts);
        SMESH_FREE(parent);
        SMESH_FREE(force);
        SMESH_FREE(lei);
        SMESH_FREE(node_a);
        SMESH_FREE(node_b);
        SMESH_FREE(h_work);
        SMESH_FREE(g_work);

        elems   = nelems;
        pts     = npts;
        parent  = nparent;
        force   = nforce;
        lei     = nlei;
        node_a  = nna;
        node_b  = nnb;
        h_work  = h_new;
        g_work  = g_new;
        n_elem  = n_elem_new;
        n_nodes = n_nodes_new;

        SMESH_FREE(rowptr);
        SMESH_FREE(colidx);
        SMESH_FREE(uid);
        SMESH_FREE(eflag);
        SMESH_FREE(emark);
        SMESH_FREE(longest);
        SMESH_FREE(long_uid);
        SMESH_FREE(mid_of);
        SMESH_FREE(split);
    }

    SMESH_FREE(h_work);
    SMESH_FREE(g_work);
    SMESH_FREE(force);
    SMESH_FREE(lei);

    count_t *pp = (count_t *)SMESH_CALLOC((size_t)(n_elements_in + 1), sizeof(count_t));
    for (ptrdiff_t e = 0; e < n_elem; ++e) {
        const ptrdiff_t p = parent[e];
        if (p >= 0 && p < n_elements_in) {
            pp[p + 1] += 1;
        }
    }
    for (ptrdiff_t e = 0; e < n_elements_in; ++e) {
        pp[e + 1] += pp[e];
    }
    idx_t *cid = (idx_t *)SMESH_ALLOC((size_t)n_elem * sizeof(idx_t));
    count_t *fill = (count_t *)SMESH_ALLOC((size_t)(n_elements_in + 1) * sizeof(count_t));
    memcpy(fill, pp, (size_t)(n_elements_in + 1) * sizeof(count_t));
    for (ptrdiff_t e = 0; e < n_elem; ++e) {
        const ptrdiff_t p = parent[e];
        if (p >= 0 && p < n_elements_in) {
            cid[fill[p]++] = (idx_t)e;
        }
    }
    SMESH_FREE(fill);

    *n_elements_out  = n_elem;
    *elements_out    = elems;
    *n_nodes_out     = n_nodes;
    *points_out      = pts;
    *parent_elem_out = parent;
    *parent_ptr_out  = pp;
    *child_id_out    = cid;
    *node_a_out      = node_a;
    *node_b_out      = node_b;
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif
