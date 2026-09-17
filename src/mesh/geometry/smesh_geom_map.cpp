#include "smesh_geom_map.hpp"

#include "smesh_sshex8.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sswedge.hpp"

namespace smesh {

struct GeomMapFlags {
    int affine;
    int axis_aligned;
};

static SMESH_INLINE geom_t abs_g(const geom_t x) { return x < geom_t(0) ? -x : x; }

static SMESH_INLINE geom_t max_g(const geom_t a, const geom_t b) { return a > b ? a : b; }

static SMESH_INLINE void gather_node(const int                                               sdim,
                                     const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                     const idx_t                                             i,
                                     geom_t                                                  p[3]) {
    p[0] = points[0][i];
    p[1] = sdim > 1 ? points[1][i] : geom_t(0);
    p[2] = sdim > 2 ? points[2][i] : geom_t(0);
}

static SMESH_INLINE void sub3(const geom_t a[3], const geom_t b[3], geom_t o[3]) {
    o[0] = a[0] - b[0];
    o[1] = a[1] - b[1];
    o[2] = a[2] - b[2];
}

static SMESH_INLINE geom_t inf_norm3(const geom_t v[3]) {
    return max_g(abs_g(v[0]), max_g(abs_g(v[1]), abs_g(v[2])));
}

static SMESH_INLINE int close3(const geom_t a[3], const geom_t b[3], const geom_t tol) {
    const int c0 = abs_g(a[0] - b[0]) <= tol;
    const int c1 = abs_g(a[1] - b[1]) <= tol;
    const int c2 = abs_g(a[2] - b[2]) <= tol;
    return c0 & c1 & c2;
}

static SMESH_INLINE int axis_bit(const geom_t e[3], const geom_t tol) {
    const geom_t ax = abs_g(e[0]);
    const geom_t ay = abs_g(e[1]);
    const geom_t az = abs_g(e[2]);
    const int    x_ok = int(ax > tol) & int(ay <= tol) & int(az <= tol);
    const int    y_ok = int(ay > tol) & int(ax <= tol) & int(az <= tol);
    const int    z_ok = int(az > tol) & int(ax <= tol) & int(ay <= tol);
    return x_ok | (y_ok << 1) | (z_ok << 2);
}

static SMESH_INLINE geom_t edge_scale(const geom_t e1[3], const geom_t e2[3], const geom_t e3[3]) {
    const geom_t s = max_g(inf_norm3(e1), max_g(inf_norm3(e2), inf_norm3(e3)));
    return s > geom_t(0) ? s : geom_t(1);
}

static SMESH_INLINE int point_matches_any(const geom_t p[3], const int nexp, const geom_t *exp, const geom_t tol) {
    int ok = 0;
    for (int k = 0; k < nexp; ++k) {
        ok |= close3(p, exp + 3 * k, tol);
    }
    return ok;
}

static SMESH_INLINE enum GeomMap map_from_flags(const GeomMapFlags f, const enum ElemType type) {
    if (!f.affine) {
        return ISOPARAMETRIC;
    }
    if (f.axis_aligned && geom_map_allows_axis_aligned(type)) {
        return AXIS_ALIGNED;
    }
    return AFFINE;
}

/// VTK HEX8 from (x,y,z) in {0,1}: binary x+2y+4z with 2↔3 and 6↔7 swapped.
static SMESH_INLINE int hex8_vtk_slot(const int x, const int y, const int z) {
    static const int vtk[8] = {0, 1, 3, 2, 4, 5, 7, 6};
    return vtk[x + 2 * y + 4 * z];
}

/// Inverse of create_cube HEX27 cartesian permutation: slot s holds lattice index
/// hex27_to_cartesian[s] = x + 3 y + 9 z.
static SMESH_INLINE int hex27_slot(const int x, const int y, const int z) {
    static const int cartesian_to_hex27[27] = {0,  8,  1,  11, 24, 9,  3,  10, 2,  16, 20, 17, 23, 26,
                                               21, 19, 22, 18, 4,  12, 5,  15, 25, 13, 7,  14, 6};
    return cartesian_to_hex27[x + 3 * y + 9 * z];
}

static SMESH_INLINE int hex_slot(const enum ElemType type, const int L, const int x, const int y, const int z) {
    if (is_semistructured_type(type)) {
        return sshex8_lidx(L, x, y, z);
    }
    if (type == HEX27) {
        return hex27_slot(x, y, z);
    }
    return hex8_vtk_slot(x, y, z);
}

static SMESH_INLINE int quad4_vtk_slot(const int x, const int y) {
    static const int vtk[4] = {0, 1, 3, 2};
    return vtk[x + 2 * y];
}

static SMESH_INLINE int quad9_slot(const int x, const int y) {
    static const int t[9] = {0, 4, 1, 7, 8, 5, 3, 6, 2};
    return t[y * 3 + x];
}

static SMESH_INLINE int quad_slot(const enum ElemType type, const int L, const int x, const int y) {
    if (is_semistructured_type(type)) {
        return ssquad4_lidx(L, x, y);
    }
    if (type == QUAD9 || type == QUADSHELL9) {
        return quad9_slot(x, y);
    }
    return quad4_vtk_slot(x, y);
}

static GeomMapFlags detect_hex_lattice(const enum ElemType                                     type,
                                       const int                                               L,
                                       const ptrdiff_t                                         ne,
                                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT els,
                                       const int                                               sdim,
                                       const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                       const geom_t                                            rel_tol) {
    GeomMapFlags flags = {1, 1};
    const geom_t invL  = geom_t(1) / geom_t(L);
    const int    i0    = hex_slot(type, L, 0, 0, 0);
    const int    i1    = hex_slot(type, L, L, 0, 0);
    const int    i3    = hex_slot(type, L, 0, L, 0);
    const int    i4    = hex_slot(type, L, 0, 0, L);

    for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t v0[3], v1[3], v3[3], v4[3], e1[3], e2[3], e3[3];
        gather_node(sdim, points, els[i0][e], v0);
        gather_node(sdim, points, els[i1][e], v1);
        gather_node(sdim, points, els[i3][e], v3);
        gather_node(sdim, points, els[i4][e], v4);
        sub3(v1, v0, e1);
        sub3(v3, v0, e2);
        sub3(v4, v0, e3);
        const geom_t tol = rel_tol * edge_scale(e1, e2, e3);
        const int    b1  = axis_bit(e1, tol);
        const int    b2  = axis_bit(e2, tol);
        const int    b3  = axis_bit(e3, tol);
        flags.axis_aligned &= int(b1 != 0) & int(b2 != 0) & int(b3 != 0) & int(b1 != b2) & int(b1 != b3) & int(b2 != b3);

        int affine = 1;
        for (int zi = 0; zi <= L; ++zi) {
            for (int yi = 0; yi <= L; ++yi) {
                for (int xi = 0; xi <= L; ++xi) {
                    geom_t pred[3], got[3];
                    pred[0] = v0[0] + geom_t(xi) * invL * e1[0] + geom_t(yi) * invL * e2[0] + geom_t(zi) * invL * e3[0];
                    pred[1] = v0[1] + geom_t(xi) * invL * e1[1] + geom_t(yi) * invL * e2[1] + geom_t(zi) * invL * e3[1];
                    pred[2] = v0[2] + geom_t(xi) * invL * e1[2] + geom_t(yi) * invL * e2[2] + geom_t(zi) * invL * e3[2];
                    gather_node(sdim, points, els[hex_slot(type, L, xi, yi, zi)][e], got);
                    affine &= close3(pred, got, tol);
                }
            }
        }
        flags.affine &= affine;
        flags.axis_aligned &= affine;
    }
    return flags;
}

static GeomMapFlags detect_quad_lattice(const enum ElemType                                     type,
                                        const int                                               L,
                                        const ptrdiff_t                                         ne,
                                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT els,
                                        const int                                               sdim,
                                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                        const geom_t                                            rel_tol) {
    GeomMapFlags flags = {1, 1};
    const geom_t invL  = geom_t(1) / geom_t(L);
    const int    i0    = quad_slot(type, L, 0, 0);
    const int    i1    = quad_slot(type, L, L, 0);
    const int    i3    = quad_slot(type, L, 0, L);
    geom_t       e3z[3] = {0, 0, 0};

    for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t v0[3], v1[3], v3[3], e1[3], e2[3];
        gather_node(sdim, points, els[i0][e], v0);
        gather_node(sdim, points, els[i1][e], v1);
        gather_node(sdim, points, els[i3][e], v3);
        sub3(v1, v0, e1);
        sub3(v3, v0, e2);
        const geom_t tol = rel_tol * edge_scale(e1, e2, e3z);
        const int    b1  = axis_bit(e1, tol);
        const int    b2  = axis_bit(e2, tol);
        flags.axis_aligned &= int(b1 != 0) & int(b2 != 0) & int(b1 != b2);

        int affine = 1;
        for (int yi = 0; yi <= L; ++yi) {
            for (int xi = 0; xi <= L; ++xi) {
                geom_t pred[3], got[3];
                pred[0] = v0[0] + geom_t(xi) * invL * e1[0] + geom_t(yi) * invL * e2[0];
                pred[1] = v0[1] + geom_t(xi) * invL * e1[1] + geom_t(yi) * invL * e2[1];
                pred[2] = v0[2] + geom_t(xi) * invL * e1[2] + geom_t(yi) * invL * e2[2];
                gather_node(sdim, points, els[quad_slot(type, L, xi, yi)][e], got);
                affine &= close3(pred, got, tol);
            }
        }
        flags.affine &= affine;
        flags.axis_aligned &= affine;
    }
    return flags;
}

static GeomMapFlags detect_tet_family(const enum ElemType                                     type,
                                      const ptrdiff_t                                         ne,
                                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT els,
                                      const int                                               sdim,
                                      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                      const geom_t                                            rel_tol) {
    GeomMapFlags flags = {1, 0};
    const int    ss    = is_semistructured_type(type);
    const int    L     = ss ? semistructured_level(type) : 1;
    const int    nxe   = elem_num_nodes(type);

    for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t v0[3], v1[3], v2[3], v3[3], e1[3], e2[3], e3[3];
        const idx_t a0 = ss ? els[sstet4_lidx(L, 0, 0, 0)][e] : els[0][e];
        const idx_t a1 = ss ? els[sstet4_lidx(L, L, 0, 0)][e] : els[1][e];
        const idx_t a2 = ss ? els[sstet4_lidx(L, 0, L, 0)][e] : els[2][e];
        const idx_t a3 = ss ? els[sstet4_lidx(L, 0, 0, L)][e] : els[3][e];
        gather_node(sdim, points, a0, v0);
        gather_node(sdim, points, a1, v1);
        gather_node(sdim, points, a2, v2);
        gather_node(sdim, points, a3, v3);
        sub3(v1, v0, e1);
        sub3(v2, v0, e2);
        sub3(v3, v0, e3);
        const geom_t tol = rel_tol * edge_scale(e1, e2, e3);

        if (ss) {
            const geom_t invL = geom_t(1) / geom_t(L);
            int          affine = 1;
            for (int z = 0; z <= L; ++z) {
                for (int y = 0; y <= L - z; ++y) {
                    for (int x = 0; x <= L - z - y; ++x) {
                        geom_t pred[3], got[3];
                        pred[0] = v0[0] + geom_t(x) * invL * e1[0] + geom_t(y) * invL * e2[0] + geom_t(z) * invL * e3[0];
                        pred[1] = v0[1] + geom_t(x) * invL * e1[1] + geom_t(y) * invL * e2[1] + geom_t(z) * invL * e3[1];
                        pred[2] = v0[2] + geom_t(x) * invL * e1[2] + geom_t(y) * invL * e2[2] + geom_t(z) * invL * e3[2];
                        gather_node(sdim, points, els[sstet4_lidx(L, x, y, z)][e], got);
                        affine &= close3(pred, got, tol);
                    }
                }
            }
            flags.affine &= affine;
            continue;
        }

        if (nxe <= 4) {
            continue;
        }

        geom_t       exp[11 * 3];
        const geom_t *va[6] = {v0, v1, v0, v0, v1, v2};
        const geom_t *vb[6] = {v1, v2, v2, v3, v3, v3};
        for (int k = 0; k < 6; ++k) {
            exp[3 * k + 0] = geom_t(0.5) * (va[k][0] + vb[k][0]);
            exp[3 * k + 1] = geom_t(0.5) * (va[k][1] + vb[k][1]);
            exp[3 * k + 2] = geom_t(0.5) * (va[k][2] + vb[k][2]);
        }
        int nexp = 6;
        if (nxe >= 15) {
            const geom_t *vs[4]       = {v0, v1, v2, v3};
            const int     faces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
            for (int f = 0; f < 4; ++f) {
                exp[3 * nexp + 0] = (vs[faces[f][0]][0] + vs[faces[f][1]][0] + vs[faces[f][2]][0]) / geom_t(3);
                exp[3 * nexp + 1] = (vs[faces[f][0]][1] + vs[faces[f][1]][1] + vs[faces[f][2]][1]) / geom_t(3);
                exp[3 * nexp + 2] = (vs[faces[f][0]][2] + vs[faces[f][1]][2] + vs[faces[f][2]][2]) / geom_t(3);
                ++nexp;
            }
            exp[3 * nexp + 0] = (v0[0] + v1[0] + v2[0] + v3[0]) / geom_t(4);
            exp[3 * nexp + 1] = (v0[1] + v1[1] + v2[1] + v3[1]) / geom_t(4);
            exp[3 * nexp + 2] = (v0[2] + v1[2] + v2[2] + v3[2]) / geom_t(4);
            ++nexp;
        }

        int affine = 1;
        for (int n = 4; n < nxe; ++n) {
            geom_t p[3];
            gather_node(sdim, points, els[n][e], p);
            affine &= point_matches_any(p, nexp, exp, tol);
        }
        flags.affine &= affine;
    }
    return flags;
}

static GeomMapFlags detect_tri_family(const enum ElemType                                     type,
                                      const ptrdiff_t                                         ne,
                                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT els,
                                      const int                                               sdim,
                                      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                      const geom_t                                            rel_tol) {
    GeomMapFlags flags = {1, 0};
    const int    nxe   = elem_num_nodes(type);
    if (nxe <= 3) {
        return flags;
    }

    geom_t e3z[3] = {0, 0, 0};
    for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t v0[3], v1[3], v2[3], e1[3], e2[3];
        gather_node(sdim, points, els[0][e], v0);
        gather_node(sdim, points, els[1][e], v1);
        gather_node(sdim, points, els[2][e], v2);
        sub3(v1, v0, e1);
        sub3(v2, v0, e2);
        const geom_t tol = rel_tol * edge_scale(e1, e2, e3z);
        geom_t       exp[3 * 3];
        exp[0] = geom_t(0.5) * (v0[0] + v1[0]);
        exp[1] = geom_t(0.5) * (v0[1] + v1[1]);
        exp[2] = geom_t(0.5) * (v0[2] + v1[2]);
        exp[3] = geom_t(0.5) * (v1[0] + v2[0]);
        exp[4] = geom_t(0.5) * (v1[1] + v2[1]);
        exp[5] = geom_t(0.5) * (v1[2] + v2[2]);
        exp[6] = geom_t(0.5) * (v0[0] + v2[0]);
        exp[7] = geom_t(0.5) * (v0[1] + v2[1]);
        exp[8] = geom_t(0.5) * (v0[2] + v2[2]);
        int affine = 1;
        for (int n = 3; n < nxe && n < 6; ++n) {
            geom_t p[3];
            gather_node(sdim, points, els[n][e], p);
            affine &= point_matches_any(p, 3, exp, tol);
        }
        flags.affine &= affine;
    }
    return flags;
}

static GeomMapFlags detect_wedge_family(const enum ElemType                                     type,
                                        const ptrdiff_t                                         ne,
                                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT els,
                                        const int                                               sdim,
                                        const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                        const geom_t                                            rel_tol) {
    GeomMapFlags flags = {1, 0};
    const int    ss    = is_semistructured_type(type);
    const int    L     = ss ? semistructured_level(type) : 1;
    const geom_t invL  = geom_t(1) / geom_t(L);

    for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t v0[3], v1[3], v2[3], v3[3], e1[3], e2[3], e3[3];
        const idx_t a0 = ss ? els[sswedge_lidx(L, 0, 0, 0)][e] : els[0][e];
        const idx_t a1 = ss ? els[sswedge_lidx(L, L, 0, 0)][e] : els[1][e];
        const idx_t a2 = ss ? els[sswedge_lidx(L, 0, L, 0)][e] : els[2][e];
        const idx_t a3 = ss ? els[sswedge_lidx(L, 0, 0, L)][e] : els[3][e];
        gather_node(sdim, points, a0, v0);
        gather_node(sdim, points, a1, v1);
        gather_node(sdim, points, a2, v2);
        gather_node(sdim, points, a3, v3);
        sub3(v1, v0, e1);
        sub3(v2, v0, e2);
        sub3(v3, v0, e3);
        const geom_t tol = rel_tol * edge_scale(e1, e2, e3);

        if (!ss) {
            geom_t v4[3], v5[3], p4[3], p5[3];
            gather_node(sdim, points, els[4][e], v4);
            gather_node(sdim, points, els[5][e], v5);
            p4[0] = v0[0] + e1[0] + e3[0];
            p4[1] = v0[1] + e1[1] + e3[1];
            p4[2] = v0[2] + e1[2] + e3[2];
            p5[0] = v0[0] + e2[0] + e3[0];
            p5[1] = v0[1] + e2[1] + e3[1];
            p5[2] = v0[2] + e2[2] + e3[2];
            flags.affine &= close3(v4, p4, tol) & close3(v5, p5, tol);
            continue;
        }

        int affine = 1;
        for (int z = 0; z <= L; ++z) {
            for (int y = 0; y <= L; ++y) {
                for (int x = 0; x <= L - y; ++x) {
                    geom_t pred[3], got[3];
                    pred[0] = v0[0] + geom_t(x) * invL * e1[0] + geom_t(y) * invL * e2[0] + geom_t(z) * invL * e3[0];
                    pred[1] = v0[1] + geom_t(x) * invL * e1[1] + geom_t(y) * invL * e2[1] + geom_t(z) * invL * e3[1];
                    pred[2] = v0[2] + geom_t(x) * invL * e1[2] + geom_t(y) * invL * e2[2] + geom_t(z) * invL * e3[2];
                    gather_node(sdim, points, els[sswedge_lidx(L, x, y, z)][e], got);
                    affine &= close3(pred, got, tol);
                }
            }
        }
        flags.affine &= affine;
    }
    return flags;
}

static GeomMapFlags detect_edge_family(const enum ElemType                                     type,
                                       const ptrdiff_t                                         ne,
                                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT els,
                                       const int                                               sdim,
                                       const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                       const geom_t                                            rel_tol) {
    GeomMapFlags flags = {1, 0};
    const int    nxe   = elem_num_nodes(type);
    if (nxe <= 2) {
        return flags;
    }
    geom_t z[3] = {0, 0, 0};
    for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t v0[3], v1[3], e1[3], mid[3], p[3];
        gather_node(sdim, points, els[0][e], v0);
        gather_node(sdim, points, els[1][e], v1);
        sub3(v1, v0, e1);
        const geom_t tol = rel_tol * edge_scale(e1, z, z);
        mid[0]           = geom_t(0.5) * (v0[0] + v1[0]);
        mid[1]           = geom_t(0.5) * (v0[1] + v1[1]);
        mid[2]           = geom_t(0.5) * (v0[2] + v1[2]);
        gather_node(sdim, points, els[2][e], p);
        flags.affine &= close3(p, mid, tol);
    }
    return flags;
}

enum GeomMap detect_geom_map(const enum ElemType                                     type,
                             const ptrdiff_t                                         nelements,
                             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                             const int                                               sdim,
                             const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                             const geom_t                                            rel_tol) {
    if (nelements <= 0 || !elements || !points || sdim < 1) {
        return ISOPARAMETRIC;
    }

    const geom_t tol = rel_tol > geom_t(0) ? rel_tol : geom_map_default_rel_tol();

    if (type == HEX8 || type == HEX27 || is_hex_ss_family(type)) {
        const int L = is_semistructured_type(type) ? semistructured_level(type) : (type == HEX27 ? 2 : 1);
        return map_from_flags(detect_hex_lattice(type, L, nelements, elements, sdim, points, tol), type);
    }
    if (type == QUAD9 || type == QUADSHELL9 || is_quad_ss_family(type)) {
        const int L = is_semistructured_type(type) ? semistructured_level(type)
                                                  : ((type == QUAD9 || type == QUADSHELL9) ? 2 : 1);
        return map_from_flags(detect_quad_lattice(type, L, nelements, elements, sdim, points, tol), type);
    }
    if (type == TET10 || type == TET15 || type == TET20 || is_tet_ss_family(type)) {
        return map_from_flags(detect_tet_family(type, nelements, elements, sdim, points, tol), type);
    }
    if (type == TRI3 || type == TRI6 || type == TRI10 || type == TRISHELL3 || type == TRISHELL6 || type == MACRO_TRI3 ||
        type == MACRO_TRISHELL3) {
        return map_from_flags(detect_tri_family(type, nelements, elements, sdim, points, tol), type);
    }
    if (is_wedge_ss_family(type)) {
        return map_from_flags(detect_wedge_family(type, nelements, elements, sdim, points, tol), type);
    }
    if (type == EDGE2 || type == EDGE3 || type == EDGESHELL2 || type == EDGESHELL3 || type == BEAM2) {
        return map_from_flags(detect_edge_family(type, nelements, elements, sdim, points, tol), type);
    }
    if (type == NODE1) {
        return AFFINE;
    }
    return ISOPARAMETRIC;
}

}  // namespace smesh
