#include "smesh_node_height.hpp"

#include "smesh_types.hpp"

#include <cmath>
#include <vector>

namespace smesh {

namespace {

enum class Shape { Unsupported, Edge, Tri, Quad, Tet, Hex, Wedge, Pyramid };

Shape shape_of(enum ElemType type) {
    switch (type) {
        case EDGE2:
        case EDGE3:
        case EDGESHELL2:
        case EDGESHELL3:
        case BEAM2: return Shape::Edge;
        case TRI3:
        case TRI6:
        case TRI10:
        case TRISHELL3:
        case TRISHELL6:
        case TRISHELL10: return Shape::Tri;
        case QUAD4:
        case QUAD9:
        case QUADSHELL4:
        case QUADSHELL9: return Shape::Quad;
        case TET4:
        case TET10:
        case TET15:
        case TET20: return Shape::Tet;
        case HEX8:
        case HEX27: return Shape::Hex;
        case WEDGE6: return Shape::Wedge;
        case PYRAMID5: return Shape::Pyramid;
        default: return Shape::Unsupported;
    }
}

void sub(const double *a, const double *b, double *r) {
    for (int d = 0; d < 3; d++) r[d] = a[d] - b[d];
}
void cross(const double *u, const double *v, double *r) {
    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];
}
double norm(const double *v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

/// 3 V / A_max of the tet (a, b, c, d).
double tet_height(const double *a, const double *b, const double *c, const double *d) {
    double u[3], v[3], w[3], n[3];
    sub(b, a, u);
    sub(c, a, v);
    sub(d, a, w);
    cross(v, w, n);
    const double vol6 = std::fabs(u[0] * n[0] + u[1] * n[1] + u[2] * n[2]);
    const double *f[4][3] = {{b, c, d}, {a, c, d}, {a, b, d}, {a, b, c}};
    double        amax2   = 0;  // (2 A)^2
    for (auto &t : f) {
        double e1[3], e2[3], cr[3];
        sub(t[1], t[0], e1);
        sub(t[2], t[0], e2);
        cross(e1, e2, cr);
        const double a2 = cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2];
        amax2           = a2 > amax2 ? a2 : amax2;
    }
    return amax2 > 0 ? vol6 / std::sqrt(amax2) : 0.0;  // 3V / A = 6V / 2A
}

/// 2 A / longest edge of the triangle (a, b, c), in 3D.
double tri_height(const double *a, const double *b, const double *c) {
    double e1[3], e2[3], e3[3], cr[3];
    sub(b, a, e1);
    sub(c, a, e2);
    sub(c, b, e3);
    cross(e1, e2, cr);
    const double area2 = norm(cr);  // 2 A
    double       lmax  = norm(e1);
    lmax               = norm(e2) > lmax ? norm(e2) : lmax;
    lmax               = norm(e3) > lmax ? norm(e3) : lmax;
    return lmax > 0 ? area2 / lmax : 0.0;
}

// HEX8: corner c and its three edge neighbours (standard ordering, 0-3 bottom, 4-7 top).
const int kHexCorner[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5}, {2, 3, 1, 6}, {3, 0, 2, 7},
                              {4, 5, 7, 0}, {5, 6, 4, 1}, {6, 7, 5, 2}, {7, 4, 6, 3}};
// WEDGE6: 0-2 bottom, 3-5 top.
const int kWedgeCorner[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4}, {2, 0, 1, 5}, {3, 5, 4, 0}, {4, 3, 5, 1}, {5, 4, 3, 2}};
// PYRAMID5: base 0-3, apex 4.
const int kPyramidCorner[4][4] = {{0, 1, 3, 4}, {1, 2, 0, 4}, {2, 3, 1, 4}, {3, 0, 2, 4}};

}  // namespace

bool node_height_supported(enum ElemType type) { return shape_of(type) != Shape::Unsupported; }

double element_min_height(enum ElemType type, const double (*p)[3]) {
    double h = INFINITY;
    switch (shape_of(type)) {
        case Shape::Edge: {
            double e[3];
            sub(p[1], p[0], e);
            return norm(e);
        }
        case Shape::Tri: return tri_height(p[0], p[1], p[2]);
        case Shape::Quad:
            for (int c = 0; c < 4; c++) {
                const double t = tri_height(p[c], p[(c + 1) % 4], p[(c + 3) % 4]);
                h              = t < h ? t : h;
            }
            return h;
        case Shape::Tet: return tet_height(p[0], p[1], p[2], p[3]);
        case Shape::Hex:
            for (const auto &c : kHexCorner) {
                const double t = tet_height(p[c[0]], p[c[1]], p[c[2]], p[c[3]]);
                h              = t < h ? t : h;
            }
            return h;
        case Shape::Wedge:
            for (const auto &c : kWedgeCorner) {
                const double t = tet_height(p[c[0]], p[c[1]], p[c[2]], p[c[3]]);
                h              = t < h ? t : h;
            }
            return h;
        case Shape::Pyramid:
            for (const auto &c : kPyramidCorner) {
                const double t = tet_height(p[c[0]], p[c[1]], p[c[2]], p[c[3]]);
                h              = t < h ? t : h;
            }
            return h;
        case Shape::Unsupported: break;
    }
    return NAN;
}

template <typename idx_t, typename geom_t>
int accumulate_node_heights(enum ElemType type, ptrdiff_t n_elements, const idx_t *const *elements, int spatial_dim,
                            const geom_t *const *points, double *h) {
    if (!node_height_supported(type)) {
        return SMESH_FAILURE;
    }
    const int           nxe = elem_num_nodes(type);
    std::vector<double> he(static_cast<size_t>(n_elements));
#pragma omp parallel for schedule(static)
    for (ptrdiff_t e = 0; e < n_elements; e++) {
        double p[27][3] = {};
        const int ncorner = nxe < 27 ? nxe : 27;
        for (int k = 0; k < ncorner; k++) {
            const idx_t n = elements[k][e];
            for (int d = 0; d < 3; d++) p[k][d] = d < spatial_dim ? double(points[d][n]) : 0.0;
        }
        he[size_t(e)] = element_min_height(type, p);
    }
    for (ptrdiff_t e = 0; e < n_elements; e++) {
        const double v = he[size_t(e)];
        for (int k = 0; k < nxe; k++) {
            const idx_t n = elements[k][e];
            h[n]          = v < h[n] ? v : h[n];
        }
    }
    return SMESH_SUCCESS;
}

template int accumulate_node_heights<i32, f32>(enum ElemType, ptrdiff_t, const i32 *const *, int, const f32 *const *,
                                               double *);
template int accumulate_node_heights<i32, f64>(enum ElemType, ptrdiff_t, const i32 *const *, int, const f64 *const *,
                                               double *);
template int accumulate_node_heights<i64, f32>(enum ElemType, ptrdiff_t, const i64 *const *, int, const f32 *const *,
                                               double *);
template int accumulate_node_heights<i64, f64>(enum ElemType, ptrdiff_t, const i64 *const *, int, const f64 *const *,
                                               double *);

}  // namespace smesh
