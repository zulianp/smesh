#include "smesh_mesh.hpp"

#include "smesh_buffer.hpp"
#include "smesh_tracer.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace smesh {
namespace {

ptrdiff_t cylinder_n_core(const ptrdiff_t ntheta) {
    const ptrdiff_t n = ntheta / 5;
    return n < 1 ? 1 : n;
}

geom_t quad_signed_area(const geom_t *const *const pts, const idx_t i0, const idx_t i1,
                        const idx_t i2, const idx_t i3) {
    const geom_t x[4] = {pts[0][i0], pts[0][i1], pts[0][i2], pts[0][i3]};
    const geom_t y[4] = {pts[1][i0], pts[1][i1], pts[1][i2], pts[1][i3]};
    geom_t       a    = 0;
    for (int k = 0; k < 4; ++k) {
        const int n = (k + 1) & 3;
        a += x[k] * y[n] - x[n] * y[k];
    }
    return a;
}

geom_t tri_signed_area(const geom_t *const *const pts, const idx_t i0, const idx_t i1,
                       const idx_t i2) {
    const geom_t x0 = pts[0][i0];
    const geom_t y0 = pts[1][i0];
    return (pts[0][i1] - x0) * (pts[1][i2] - y0) - (pts[0][i2] - x0) * (pts[1][i1] - y0);
}

void emit_hex_prism(idx_t *const *const hex, const ptrdiff_t e, idx_t b0, idx_t b1, idx_t b2,
                    idx_t b3, const ptrdiff_t n_plane, const geom_t *const *const pts) {
    if (quad_signed_area(pts, b0, b1, b2, b3) < 0) {
        const idx_t t = b1;
        b1            = b3;
        b3            = t;
    }
    hex[0][e] = b0;
    hex[1][e] = b1;
    hex[2][e] = b2;
    hex[3][e] = b3;
    hex[4][e] = static_cast<idx_t>(b0 + n_plane);
    hex[5][e] = static_cast<idx_t>(b1 + n_plane);
    hex[6][e] = static_cast<idx_t>(b2 + n_plane);
    hex[7][e] = static_cast<idx_t>(b3 + n_plane);
}

void emit_wedge_prism(idx_t *const *const wedge, const ptrdiff_t e, idx_t b0, idx_t b1, idx_t b2,
                      const ptrdiff_t n_plane, const geom_t *const *const pts) {
    if (tri_signed_area(pts, b0, b1, b2) < 0) {
        const idx_t t = b1;
        b1            = b2;
        b2            = t;
    }
    wedge[0][e] = b0;
    wedge[1][e] = b1;
    wedge[2][e] = b2;
    wedge[3][e] = static_cast<idx_t>(b0 + n_plane);
    wedge[4][e] = static_cast<idx_t>(b1 + n_plane);
    wedge[5][e] = static_cast<idx_t>(b2 + n_plane);
}

idx_t core_node(const ptrdiff_t xi, const ptrdiff_t yi, const ptrdiff_t k, const ptrdiff_t n_core,
                const ptrdiff_t n_plane) {
    return static_cast<idx_t>(k * n_plane + yi * (n_core + 1) + xi);
}

idx_t polar_node(const ptrdiff_t l, const ptrdiff_t i, const ptrdiff_t k, const ptrdiff_t ntheta,
                 const ptrdiff_t n_core_plane, const ptrdiff_t n_plane) {
    return static_cast<idx_t>(k * n_plane + n_core_plane + l * ntheta + i);
}

idx_t square_side_node(const int side, const ptrdiff_t t, const ptrdiff_t k, const ptrdiff_t n_core,
                       const ptrdiff_t n_plane) {
    ptrdiff_t xi = 0;
    ptrdiff_t yi = 0;
    if (side == 0) {
        xi = n_core;
        yi = t;
    } else if (side == 1) {
        xi = n_core - t;
        yi = n_core;
    } else if (side == 2) {
        xi = 0;
        yi = n_core - t;
    } else {
        xi = t;
        yi = 0;
    }
    return core_node(xi, yi, k, n_core, n_plane);
}

}  // namespace

std::shared_ptr<Mesh> Mesh::create_hex_dominant_cylinder(const std::shared_ptr<Communicator> &comm,
                                                         const geom_t                         radius,
                                                         const geom_t                         height,
                                                         const ptrdiff_t                      nr,
                                                         const ptrdiff_t                      ntheta,
                                                         const ptrdiff_t                      nz,
                                                         const geom_t                         zmin) {
    SMESH_TRACE_SCOPE("Mesh::create_hex_dominant_cylinder");

    if (radius <= 0 || height <= 0 || nr < 1 || nz < 1 || ntheta < 8 || (ntheta % 4) != 0) {
        SMESH_ERROR(
                "create_hex_dominant_cylinder: require radius>0, height>0, nr>=1, nz>=1, "
                "ntheta>=8 and ntheta %% 4 == 0 (got radius=%g height=%g nr=%ld ntheta=%ld nz=%ld)\n",
                (double)radius,
                (double)height,
                (long)nr,
                (long)ntheta,
                (long)nz);
        return nullptr;
    }

    const ptrdiff_t n_core = cylinder_n_core(ntheta);
    const ptrdiff_t n_inner_side = n_core + 1;
    const ptrdiff_t n_arc_side   = ntheta / 4 + 1;
    const ptrdiff_t n_quad_side  = (n_inner_side < n_arc_side ? n_inner_side : n_arc_side) - 1;
    const ptrdiff_t n_tri_side =
            n_arc_side > n_inner_side ? n_arc_side - n_inner_side : n_inner_side - n_arc_side;

    const ptrdiff_t n_hex_core       = n_core * n_core * nz;
    const ptrdiff_t n_hex_annulus    = nr * ntheta * nz;
    const ptrdiff_t n_hex_transition = 4 * n_quad_side * nz;
    const ptrdiff_t n_hex            = n_hex_core + n_hex_annulus + n_hex_transition;
    const ptrdiff_t n_wedge          = 4 * n_tri_side * nz;

    const ptrdiff_t n_core_plane  = (n_core + 1) * (n_core + 1);
    const ptrdiff_t n_polar_plane = (nr + 1) * ntheta;
    const ptrdiff_t n_plane       = n_core_plane + n_polar_plane;
    const ptrdiff_t n_nodes       = n_plane * (nz + 1);

    const double R = static_cast<double>(radius);
    double r_in    = R / (1.0 + 2.0 * M_PI * static_cast<double>(nr) / static_cast<double>(ntheta));
    if (r_in < 0.3 * R) {
        r_in = 0.3 * R;
    }
    if (r_in > 0.65 * R) {
        r_in = 0.65 * R;
    }

    double a = static_cast<double>(n_core) * M_PI * r_in / static_cast<double>(ntheta);
    const double a_max = 0.9 * r_in / std::sqrt(2.0);
    if (a > a_max) {
        a = a_max;
    }

    auto points_buf = create_host_buffer<geom_t>(3, static_cast<size_t>(n_nodes));
    auto hex_buf    = create_host_buffer<idx_t>(8, static_cast<size_t>(n_hex));
    auto wedge_buf =
            n_wedge > 0 ? create_host_buffer<idx_t>(6, static_cast<size_t>(n_wedge)) : nullptr;

    auto pts = points_buf->data();
    auto hex = hex_buf->data();
    idx_t *const *wedge = wedge_buf ? wedge_buf->data() : nullptr;

    const double inv_n_core = 1.0 / static_cast<double>(n_core);
    const double dtheta     = 2.0 * M_PI / static_cast<double>(ntheta);
    const double dr         = (R - r_in) / static_cast<double>(nr);
    const double dz         = static_cast<double>(height) / static_cast<double>(nz);
    const double z0         = static_cast<double>(zmin);

    for (ptrdiff_t k = 0; k <= nz; ++k) {
        const geom_t z = static_cast<geom_t>(z0 + dz * static_cast<double>(k));
        for (ptrdiff_t yi = 0; yi <= n_core; ++yi) {
            const geom_t y = static_cast<geom_t>(-a + 2.0 * a * static_cast<double>(yi) * inv_n_core);
            for (ptrdiff_t xi = 0; xi <= n_core; ++xi) {
                const idx_t n = core_node(xi, yi, k, n_core, n_plane);
                pts[0][n]     = static_cast<geom_t>(-a + 2.0 * a * static_cast<double>(xi) * inv_n_core);
                pts[1][n]     = y;
                pts[2][n]     = z;
            }
        }
        for (ptrdiff_t l = 0; l <= nr; ++l) {
            const double r = r_in + dr * static_cast<double>(l);
            for (ptrdiff_t i = 0; i < ntheta; ++i) {
                const double ang = -M_PI / 4.0 + dtheta * static_cast<double>(i);
                const idx_t  n   = polar_node(l, i, k, ntheta, n_core_plane, n_plane);
                pts[0][n]        = static_cast<geom_t>(r * std::cos(ang));
                pts[1][n]        = static_cast<geom_t>(r * std::sin(ang));
                pts[2][n]        = z;
            }
        }
    }

    ptrdiff_t eh = 0;
    for (ptrdiff_t k = 0; k < nz; ++k) {
        for (ptrdiff_t yi = 0; yi < n_core; ++yi) {
            for (ptrdiff_t xi = 0; xi < n_core; ++xi) {
                const idx_t i0 = core_node(xi, yi, k, n_core, n_plane);
                const idx_t i1 = core_node(xi + 1, yi, k, n_core, n_plane);
                const idx_t i2 = core_node(xi + 1, yi + 1, k, n_core, n_plane);
                const idx_t i3 = core_node(xi, yi + 1, k, n_core, n_plane);
                emit_hex_prism(hex, eh++, i0, i1, i2, i3, n_plane, pts);
            }
        }
    }

    for (ptrdiff_t k = 0; k < nz; ++k) {
        for (ptrdiff_t l = 0; l < nr; ++l) {
            for (ptrdiff_t i = 0; i < ntheta; ++i) {
                const ptrdiff_t i1 = (i + 1) % ntheta;
                const idx_t     b0 = polar_node(l, i1, k, ntheta, n_core_plane, n_plane);
                const idx_t     b1 = polar_node(l, i, k, ntheta, n_core_plane, n_plane);
                const idx_t     b2 = polar_node(l + 1, i, k, ntheta, n_core_plane, n_plane);
                const idx_t     b3 = polar_node(l + 1, i1, k, ntheta, n_core_plane, n_plane);
                emit_hex_prism(hex, eh++, b0, b1, b2, b3, n_plane, pts);
            }
        }
    }

    ptrdiff_t ew = 0;
    const ptrdiff_t nq = ntheta / 4;
    for (ptrdiff_t k = 0; k < nz; ++k) {
        for (int side = 0; side < 4; ++side) {
            ptrdiff_t si = 0;
            ptrdiff_t sj = 0;
            while (si < n_inner_side - 1 || sj < n_arc_side - 1) {
                const ptrdiff_t rem_i = (n_inner_side - 1) - si;
                const ptrdiff_t rem_j = (n_arc_side - 1) - sj;
                const idx_t inner0 =
                        square_side_node(side, si, k, n_core, n_plane);
                const ptrdiff_t p0 = (static_cast<ptrdiff_t>(side) * nq + sj) % ntheta;
                const idx_t     outer0 =
                        polar_node(0, p0, k, ntheta, n_core_plane, n_plane);

                if (si == n_inner_side - 1) {
                    const ptrdiff_t p1 = (static_cast<ptrdiff_t>(side) * nq + sj + 1) % ntheta;
                    const idx_t     outer1 =
                            polar_node(0, p1, k, ntheta, n_core_plane, n_plane);
                    emit_wedge_prism(wedge, ew++, inner0, outer0, outer1, n_plane, pts);
                    ++sj;
                } else if (sj == n_arc_side - 1) {
                    const idx_t inner1 = square_side_node(side, si + 1, k, n_core, n_plane);
                    emit_wedge_prism(wedge, ew++, inner0, inner1, outer0, n_plane, pts);
                    ++si;
                } else if (rem_j > rem_i) {
                    const ptrdiff_t p1 = (static_cast<ptrdiff_t>(side) * nq + sj + 1) % ntheta;
                    const idx_t     outer1 =
                            polar_node(0, p1, k, ntheta, n_core_plane, n_plane);
                    emit_wedge_prism(wedge, ew++, inner0, outer0, outer1, n_plane, pts);
                    ++sj;
                } else if (rem_i > rem_j) {
                    const idx_t inner1 = square_side_node(side, si + 1, k, n_core, n_plane);
                    emit_wedge_prism(wedge, ew++, inner0, inner1, outer0, n_plane, pts);
                    ++si;
                } else {
                    const idx_t inner1 = square_side_node(side, si + 1, k, n_core, n_plane);
                    const ptrdiff_t p1 = (static_cast<ptrdiff_t>(side) * nq + sj + 1) % ntheta;
                    const idx_t     outer1 =
                            polar_node(0, p1, k, ntheta, n_core_plane, n_plane);
                    emit_hex_prism(hex, eh++, inner0, inner1, outer1, outer0, n_plane, pts);
                    ++si;
                    ++sj;
                }
            }
        }
    }

    if (eh != n_hex) {
        SMESH_ERROR("create_hex_dominant_cylinder: hex count mismatch %ld vs %ld\n",
                    (long)eh,
                    (long)n_hex);
    }
    if (ew != n_wedge) {
        SMESH_ERROR("create_hex_dominant_cylinder: wedge count mismatch %ld vs %ld\n",
                    (long)ew,
                    (long)n_wedge);
    }

    std::vector<std::shared_ptr<Block>> blocks;
    blocks.push_back(std::make_shared<Block>("hex", HEX8, hex_buf));
    if (n_wedge > 0) {
        blocks.push_back(std::make_shared<Block>("wedge", WEDGE6, wedge_buf));
    }
    return std::make_shared<Mesh>(comm, blocks, points_buf);
}

}  // namespace smesh
