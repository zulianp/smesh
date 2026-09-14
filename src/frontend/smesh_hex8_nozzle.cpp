// Axisymmetric HEX8 nozzle: a bore of piecewise-linear radius that may open, at one axial
// station, into a wider pipe -- the geometry of the FDA benchmark nozzle (Hariharan et al.,
// J. Biomech. Eng. 133(4), 2011), and of any contraction / throat / sudden-expansion device.
//
// The axis is x. Every cross-section is a butterfly O-grid: a square core of n_core x n_core
// cells, n_bore ring layers from the core to the bore wall, and -- downstream of the
// expansion only -- n_outer ring layers from the bore radius out to the expanded radius.
//
// Why a butterfly and not one warped box. A single structured block mapped onto a disc puts
// the four corners of the square on a smooth circle, so the four cells at those corners
// carry an interior angle that tends to 180 degrees. That is topological, not a matter of
// choosing a better map, and it gets worse under refinement: measured, the minimum scaled
// Jacobian of the cross-section is 0.22, 0.12 and 0.06 at 8, 16 and 32 cells across. The
// butterfly's worst angle is the 135 degrees at the core corners, independent of resolution.
//
// Why the sudden expansion is a carved notch. The bore's cells continue straight through the
// expansion plane at the bore radius, and the outer ring exists everywhere in the index grid
// but is kept only downstream. So the expansion face -- the annulus between the bore and the
// expanded radius at the expansion station -- is an exterior face of the kept elements, found
// by any topological skin, and the jet's shear layer lies on a mesh line. Dropping the
// upstream outer-ring elements and compacting the node numbering is exactly what
// create_hex8_lshape does for the backward-facing step, and for the same reasons: the test is
// on integer indices, and a node is kept only if a kept element uses it, so the result is
// orphan-free by construction. See smesh_lshape.cpp for why orphans are not cosmetic.
//
// The elements are not affine. A cross-section of a round pipe cannot be tiled by
// parallelograms, and the cone scales every cell along its length, so a consumer that assumes
// a constant Jacobian per element will be wrong on this mesh.

#include "smesh_mesh.hpp"

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"

#include <cmath>
#include <vector>

namespace smesh {

    std::shared_ptr<Mesh> Mesh::create_hex8_nozzle(const std::shared_ptr<Communicator> &comm,
                                                   const std::vector<geom_t>           &x_breaks,
                                                   const std::vector<geom_t>           &bore_radius,
                                                   const std::vector<ptrdiff_t>        &n_axial,
                                                   const ptrdiff_t                      expansion,
                                                   const geom_t                         expanded_radius,
                                                   const ptrdiff_t                      n_core,
                                                   const ptrdiff_t                      n_bore,
                                                   const ptrdiff_t                      n_outer,
                                                   const geom_t                         core_fraction) {
        const ptrdiff_t n_segments = (ptrdiff_t)n_axial.size();
        if (n_segments < 1 || (ptrdiff_t)x_breaks.size() != n_segments + 1 ||
            (ptrdiff_t)bore_radius.size() != n_segments + 1) {
            SMESH_ERROR("create_hex8_nozzle: need S >= 1 segments, S+1 breaks and S+1 radii");
            return nullptr;
        }
        for (ptrdiff_t s = 0; s < n_segments; ++s) {
            if (n_axial[(size_t)s] < 1 || !(x_breaks[(size_t)s + 1] > x_breaks[(size_t)s])) {
                SMESH_ERROR("create_hex8_nozzle: segment %td needs n_axial >= 1 and increasing x", s);
                return nullptr;
            }
        }
        for (ptrdiff_t k = 0; k <= n_segments; ++k) {
            if (!(bore_radius[(size_t)k] > 0)) {
                SMESH_ERROR("create_hex8_nozzle: bore radius must be positive at break %td", k);
                return nullptr;
            }
        }
        // Even, so that the core has node lines on both symmetry planes and therefore a
        // column of nodes on the axis. Centreline data is what this geometry is compared on,
        // and interpolating it off-axis would be a second approximation on top of the scheme.
        if (n_core < 2 || n_core % 2 != 0 || n_bore < 1 || !(core_fraction > 0) || !(core_fraction < 1)) {
            SMESH_ERROR("create_hex8_nozzle: need an even n_core >= 2, n_bore >= 1 and "
                        "0 < core_fraction < 1");
            return nullptr;
        }
        const bool has_expansion = expansion >= 0 && expansion < n_segments;
        if (has_expansion) {
            if (n_outer < 1) {
                SMESH_ERROR("create_hex8_nozzle: an expansion needs n_outer >= 1");
                return nullptr;
            }
            for (ptrdiff_t k = expansion; k <= n_segments; ++k) {
                if (!(bore_radius[(size_t)k] < expanded_radius)) {
                    SMESH_ERROR("create_hex8_nozzle: downstream of the expansion the bore "
                                "(break %td) must be narrower than the expanded radius", k);
                    return nullptr;
                }
            }
        }

        // ---- the cross-section, in (y, z) ----
        //
        // Core nodes (i, j) first, then one loop of 4 n_core nodes per ring layer l >= 1.
        // Loop 0 is the core boundary and is not stored twice: loop_node(0, m) returns the
        // core node it coincides with.
        const ptrdiff_t n       = n_core;
        const ptrdiff_t n_perim = 4 * n;
        const ptrdiff_t n_loops = n_bore + (has_expansion ? n_outer : 0);  // beyond loop 0
        const ptrdiff_t n_core_nodes = (n + 1) * (n + 1);
        const ptrdiff_t n2d          = n_core_nodes + n_loops * n_perim;

        auto core_node = [n](const ptrdiff_t i, const ptrdiff_t j) { return i + j * (n + 1); };
        // Counter-clockwise from the corner (1, -1).
        auto perim_ij = [n](const ptrdiff_t m, ptrdiff_t &i, ptrdiff_t &j) {
            if (m < n) { i = n; j = m; }
            else if (m < 2 * n) { i = 2 * n - m; j = n; }
            else if (m < 3 * n) { i = 0; j = 3 * n - m; }
            else { i = m - 3 * n; j = 0; }
        };
        auto loop_node = [&](const ptrdiff_t l, const ptrdiff_t m) -> ptrdiff_t {
            const ptrdiff_t mm = ((m % n_perim) + n_perim) % n_perim;
            if (l == 0) {
                ptrdiff_t i, j;
                perim_ij(mm, i, j);
                return core_node(i, j);
            }
            return n_core_nodes + (l - 1) * n_perim + mm;
        };

        // Reference position of every 2D node, split into the part that scales with the bore
        // radius and the part that spans the outer ring, so a plane is placed by
        //     p = rb * bore_part + (R - rb) * outer_part.
        std::vector<double> by((size_t)n2d, 0), bz((size_t)n2d, 0), oy((size_t)n2d, 0), oz((size_t)n2d, 0);
        const double c = (double)core_fraction;
        for (ptrdiff_t j = 0; j <= n; ++j)
            for (ptrdiff_t i = 0; i <= n; ++i) {
                const ptrdiff_t v = core_node(i, j);
                by[(size_t)v]     = c * (-1.0 + 2.0 * (double)i / (double)n);
                bz[(size_t)v]     = c * (-1.0 + 2.0 * (double)j / (double)n);
            }
        for (ptrdiff_t l = 1; l <= n_loops; ++l)
            for (ptrdiff_t m = 0; m < n_perim; ++m) {
                ptrdiff_t i, j;
                perim_ij(m, i, j);
                const double sq_y = -1.0 + 2.0 * (double)i / (double)n;
                const double sq_z = -1.0 + 2.0 * (double)j / (double)n;
                // Equiangular on the circle: the square corner (1,-1) goes to -45 degrees and
                // a side midpoint to 0, so every ring cell subtends the same angle.
                const double th  = -0.25 * M_PI + 0.5 * M_PI * (double)m / (double)n;
                const double cy  = std::cos(th), cz = std::sin(th);
                const ptrdiff_t v = loop_node(l, m);
                if (l <= n_bore) {
                    const double s = (double)l / (double)n_bore;
                    by[(size_t)v]  = (1.0 - s) * c * sq_y + s * cy;
                    bz[(size_t)v]  = (1.0 - s) * c * sq_z + s * cz;
                } else {
                    const double s = (double)(l - n_bore) / (double)n_outer;
                    by[(size_t)v]  = cy;
                    bz[(size_t)v]  = cz;
                    oy[(size_t)v]  = s * cy;
                    oz[(size_t)v]  = s * cz;
                }
            }

        // 2D quads, counter-clockwise in (y, z) so that extrusion along +x gives a positive
        // Jacobian in the corner order below. `outer` marks the quads that exist only
        // downstream of the expansion.
        struct Quad {
            ptrdiff_t v[4];
            bool      outer;
        };
        std::vector<Quad> quads;
        quads.reserve((size_t)(n * n + n_loops * n_perim));
        for (ptrdiff_t j = 0; j < n; ++j)
            for (ptrdiff_t i = 0; i < n; ++i)
                quads.push_back({{core_node(i, j), core_node(i + 1, j), core_node(i + 1, j + 1),
                                  core_node(i, j + 1)},
                                 false});
        for (ptrdiff_t l = 0; l < n_loops; ++l)
            for (ptrdiff_t m = 0; m < n_perim; ++m)
                quads.push_back({{loop_node(l, m), loop_node(l + 1, m), loop_node(l + 1, m + 1),
                                  loop_node(l, m + 1)},
                                 l >= n_bore});

        // ---- the axial planes ----
        ptrdiff_t n_cells_x = 0;
        for (auto na : n_axial) n_cells_x += na;
        const ptrdiff_t      n_planes = n_cells_x + 1;
        std::vector<double>  plane_x((size_t)n_planes), plane_rb((size_t)n_planes);
        std::vector<ptrdiff_t> cell_segment((size_t)n_cells_x);
        {
            ptrdiff_t k = 0;
            for (ptrdiff_t s = 0; s < n_segments; ++s) {
                const ptrdiff_t na = n_axial[(size_t)s];
                for (ptrdiff_t a = 0; a < na; ++a, ++k) {
                    const double t     = (double)a / (double)na;
                    plane_x[(size_t)k]  = (1 - t) * x_breaks[(size_t)s] + t * x_breaks[(size_t)s + 1];
                    plane_rb[(size_t)k] = (1 - t) * bore_radius[(size_t)s] + t * bore_radius[(size_t)s + 1];
                    cell_segment[(size_t)k] = s;
                }
            }
            // Written from the break itself, not accumulated, so the end planes are exact.
            plane_x[(size_t)k]  = x_breaks[(size_t)n_segments];
            plane_rb[(size_t)k] = bore_radius[(size_t)n_segments];
        }

        auto kept = [&](const Quad &q, const ptrdiff_t k) {
            return !q.outer || cell_segment[(size_t)k] >= expansion;
        };

        // ---- pass 1: which grid nodes survive ----
        const ptrdiff_t   n_grid_nodes = n_planes * n2d;
        std::vector<char> used((size_t)n_grid_nodes, 0);
        ptrdiff_t         n_kept_elements = 0;
        for (ptrdiff_t k = 0; k < n_cells_x; ++k)
            for (const auto &q : quads) {
                if (!kept(q, k)) continue;
                ++n_kept_elements;
                for (int a = 0; a < 4; ++a) {
                    used[(size_t)(k * n2d + q.v[a])]     = 1;
                    used[(size_t)((k + 1) * n2d + q.v[a])] = 1;
                }
            }

        std::vector<idx_t> old2new((size_t)n_grid_nodes, (idx_t)-1);
        ptrdiff_t          n_kept_nodes = 0;
        for (ptrdiff_t g = 0; g < n_grid_nodes; ++g)
            if (used[(size_t)g]) old2new[(size_t)g] = (idx_t)n_kept_nodes++;

        auto points_buffer   = create_host_buffer<geom_t>(3, n_kept_nodes);
        auto elements_buffer = create_host_buffer<idx_t>(8, n_kept_elements);
        auto points          = points_buffer->data();
        auto elements        = elements_buffer->data();

        const double R = has_expansion ? (double)expanded_radius : 0.0;
        for (ptrdiff_t k = 0; k < n_planes; ++k) {
            const double rb = plane_rb[(size_t)k];
            for (ptrdiff_t v = 0; v < n2d; ++v) {
                const idx_t nn = old2new[(size_t)(k * n2d + v)];
                if (nn < 0) continue;
                points[0][nn] = (geom_t)plane_x[(size_t)k];
                points[1][nn] = (geom_t)(rb * by[(size_t)v] + (R - rb) * oy[(size_t)v]);
                points[2][nn] = (geom_t)(rb * bz[(size_t)v] + (R - rb) * oz[(size_t)v]);
            }
        }

        // Same corner order as mesh_fill_hex8_cube, with the quad in place of the (x, y) face
        // and the axial direction in place of z.
        ptrdiff_t e = 0;
        for (ptrdiff_t k = 0; k < n_cells_x; ++k)
            for (const auto &q : quads) {
                if (!kept(q, k)) continue;
                for (int a = 0; a < 4; ++a) {
                    elements[a][e]     = old2new[(size_t)(k * n2d + q.v[a])];
                    elements[a + 4][e] = old2new[(size_t)((k + 1) * n2d + q.v[a])];
                }
                ++e;
            }
        SMESH_ASSERT(e == n_kept_elements);

        auto block = std::make_shared<Block>();
        block->set_name("fluid");
        block->set_element_type(HEX8);
        block->set_elements(elements_buffer);

        return std::make_shared<Mesh>(comm, std::vector<std::shared_ptr<Block>>{block},
                                      points_buffer);
    }

}  // namespace smesh
