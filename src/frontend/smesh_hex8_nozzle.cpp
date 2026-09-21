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
//
// The map is written once, below, as a function of CONTINUOUS grid parameters -- (i, j) across
// the core, (l, m) around a ring, a plane coordinate along x -- and the integer grid is just
// that map sampled at integers. That is what lets warp_semistructured_hex8_nozzle place the
// micro nodes of a refined mesh on the nozzle itself instead of on the chords between macro
// corners: they are the same map sampled at the fractions the lattice sits at.

#include "smesh_mesh.hpp"

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_hex8_nozzle_grid.hpp"
#include "smesh_sshex8.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_create.hpp"
#endif

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
                                                   const geom_t                         core_fraction,
                                                   const geom_t                         radial_grading,
                                                   const geom_t                         axial_grading) {
        if (!nozzle_arguments_valid(x_breaks, bore_radius, n_axial, expansion, expanded_radius, n_core, n_bore,
                                    n_outer, core_fraction, radial_grading, axial_grading, "create_hex8_nozzle"))
            return nullptr;
#ifdef SMESH_ENABLE_MPI
        if (comm && comm->size() > 1) {
            int       nxe = 0, sdim = 0;
            ptrdiff_t n_local_e = 0, n_global_e = 0, n_local_n = 0, n_global_n = 0;
            idx_t  **elems  = nullptr;
            geom_t **points = nullptr;
            if (hex8_nozzle_create_distributed<idx_t, geom_t>(comm->get(),
                                                              x_breaks.data(),
                                                              bore_radius.data(),
                                                              n_axial.data(),
                                                              (ptrdiff_t)n_axial.size(),
                                                              expansion,
                                                              expanded_radius,
                                                              n_core,
                                                              n_bore,
                                                              n_outer,
                                                              core_fraction,
                                                              radial_grading,
                                                              axial_grading,
                                                              &nxe,
                                                              &n_local_e,
                                                              &n_global_e,
                                                              &elems,
                                                              &sdim,
                                                              &n_local_n,
                                                              &n_global_n,
                                                              &points) != SMESH_SUCCESS) {
                return nullptr;
            }
            auto mesh = Mesh::wrap_create_parallel(comm, HEX8, nxe, n_local_e, n_global_e, elems, sdim, n_local_n,
                                                   n_global_n, points, ISOPARAMETRIC);
            if (mesh && mesh->n_blocks() > 0) {
                mesh->block(0)->set_name("fluid");
            }
            return mesh;
        }
#endif
        const ptrdiff_t n_segments    = (ptrdiff_t)n_axial.size();
        const bool      has_expansion = expansion >= 0 && expansion < n_segments;

        const NozzleSection sec(n_core, n_bore, n_outer, has_expansion, (double)core_fraction,
                                (double)radial_grading);
        const ptrdiff_t     n2d = sec.n2d;
        std::vector<double> by, bz, oy, oz;
        sec.reference(by, bz, oy, oz);

        // 2D quads, counter-clockwise in (y, z) so that extrusion along +x gives a positive
        // Jacobian in the corner order below. `outer` marks the quads that exist only
        // downstream of the expansion.
        std::vector<NozzleQuad> quads;
        nozzle_build_quads(sec, quads);

        const NozzlePlanes pl(x_breaks, bore_radius, n_axial, (double)axial_grading);
        const ptrdiff_t    n_cells_x = pl.n_cells_x;
        const ptrdiff_t    n_planes  = pl.n_planes;

        auto kept = [&](const NozzleQuad &q, const ptrdiff_t k) {
            return !q.outer || pl.cell_segment[(size_t)k] >= expansion;
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
            const double rb = pl.plane_rb[(size_t)k];
            for (ptrdiff_t v = 0; v < n2d; ++v) {
                const idx_t nn = old2new[(size_t)(k * n2d + v)];
                if (nn < 0) continue;
                points[0][nn] = (geom_t)pl.plane_x[(size_t)k];
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
        block->set_geom_map(ISOPARAMETRIC);

        return std::make_shared<Mesh>(comm, std::vector<std::shared_ptr<Block>>{block},
                                      points_buffer);
    }

    // Micro nodes on the nozzle, not on the chords between macro corners.
    //
    // to_semistructured places each micro node by trilinear interpolation of its macro element's
    // corners. The corners are on the nozzle; the interpolated nodes lie on the polygon they
    // span, which cuts the throat's area by what the macro mesh's ring cells cut, however fine
    // the lattice. Measured on the FDA nozzle at throat Re 500, macro core 2 at level 4 (116,212
    // dof) carried a 3.7% inflow-flux deficit and a throat centreline 4-8% above the PIV data,
    // and only a finer MACRO mesh helped -- which costs the multigrid hierarchy its depth, and a
    // coarsest level of 15,740 dof makes the dense coarse LU the whole run.
    //
    // Here each macro element's corners are matched back to their grid parameters, and every
    // micro node is placed by the generator's own map at the trilinear interpolation of those
    // PARAMETERS. A lattice of level L is then exactly create_hex8_nozzle at L times the
    // resolution, node for node, and every coarser level of the hierarchy, whose nodes are a
    // subset, is on the nozzle too.
    int Mesh::warp_semistructured_hex8_nozzle(const std::shared_ptr<Mesh>  &sshex,
                                              const std::vector<geom_t>    &x_breaks,
                                              const std::vector<geom_t>    &bore_radius,
                                              const std::vector<ptrdiff_t> &n_axial,
                                              const ptrdiff_t               expansion,
                                              const geom_t                  expanded_radius,
                                              const ptrdiff_t               n_core,
                                              const ptrdiff_t               n_bore,
                                              const ptrdiff_t               n_outer,
                                              const geom_t                  core_fraction,
                                              const geom_t                  radial_grading,
                                              const geom_t                  axial_grading) {
        if (!sshex || !nozzle_arguments_valid(x_breaks, bore_radius, n_axial, expansion, expanded_radius, n_core,
                                              n_bore, n_outer, core_fraction, radial_grading, axial_grading,
                                              "warp_semistructured_hex8_nozzle"))
            return SMESH_FAILURE;
        const int nxe = sshex->n_nodes_per_element(0);
        int       L   = (int)std::lround(std::cbrt((double)nxe)) - 1;
        if (L < 1 || (L + 1) * (L + 1) * (L + 1) != nxe) {
            SMESH_ERROR("warp_semistructured_hex8_nozzle: %d nodes per element is not a hex lattice", nxe);
            return SMESH_FAILURE;
        }
        const ptrdiff_t n_segments    = (ptrdiff_t)n_axial.size();
        const bool      has_expansion = expansion >= 0 && expansion < n_segments;
        const NozzleSection sec(n_core, n_bore, n_outer, has_expansion, (double)core_fraction,
                                (double)radial_grading);
        const NozzlePlanes  pl(x_breaks, bore_radius, n_axial, (double)axial_grading);
        std::vector<double> by, bz, oy, oz;
        sec.reference(by, bz, oy, oz);
        const double R = has_expansion ? (double)expanded_radius : 0.0;

        const ptrdiff_t          nelements = sshex->n_elements(0);
        const auto *const       *el        = sshex->elements(0)->data();
        const auto *const        pts       = sshex->points()->data();
        const int corner[8] = {sshex8_lidx(L, 0, 0, 0), sshex8_lidx(L, L, 0, 0), sshex8_lidx(L, L, L, 0),
                               sshex8_lidx(L, 0, L, 0), sshex8_lidx(L, 0, 0, L), sshex8_lidx(L, L, 0, L),
                               sshex8_lidx(L, L, L, L), sshex8_lidx(L, 0, L, L)};

        auto nearest_plane = [&](const double x) {
            ptrdiff_t best = 0;
            for (ptrdiff_t k = 1; k < pl.n_planes; ++k)
                if (std::fabs(pl.plane_x[(size_t)k] - x) < std::fabs(pl.plane_x[(size_t)best] - x)) best = k;
            return best;
        };
        // Nearest 2D grid node of a macro corner in plane k, with its distance.
        auto nearest_node = [&](const ptrdiff_t k, const double y, const double z, double &dist) {
            const double rb   = pl.plane_rb[(size_t)k];
            ptrdiff_t    best = -1;
            dist              = 1e300;
            for (ptrdiff_t v = 0; v < sec.n2d; ++v) {
                const double dy = rb * by[(size_t)v] + (R - rb) * oy[(size_t)v] - y;
                const double dz = rb * bz[(size_t)v] + (R - rb) * oz[(size_t)v] - z;
                const double d2 = dy * dy + dz * dz;
                if (d2 < dist) {
                    dist = d2;
                    best = v;
                }
            }
            dist = std::sqrt(dist);
            return best;
        };

        for (ptrdiff_t e = 0; e < nelements; ++e) {
            const idx_t  g0 = el[corner[0]][e], g1 = el[corner[1]][e], g3 = el[corner[3]][e], g4 = el[corner[4]][e];
            const ptrdiff_t k0 = nearest_plane((double)pts[0][g0]);
            const ptrdiff_t k4 = nearest_plane((double)pts[0][g4]);
            if (k4 != k0 + 1 || nearest_plane((double)pts[0][g1]) != k0 || nearest_plane((double)pts[0][g3]) != k0) {
                SMESH_ERROR("warp_semistructured_hex8_nozzle: macro element %td does not span one axial cell "
                            "of this nozzle; were the arguments the ones the mesh was created with?", e);
                return SMESH_FAILURE;
            }
            const double tol = 1e-4 * pl.plane_rb[(size_t)k0];
            double       d0, d1, d3;
            const ptrdiff_t v0 = nearest_node(k0, (double)pts[1][g0], (double)pts[2][g0], d0);
            const ptrdiff_t v1 = nearest_node(k0, (double)pts[1][g1], (double)pts[2][g1], d1);
            const ptrdiff_t v3 = nearest_node(k0, (double)pts[1][g3], (double)pts[2][g3], d3);
            if (d0 > tol || d1 > tol || d3 > tol) {
                SMESH_ERROR("warp_semistructured_hex8_nozzle: macro element %td has a corner off the nozzle's "
                            "grid (%g, %g, %g from it)", e, d0, d1, d3);
                return SMESH_FAILURE;
            }

            // Core element: its lower quad is (i, j), (i+1, j), (i+1, j+1), (i, j+1). Ring
            // element: (l, m), (l+1, m), (l+1, m+1), (l, m+1), whose v1 is always a loop node.
            bool      is_core = false;
            double    a0 = 0, b0 = 0;
            ptrdiff_t l_cell = 0;
            if (v1 < sec.n_core_nodes) {
                const ptrdiff_t i0 = v0 % (sec.n + 1), j0 = v0 / (sec.n + 1);
                if (v0 >= sec.n_core_nodes || v1 != sec.core_node(i0 + 1, j0) || v3 != sec.core_node(i0, j0 + 1)) {
                    SMESH_ERROR("warp_semistructured_hex8_nozzle: core macro element %td is not in grid order", e);
                    return SMESH_FAILURE;
                }
                is_core = true;
                a0      = (double)i0;
                b0      = (double)j0;
            } else {
                const ptrdiff_t l1 = (v1 - sec.n_core_nodes) / sec.n_perim + 1;
                const ptrdiff_t m1 = (v1 - sec.n_core_nodes) % sec.n_perim;
                if (v0 != sec.loop_node(l1 - 1, m1) || v3 != sec.loop_node(l1 - 1, m1 + 1)) {
                    SMESH_ERROR("warp_semistructured_hex8_nozzle: ring macro element %td is not in grid order", e);
                    return SMESH_FAILURE;
                }
                l_cell = l1 - 1;
                a0     = (double)l_cell;
                b0     = (double)m1;
            }
            const bool outer = !is_core && l_cell >= sec.n_bore;

            for (int zi = 0; zi <= L; ++zi) {
                const double t = (double)zi / (double)L;
                // Through the map at the fractional plane coordinate rather than by
                // interpolating the two bounding planes: with grading on, the chord between
                // them is not the nozzle.
                double px, rb;
                pl.at((double)k0 + t, px, rb);
                for (int yi = 0; yi <= L; ++yi)
                    for (int xi = 0; xi <= L; ++xi) {
                        const double a = a0 + (double)xi / (double)L;
                        const double b = b0 + (double)yi / (double)L;
                        double       qy, qz, ry, rz;
                        if (is_core)
                            sec.core(a, b, qy, qz, ry, rz);
                        else
                            sec.ring(a, b, outer, qy, qz, ry, rz);
                        const idx_t node = el[sshex8_lidx(L, xi, yi, zi)][e];
                        pts[0][node]     = (geom_t)px;
                        pts[1][node]     = (geom_t)(rb * qy + (R - rb) * ry);
                        pts[2][node]     = (geom_t)(rb * qz + (R - rb) * rz);
                    }
            }
        }
        for (size_t b = 0; b < sshex->n_blocks(); ++b) {
            sshex->block(b)->set_geom_map(ISOPARAMETRIC);
        }
        return SMESH_SUCCESS;
    }

}  // namespace smesh
