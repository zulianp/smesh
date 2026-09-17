// L-shaped (backward-facing step) HEX8 mesh.
//
// The domain is the box [0,xmax] x [0,ymax] x [0,zmax] with the notch
// [0,step_x] x [0,step_y] x [0,zmax] removed, i.e.
//
//     Omega = ( ([0,xmax] x [step_y,ymax]) union ([step_x,xmax] x [0,step_y]) ) x [0,zmax]
//
// which is the geometry of the 3D backward-facing step in Farrell, Mitchell & Wechsung,
// SIAM J. Sci. Comput. 41(5), A3073-A3096, with the defaults below reproducing their
// Omega = ( ([0,10] x [1,2]) union ([1,10] x [0,1]) ) x [0,1] exactly.
//
// Built by generating the full structured grid and dropping the notch elements, then
// compacting the node numbering. Three details in that sentence are load-bearing.
//
// The notch test is on integer indices, not element centroids. A float test on a grid
// where the step does not fall on a cell boundary produces a ragged step that still looks
// plausible in a viewer, and the resolution is then silently not what was asked for. The
// integer test cannot do that, and the constructor rejects a resolution that would not put
// a grid line on the step.
//
// The compaction is done here rather than by Mesh::renumber_nodes(), which is a permutation
// and not a compaction: it allocates the same node count and asserts that every node is
// referenced. That assert is compiled out under NDEBUG, so an orphan node writes to
// new_points[d][-1] -- heap corruption in a release build.
//
// Orphans are not cosmetic further downstream either. sshex8_generate_elements starts its
// index_base at the macro mesh's node count, so an orphan consumes a semi-structured node
// id whose coordinates are never written by sshex8_fill_points. It stays at (0,0,0) from
// the calloc, which then satisfies the CVFEM driver's on_plane(x, 0, Lx) test and becomes a
// spurious interior Dirichlet node -- a defect that produces a converged, plausible, wrong
// answer rather than a failure.

#include "smesh_mesh.hpp"

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_create.hpp"
#endif

#include <cmath>
#include <vector>

namespace smesh {

    std::shared_ptr<Mesh> Mesh::create_hex8_lshape(const std::shared_ptr<Communicator> &comm,
                                                   const ptrdiff_t                      nx,
                                                   const ptrdiff_t                      ny,
                                                   const ptrdiff_t                      nz,
                                                   const geom_t                         xmax,
                                                   const geom_t                         ymax,
                                                   const geom_t                         zmax,
                                                   const geom_t                         step_x,
                                                   const geom_t                         step_y) {
        if (nx < 1 || ny < 1 || nz < 1) {
            SMESH_ERROR("create_hex8_lshape: resolution must be positive (%td, %td, %td)", nx, ny, nz);
            return nullptr;
        }
        if (!(xmax > 0) || !(ymax > 0) || !(zmax > 0)) {
            SMESH_ERROR("create_hex8_lshape: extents must be positive");
            return nullptr;
        }
        if (!(step_x > 0) || !(step_y > 0) || step_x >= xmax || step_y >= ymax) {
            SMESH_ERROR("create_hex8_lshape: the step must lie strictly inside the box");
            return nullptr;
        }

        // The step has to land on a grid line, or the notch is not the notch that was asked
        // for. Reject rather than round: a mesh whose step is half a cell out of place is
        // exactly the kind of error that shows up later as an unexplained discrepancy.
        const double  fx  = (double)step_x / (double)xmax * (double)nx;
        const double  fy  = (double)step_y / (double)ymax * (double)ny;
        const ptrdiff_t nxs = (ptrdiff_t)(fx + 0.5);
        const ptrdiff_t nys = (ptrdiff_t)(fy + 0.5);
        if (std::fabs(fx - (double)nxs) > 1e-9 || std::fabs(fy - (double)nys) > 1e-9) {
            SMESH_ERROR("create_hex8_lshape: the step must fall on a grid line; "
                        "step_x/xmax*nx = %g and step_y/ymax*ny = %g must both be integers",
                        fx, fy);
            return nullptr;
        }

#ifdef SMESH_ENABLE_MPI
        if (comm && comm->size() > 1) {
            int       nxe = 0, sdim = 0;
            ptrdiff_t n_local_e = 0, n_global_e = 0, n_local_n = 0, n_global_n = 0;
            idx_t  **elems  = nullptr;
            geom_t **points = nullptr;
            if (hex8_lshape_create_distributed<idx_t, geom_t>(comm->get(),
                                                              nx,
                                                              ny,
                                                              nz,
                                                              xmax,
                                                              ymax,
                                                              zmax,
                                                              nxs,
                                                              nys,
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
                                                   n_global_n, points, AXIS_ALIGNED);
            if (mesh && mesh->n_blocks() > 0) {
                mesh->block(0)->set_name("fluid");
            }
            return mesh;
        }
#endif

        const ptrdiff_t ldz = (ny + 1) * (nx + 1);
        const ptrdiff_t ldy = (nx + 1);
        const ptrdiff_t ldx = 1;
        const ptrdiff_t n_grid_nodes = (nx + 1) * (ny + 1) * (nz + 1);

        auto in_notch = [nxs, nys](const ptrdiff_t xi, const ptrdiff_t yi) {
            return xi < nxs && yi < nys;
        };

        // Pass 1: which grid nodes survive. A node is kept iff some kept element uses it,
        // which is what makes the result orphan-free by construction rather than by check.
        std::vector<char> used((size_t)n_grid_nodes, 0);
        ptrdiff_t         n_kept_elements = 0;
        for (ptrdiff_t zi = 0; zi < nz; zi++)
            for (ptrdiff_t yi = 0; yi < ny; yi++)
                for (ptrdiff_t xi = 0; xi < nx; xi++) {
                    if (in_notch(xi, yi)) continue;
                    ++n_kept_elements;
                    for (int dz = 0; dz < 2; ++dz)
                        for (int dy = 0; dy < 2; ++dy)
                            for (int dx = 0; dx < 2; ++dx)
                                used[(size_t)((xi + dx) * ldx + (yi + dy) * ldy + (zi + dz) * ldz)] = 1;
                }

        std::vector<idx_t> old2new((size_t)n_grid_nodes, (idx_t)-1);
        ptrdiff_t          n_kept_nodes = 0;
        for (ptrdiff_t i = 0; i < n_grid_nodes; ++i)
            if (used[(size_t)i]) old2new[(size_t)i] = (idx_t)n_kept_nodes++;

        auto points_buffer   = create_host_buffer<geom_t>(3, n_kept_nodes);
        auto elements_buffer = create_host_buffer<idx_t>(8, n_kept_elements);
        auto points          = points_buffer->data();
        auto elements        = elements_buffer->data();

        const double hx = (double)xmax / (double)nx;
        const double hy = (double)ymax / (double)ny;
        const double hz = (double)zmax / (double)nz;

        // Kept nodes are visited in the same lexicographic order as the full grid, so the
        // compacted numbering stays monotone in (zi, yi, xi) and locality is preserved.
        for (ptrdiff_t zi = 0; zi <= nz; zi++)
            for (ptrdiff_t yi = 0; yi <= ny; yi++)
                for (ptrdiff_t xi = 0; xi <= nx; xi++) {
                    const ptrdiff_t g = xi * ldx + yi * ldy + zi * ldz;
                    const idx_t     n = old2new[(size_t)g];
                    if (n < 0) continue;
                    points[0][n] = (geom_t)(xi * hx);
                    points[1][n] = (geom_t)(yi * hy);
                    points[2][n] = (geom_t)(zi * hz);
                }

        // Same corner ordering as mesh_fill_hex8_cube: 0:(0,0,0) 1:(1,0,0) 2:(1,1,0)
        // 3:(0,1,0), then the same four shifted by one in z.
        ptrdiff_t e = 0;
        for (ptrdiff_t zi = 0; zi < nz; zi++)
            for (ptrdiff_t yi = 0; yi < ny; yi++)
                for (ptrdiff_t xi = 0; xi < nx; xi++) {
                    if (in_notch(xi, yi)) continue;
                    const ptrdiff_t c[8] = {(xi + 0) * ldx + (yi + 0) * ldy + (zi + 0) * ldz,
                                            (xi + 1) * ldx + (yi + 0) * ldy + (zi + 0) * ldz,
                                            (xi + 1) * ldx + (yi + 1) * ldy + (zi + 0) * ldz,
                                            (xi + 0) * ldx + (yi + 1) * ldy + (zi + 0) * ldz,
                                            (xi + 0) * ldx + (yi + 0) * ldy + (zi + 1) * ldz,
                                            (xi + 1) * ldx + (yi + 0) * ldy + (zi + 1) * ldz,
                                            (xi + 1) * ldx + (yi + 1) * ldy + (zi + 1) * ldz,
                                            (xi + 0) * ldx + (yi + 1) * ldy + (zi + 1) * ldz};
                    for (int a = 0; a < 8; ++a) elements[a][e] = old2new[(size_t)c[a]];
                    ++e;
                }
        SMESH_ASSERT(e == n_kept_elements);

        // Built through the public block constructor rather than by writing Mesh::Impl
        // directly: Impl is opaque outside smesh_mesh.cpp, and there is no reason for a
        // generator to need private access.
        auto block = std::make_shared<Block>();
        block->set_name("fluid");
        block->set_element_type(HEX8);
        block->set_elements(elements_buffer);
        block->set_geom_map(AXIS_ALIGNED);

        return std::make_shared<Mesh>(comm, std::vector<std::shared_ptr<Block>>{block},
                                      points_buffer);
    }

}  // namespace smesh
