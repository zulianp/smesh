#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "smesh_fff.hpp"
#include "smesh_jacobians.hpp"
#include "smesh_kernel_data.hpp"
#include "smesh_mesh.hpp"
#include "smesh_mesh_reorder.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_ops.hpp"
#include "smesh_parametrization.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sfc.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#include <limits>
#endif

using namespace smesh;

namespace {

constexpr geom_t kTol = sizeof(geom_t) == 8 ? static_cast<geom_t>(1e-10)
                                            : static_cast<geom_t>(1e-5);

std::shared_ptr<Mesh> split_first_half(const std::shared_ptr<Mesh> &mesh) {
    auto out = mesh->clone();
    const ptrdiff_t n = out->n_elements(0);
    const ptrdiff_t n_split = n / 2;
    auto parents = create_host_buffer<element_idx_t>(static_cast<size_t>(n_split));
    for (ptrdiff_t i = 0; i < n_split; ++i) {
        parents->data()[i] = static_cast<element_idx_t>(i);
    }
    if (out->split_block(parents, "part0", 0) != SMESH_SUCCESS) {
        return nullptr;
    }
    return out;
}

int check_multiblock_sfc(Mesh &mesh, const char *ordering) {
    const size_t n_blocks = mesh.n_blocks();
    SMESH_TEST_ASSERT(n_blocks > 1);
    const ptrdiff_t n_nodes = mesh.n_nodes();
    const ptrdiff_t n_total = mesh.n_elements();
    std::vector<ptrdiff_t> n_e(n_blocks);
    for (size_t b = 0; b < n_blocks; ++b) {
        n_e[b] = mesh.n_elements(static_cast<block_idx_t>(b));
    }

    SMESH_TEST_EQ(SFC(ordering).reorder(mesh), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh.n_blocks(), n_blocks);
    SMESH_TEST_EQ(mesh.n_nodes(), n_nodes);
    SMESH_TEST_EQ(mesh.n_elements(), n_total);
    for (size_t b = 0; b < n_blocks; ++b) {
        SMESH_TEST_EQ(mesh.n_elements(static_cast<block_idx_t>(b)), n_e[b]);
    }

    std::vector<char> seen(static_cast<size_t>(n_nodes), 0);
    for (size_t b = 0; b < n_blocks; ++b) {
        const block_idx_t bid = static_cast<block_idx_t>(b);
        const int nxe = mesh.n_nodes_per_element(bid);
        const ptrdiff_t ne = mesh.n_elements(bid);
        idx_t *const *const elems = mesh.elements(bid)->data();
        for (ptrdiff_t e = 0; e < ne; ++e) {
            for (int d = 0; d < nxe; ++d) {
                const idx_t n = elems[d][e];
                SMESH_TEST_ASSERT(n >= 0 && n < n_nodes);
                seen[static_cast<size_t>(n)] = 1;
            }
        }
    }
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        SMESH_TEST_EQ(static_cast<int>(seen[static_cast<size_t>(i)]), 1);
    }

    auto bary = create_host_buffer<geom_t>(3, static_cast<size_t>(n_total));
    geom_t **d_b = bary->data();
    const int sdim = mesh.spatial_dimension();
    geom_t *const *const pts = mesh.points()->data();
    ptrdiff_t off = 0;
    std::vector<ptrdiff_t> offset(n_blocks + 1, 0);
    for (size_t b = 0; b < n_blocks; ++b) {
        offset[b] = off;
        const block_idx_t bid = static_cast<block_idx_t>(b);
        const ptrdiff_t ne = mesh.n_elements(bid);
        if (ne > 0) {
            geom_t *slice[3] = {d_b[0] + off, d_b[1] + off, d_b[2] + off};
            barycenters(mesh.n_nodes_per_element(bid), ne, mesh.elements(bid)->data(),
                        sdim, pts, slice);
        }
        off += ne;
    }
    offset[n_blocks] = off;
    if (sdim < 3) {
        std::memset(d_b[2], 0, static_cast<size_t>(n_total) * sizeof(geom_t));
    }
    auto enc = create_host_buffer<u32>(static_cast<size_t>(n_total));
    if (std::strcmp(ordering, "morton3") == 0) {
        SMESH_TEST_EQ(encode_morton3<geom_t>(n_total, d_b[0], d_b[1], d_b[2], enc->data()),
                      SMESH_SUCCESS);
    } else if (std::strcmp(ordering, "hilbert3") == 0) {
        SMESH_TEST_EQ(encode_hilbert3<geom_t>(n_total, d_b[0], d_b[1], d_b[2], enc->data()),
                      SMESH_SUCCESS);
    } else {
        return SMESH_TEST_FAILURE;
    }
    const u32 *keys = enc->data();
    for (size_t b = 0; b < n_blocks; ++b) {
        for (ptrdiff_t e = offset[b] + 1; e < offset[b + 1]; ++e) {
            SMESH_TEST_ASSERT(keys[e - 1] <= keys[e]);
        }
    }
    return SMESH_TEST_SUCCESS;
}

#ifdef SMESH_ENABLE_MPI
int encode_keys(const char *ordering, const ptrdiff_t n, geom_t **d_b,
                const geom_t *bmin, const geom_t *bmax, u32 *keys) {
    if (std::strcmp(ordering, "morton3") == 0) {
        return encode_morton3<geom_t>(n, d_b[0], d_b[1], d_b[2], bmin[0], bmax[0],
                                      bmin[1], bmax[1], bmin[2], bmax[2], keys);
    }
    if (std::strcmp(ordering, "hilbert3") == 0) {
        return encode_hilbert3<geom_t>(n, d_b[0], d_b[1], d_b[2], bmin[0], bmax[0],
                                       bmin[1], bmax[1], bmin[2], bmax[2], keys);
    }
    return SMESH_FAILURE;
}

int check_distributed_sfc(Mesh &mesh, const char *ordering) {
    SMESH_TEST_ASSERT(mesh.is_distributed());
    const size_t n_blocks = mesh.n_blocks();
    const ptrdiff_t n_nodes = mesh.n_nodes();
    const ptrdiff_t n_total = mesh.n_elements();
    auto dist = mesh.distributed();
    const ptrdiff_t n_owned_all = dist->n_elements_owned();
    const ptrdiff_t n_ghosts_all = dist->n_elements_ghosts();
    const ptrdiff_t n_nodes_owned = dist->n_nodes_owned();
    std::vector<ptrdiff_t> n_e(n_blocks), n_owned(n_blocks), n_shared(n_blocks),
        n_ghosts(n_blocks);
    for (size_t b = 0; b < n_blocks; ++b) {
        auto block = mesh.block(b);
        n_e[b] = block->n_elements();
        n_owned[b] = block->n_elements_owned();
        n_shared[b] = block->n_elements_shared();
        n_ghosts[b] = block->n_elements_ghosts();
    }

    SMESH_TEST_EQ(SFC(ordering).reorder(mesh), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh.n_blocks(), n_blocks);
    SMESH_TEST_EQ(mesh.n_nodes(), n_nodes);
    SMESH_TEST_EQ(mesh.n_elements(), n_total);
    SMESH_TEST_EQ(mesh.distributed()->n_elements_owned(), n_owned_all);
    SMESH_TEST_EQ(mesh.distributed()->n_elements_ghosts(), n_ghosts_all);
    SMESH_TEST_EQ(mesh.distributed()->n_nodes_owned(), n_nodes_owned);
    for (size_t b = 0; b < n_blocks; ++b) {
        auto block = mesh.block(b);
        SMESH_TEST_EQ(block->n_elements(), n_e[b]);
        SMESH_TEST_EQ(block->n_elements_owned(), n_owned[b]);
        SMESH_TEST_EQ(block->n_elements_shared(), n_shared[b]);
        SMESH_TEST_EQ(block->n_elements_ghosts(), n_ghosts[b]);
        SMESH_TEST_EQ(block->n_elements_owned() + block->n_elements_ghosts(),
                      block->n_elements());
    }

    for (size_t b = 0; b < n_blocks; ++b) {
        auto block = mesh.block(b);
        const int nxe = block->n_nodes_per_element();
        idx_t *const *const elems = block->elements()->data();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            for (int d = 0; d < nxe; ++d) {
                SMESH_TEST_ASSERT(elems[d][e] >= 0);
                SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(elems[d][e]) < n_nodes);
            }
        }
        if (n_owned[b] > 0) {
            SMESH_TEST_ASSERT(block->element_mapping() != nullptr);
            SMESH_TEST_EQ(static_cast<ptrdiff_t>(block->element_mapping()->size()),
                          n_owned[b]);
        }
    }

    if (n_total > 0) {
        auto bary = create_host_buffer<geom_t>(3, static_cast<size_t>(n_total));
        geom_t **d_b = bary->data();
        const int sdim = mesh.spatial_dimension();
        geom_t *const *const pts = mesh.points()->data();
        ptrdiff_t off = 0;
        std::vector<ptrdiff_t> offset(n_blocks + 1, 0);
        for (size_t b = 0; b < n_blocks; ++b) {
            offset[b] = off;
            const block_idx_t bid = static_cast<block_idx_t>(b);
            const ptrdiff_t ne = mesh.n_elements(bid);
            if (ne > 0) {
                geom_t *slice[3] = {d_b[0] + off, d_b[1] + off, d_b[2] + off};
                barycenters(mesh.n_nodes_per_element(bid), ne, mesh.elements(bid)->data(),
                            sdim, pts, slice);
            }
            off += ne;
        }
        offset[n_blocks] = off;
        if (sdim < 3) {
            std::memset(d_b[2], 0, static_cast<size_t>(n_total) * sizeof(geom_t));
        }
        geom_t lmin[3] = {std::numeric_limits<geom_t>::max(),
                          std::numeric_limits<geom_t>::max(),
                          std::numeric_limits<geom_t>::max()};
        geom_t lmax[3] = {std::numeric_limits<geom_t>::lowest(),
                          std::numeric_limits<geom_t>::lowest(),
                          std::numeric_limits<geom_t>::lowest()};
        for (int d = 0; d < 3; ++d) {
            for (ptrdiff_t i = 0; i < n_total; ++i) {
                lmin[d] = std::min(lmin[d], d_b[d][i]);
                lmax[d] = std::max(lmax[d], d_b[d][i]);
            }
        }
        geom_t gmin[3], gmax[3];
        SMESH_MPI_CATCH(MPI_Allreduce(lmin, gmin, 3, mpi_type<geom_t>(), MPI_MIN,
                                      mesh.comm()->get()));
        SMESH_MPI_CATCH(MPI_Allreduce(lmax, gmax, 3, mpi_type<geom_t>(), MPI_MAX,
                                      mesh.comm()->get()));
        auto enc = create_host_buffer<u32>(static_cast<size_t>(n_total));
        SMESH_TEST_EQ(encode_keys(ordering, n_total, d_b, gmin, gmax, enc->data()),
                      SMESH_SUCCESS);
        const u32 *keys = enc->data();
        for (size_t b = 0; b < n_blocks; ++b) {
            const ptrdiff_t lo = offset[b];
            const ptrdiff_t n_ons = n_owned[b] - n_shared[b];
            const ptrdiff_t segs[4] = {lo, lo + n_ons, lo + n_owned[b], offset[b + 1]};
            for (int s = 0; s < 3; ++s) {
                for (ptrdiff_t e = segs[s] + 1; e < segs[s + 1]; ++e) {
                    SMESH_TEST_ASSERT(keys[e - 1] <= keys[e]);
                }
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}
#endif

}  // namespace

static int test_quad4_jacobian_fff_analytic() {
    auto mesh = Mesh::create_quad4_square(Communicator::self(), 1, 1);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(mesh->spatial_dimension(), 2);
    SMESH_TEST_EQ(jacobian_adjugate_components(QUAD4), 4);
    SMESH_TEST_EQ(fff_components(QUAD4), 3);

    auto jac = JacobianAdjugateAndDeterminant::create_SoA(mesh, MEMORY_SPACE_HOST, 0);
    SMESH_TEST_ASSERT(jac != nullptr);
    SMESH_TEST_EQ(static_cast<int>(jac->jacobian_adjugate_SoA()->extent(0)), 4);
    SMESH_TEST_APPROXEQ(jac->jacobian_determinant()->data()[0], static_cast<geom_t>(1), kTol);
    jacobian_t **const adj = jac->jacobian_adjugate_SoA()->data();
    SMESH_TEST_APPROXEQ(adj[0][0], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(adj[1][0], static_cast<jacobian_t>(0), kTol);
    SMESH_TEST_APPROXEQ(adj[2][0], static_cast<jacobian_t>(0), kTol);
    SMESH_TEST_APPROXEQ(adj[3][0], static_cast<jacobian_t>(1), kTol);

    auto fff = FFF::create_SoA(mesh, MEMORY_SPACE_HOST, 0);
    SMESH_TEST_ASSERT(fff != nullptr);
    SMESH_TEST_EQ(static_cast<int>(fff->fff_SoA()->extent(0)), 3);
    jacobian_t **const f = fff->fff_SoA()->data();
    SMESH_TEST_APPROXEQ(f[0][0], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(f[1][0], static_cast<jacobian_t>(0), kTol);
    SMESH_TEST_APPROXEQ(f[2][0], static_cast<jacobian_t>(1), kTol);
    return SMESH_TEST_SUCCESS;
}

static int test_tri3_jacobian_fff_analytic() {
    auto mesh = Mesh::create_tri3_square(Communicator::self(), 1, 1);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(mesh->spatial_dimension(), 2);
    SMESH_TEST_EQ(jacobian_adjugate_components(TRI3), 4);
    SMESH_TEST_EQ(fff_components(TRI3), 3);
    SMESH_TEST_EQ(mesh->n_elements(0), static_cast<ptrdiff_t>(2));

    auto jac = JacobianAdjugateAndDeterminant::create_SoA(mesh, MEMORY_SPACE_HOST, 0);
    SMESH_TEST_ASSERT(jac != nullptr);
    SMESH_TEST_EQ(static_cast<int>(jac->jacobian_adjugate_SoA()->extent(0)), 4);
    const geom_t *const det = jac->jacobian_determinant()->data();
    jacobian_t **const  adj = jac->jacobian_adjugate_SoA()->data();
    SMESH_TEST_APPROXEQ(det[0], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(adj[0][0], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(adj[1][0], static_cast<jacobian_t>(0), kTol);
    SMESH_TEST_APPROXEQ(adj[2][0], static_cast<jacobian_t>(0), kTol);
    SMESH_TEST_APPROXEQ(adj[3][0], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(det[1], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(adj[0][1], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(adj[1][1], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(adj[2][1], static_cast<jacobian_t>(-1), kTol);
    SMESH_TEST_APPROXEQ(adj[3][1], static_cast<jacobian_t>(0), kTol);

    auto fff = FFF::create_SoA(mesh, MEMORY_SPACE_HOST, 0);
    SMESH_TEST_ASSERT(fff != nullptr);
    jacobian_t **const f = fff->fff_SoA()->data();
    SMESH_TEST_APPROXEQ(f[0][0], static_cast<jacobian_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(f[1][0], static_cast<jacobian_t>(0), kTol);
    SMESH_TEST_APPROXEQ(f[2][0], static_cast<jacobian_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(f[0][1], static_cast<jacobian_t>(1), kTol);
    SMESH_TEST_APPROXEQ(f[1][1], static_cast<jacobian_t>(-0.5), kTol);
    SMESH_TEST_APPROXEQ(f[2][1], static_cast<jacobian_t>(0.5), kTol);
    return SMESH_TEST_SUCCESS;
}

static int test_sfc_unpadded_quad_square() {
    auto mesh = Mesh::create_quad4_square(Communicator::self(), 4, 3);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(mesh->spatial_dimension(), 2);
    const ptrdiff_t n_e = mesh->n_elements();
    const ptrdiff_t n_n = mesh->n_nodes();
    SFC             sfc("morton3");
    SMESH_TEST_EQ(sfc.reorder(*mesh), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh->n_elements(), n_e);
    SMESH_TEST_EQ(mesh->n_nodes(), n_n);
    SMESH_TEST_EQ(mesh->spatial_dimension(), 2);
    return SMESH_TEST_SUCCESS;
}

static int test_serial_unpadded_extrude() {
    auto quad = Mesh::create_quad4_square(Communicator::self(), 2, 2);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_EQ(quad->spatial_dimension(), 2);
    auto hex = extrude(quad, 1.0, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    SMESH_TEST_EQ(hex->element_type(0), HEX8);
    SMESH_TEST_EQ(hex->spatial_dimension(), 3);
    SMESH_TEST_EQ(hex->n_nodes(), quad->n_nodes() * 3);
    SMESH_TEST_EQ(hex->n_elements(), quad->n_elements() * 2);

    auto tri = Mesh::create_tri3_square(Communicator::self(), 2, 2);
    auto wedge = extrude(tri, 1.0, 2);
    SMESH_TEST_ASSERT(wedge != nullptr);
    SMESH_TEST_EQ(wedge->element_type(0), WEDGE6);
    SMESH_TEST_EQ(wedge->n_nodes(), tri->n_nodes() * 3);
    SMESH_TEST_EQ(wedge->n_elements(), tri->n_elements() * 2);
    return SMESH_TEST_SUCCESS;
}

static int test_circle_quad_square_boundary() {
    auto mesh = Mesh::create_quad4_square(Communicator::self(), 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss = skin_sideset(mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    auto ids = create_nodeset_from_sideset(mesh, ss);
    SMESH_TEST_ASSERT(ids != nullptr);
    auto ns = Nodeset::create(mesh->comm(), ids);
    SMESH_TEST_ASSERT(ns != nullptr);
    SMESH_TEST_ASSERT(ns->size() > 0);

    const ptrdiff_t n = mesh->n_nodes();
    auto            before = create_host_buffer<geom_t>(2, static_cast<size_t>(n));
    std::memcpy(before->data()[0], mesh->points()->data()[0], static_cast<size_t>(n) * sizeof(geom_t));
    std::memcpy(before->data()[1], mesh->points()->data()[1], static_cast<size_t>(n) * sizeof(geom_t));

    std::vector<char> on_skin(static_cast<size_t>(n), 0);
    const idx_t *const nids = ns->nodes()->data();
    for (ptrdiff_t i = 0; i < ns->size(); ++i) {
        on_skin[static_cast<size_t>(nids[i])] = 1;
    }

    const geom_t cx = static_cast<geom_t>(0.5);
    const geom_t cy = static_cast<geom_t>(0.5);
    const geom_t r  = static_cast<geom_t>(0.5);
    auto         p  = CircleParametrization::create(ns, cx, cy, 99, 1, 0, 0, r);
    SMESH_TEST_ASSERT(p != nullptr);
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);

    geom_t **const pts = mesh->points()->data();
    geom_t **const old = before->data();
    for (ptrdiff_t i = 0; i < n; ++i) {
        if (!on_skin[static_cast<size_t>(i)]) {
            SMESH_TEST_APPROXEQ(pts[0][i], old[0][i], kTol);
            SMESH_TEST_APPROXEQ(pts[1][i], old[1][i], kTol);
            continue;
        }
        const geom_t dx = pts[0][i] - cx;
        const geom_t dy = pts[1][i] - cy;
        SMESH_TEST_APPROXEQ(std::sqrt(dx * dx + dy * dy), r, kTol);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_promote_quad4_to_quad9() {
    auto q4 = Mesh::create_quad4_square(Communicator::self(), 1, 1);
    SMESH_TEST_ASSERT(q4 != nullptr);
    auto q9 = promote_to(QUAD9, q4);
    SMESH_TEST_ASSERT(q9 != nullptr);
    SMESH_TEST_EQ(q9->element_type(0), QUAD9);
    SMESH_TEST_EQ(q9->n_elements(), static_cast<ptrdiff_t>(1));
    SMESH_TEST_ASSERT(q9->n_nodes() >= 9);
    SMESH_TEST_EQ(q9->n_nodes_per_element(0), 9);

    idx_t **const  e = q9->elements(0)->data();
    geom_t **const p = q9->points()->data();
    SMESH_TEST_APPROXEQ(p[0][e[0][0]], static_cast<geom_t>(0), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[0][0]], static_cast<geom_t>(0), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[1][0]], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[1][0]], static_cast<geom_t>(0), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[2][0]], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[2][0]], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[3][0]], static_cast<geom_t>(0), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[3][0]], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[4][0]], static_cast<geom_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[4][0]], static_cast<geom_t>(0), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[5][0]], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[5][0]], static_cast<geom_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[6][0]], static_cast<geom_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[6][0]], static_cast<geom_t>(1), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[7][0]], static_cast<geom_t>(0), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[7][0]], static_cast<geom_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(p[0][e[8][0]], static_cast<geom_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(p[1][e[8][0]], static_cast<geom_t>(0.5), kTol);

    auto q4b = Mesh::create_quad4_square(Communicator::self(), 2, 2);
    auto q9b = promote_to(QUAD9, q4b);
    SMESH_TEST_ASSERT(q9b != nullptr);
    SMESH_TEST_ASSERT(q9b->n_nodes() > q4b->n_nodes());
    std::vector<char> seen(static_cast<size_t>(q9b->n_nodes()), 0);
    idx_t **const     e2 = q9b->elements(0)->data();
    for (ptrdiff_t el = 0; el < q9b->n_elements(); ++el) {
        const idx_t c = e2[8][el];
        SMESH_TEST_ASSERT(c >= 0 && static_cast<ptrdiff_t>(c) < q9b->n_nodes());
        SMESH_TEST_ASSERT(!seen[static_cast<size_t>(c)]);
        seen[static_cast<size_t>(c)] = 1;
    }
    return SMESH_TEST_SUCCESS;
}

static int test_promote_trishell3_to_trishell6() {
    auto t3 = Mesh::create_tri3_square(Communicator::self(), 1, 1);
    SMESH_TEST_ASSERT(t3 != nullptr);
    t3->set_element_type(0, TRISHELL3);
    auto t6 = promote_to(TRISHELL6, t3);
    SMESH_TEST_ASSERT(t6 != nullptr);
    SMESH_TEST_EQ(t6->element_type(0), TRISHELL6);
    SMESH_TEST_EQ(t6->n_nodes_per_element(0), 6);
    SMESH_TEST_EQ(t6->n_elements(), t3->n_elements());
    SMESH_TEST_EQ(t6->n_nodes(), static_cast<ptrdiff_t>(9));

    idx_t **const  e = t6->elements(0)->data();
    geom_t **const p = t6->points()->data();
    for (ptrdiff_t el = 0; el < t6->n_elements(); ++el) {
        for (int m = 0; m < 3; ++m) {
            const int a = m;
            const int b = (m + 1) % 3;
            const idx_t mid = e[3 + m][el];
            SMESH_TEST_APPROXEQ(p[0][mid], (p[0][e[a][el]] + p[0][e[b][el]]) / geom_t(2), kTol);
            SMESH_TEST_APPROXEQ(p[1][mid], (p[1][e[a][el]] + p[1][e[b][el]]) / geom_t(2), kTol);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_ssquad_square_factory() {
    auto quad = Mesh::create_quad4_square(Communicator::self(), 2, 3);
    SMESH_TEST_ASSERT(quad != nullptr);
    auto via_ss = to_semistructured(2, quad, false, false);
    SMESH_TEST_ASSERT(via_ss != nullptr);
    auto factory = Mesh::create_semistructured_quad_square(Communicator::self(), 2, 2, 3);
    SMESH_TEST_ASSERT(factory != nullptr);
    SMESH_TEST_EQ(factory->element_type(0), proteus_quad_type(2));
    SMESH_TEST_EQ(factory->element_type(0), via_ss->element_type(0));
    SMESH_TEST_EQ(factory->n_nodes(), via_ss->n_nodes());
    SMESH_TEST_EQ(factory->n_elements(), via_ss->n_elements());
    SMESH_TEST_EQ(factory->spatial_dimension(), 2);

    auto from_square = Mesh::create_square(Communicator::self(), PROTEUS_QUAD9, 2, 3);
    SMESH_TEST_ASSERT(from_square != nullptr);
    SMESH_TEST_EQ(from_square->element_type(0), proteus_quad_type(2));
    SMESH_TEST_EQ(from_square->n_nodes(), factory->n_nodes());
    SMESH_TEST_EQ(from_square->n_elements(), factory->n_elements());
    return SMESH_TEST_SUCCESS;
}

static int test_serial_trishell3_extrude() {
    auto tri = Mesh::create_tri3_square(Communicator::self(), 2, 2);
    SMESH_TEST_ASSERT(tri != nullptr);
    tri->set_element_type(0, TRISHELL3);
    auto wedge = extrude(tri, 1.0, 2);
    SMESH_TEST_ASSERT(wedge != nullptr);
    SMESH_TEST_EQ(wedge->element_type(0), WEDGE6);
    SMESH_TEST_EQ(wedge->n_elements(), tri->n_elements() * 2);
    SMESH_TEST_EQ(wedge->n_nodes(), tri->n_nodes() * 3);
    return SMESH_TEST_SUCCESS;
}

static int test_mpi_unpadded_extrude() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;

    auto q_serial = Mesh::create_quad4_square(Communicator::self(), nx, ny);
    SMESH_TEST_EQ(q_serial->spatial_dimension(), 2);
    auto q_ex_s = extrude(q_serial, 1.0, 2);
    SMESH_TEST_ASSERT(q_ex_s != nullptr);
    auto q_par = Mesh::create_quad4_square(comm, nx, ny);
    SMESH_TEST_ASSERT(q_par != nullptr);
    SMESH_TEST_EQ(q_par->spatial_dimension(), 2);
    auto q_ex = extrude(q_par, 1.0, 2);
    SMESH_TEST_ASSERT(q_ex != nullptr);
    SMESH_TEST_EQ(q_ex->element_type(0), HEX8);
    SMESH_TEST_EQ(q_ex->distributed()->n_nodes_global(), q_ex_s->n_nodes());
    SMESH_TEST_EQ(q_ex->distributed()->n_elements_global(), q_ex_s->n_elements());

    auto t_serial = Mesh::create_tri3_square(Communicator::self(), nx, ny);
    auto t_ex_s   = extrude(t_serial, 1.0, 2);
    auto t_par    = Mesh::create_tri3_square(comm, nx, ny);
    auto t_ex     = extrude(t_par, 1.0, 2);
    SMESH_TEST_ASSERT(t_ex != nullptr);
    SMESH_TEST_EQ(t_ex->element_type(0), WEDGE6);
    SMESH_TEST_EQ(t_ex->distributed()->n_nodes_global(), t_ex_s->n_nodes());
    SMESH_TEST_EQ(t_ex->distributed()->n_elements_global(), t_ex_s->n_elements());
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_trishell3_extrude() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;

    auto serial = Mesh::create_tri3_square(Communicator::self(), nx, ny);
    serial->set_element_type(0, TRISHELL3);
    auto serial_ex = extrude(serial, 1.0, 2);
    SMESH_TEST_ASSERT(serial_ex != nullptr);
    SMESH_TEST_EQ(serial_ex->element_type(0), WEDGE6);

    auto par = Mesh::create_tri3_square(comm, nx, ny);
    SMESH_TEST_ASSERT(par != nullptr);
    par->set_element_type(0, TRISHELL3);
    auto par_ex = extrude(par, 1.0, 2);
    SMESH_TEST_ASSERT(par_ex != nullptr);
    SMESH_TEST_EQ(par_ex->element_type(0), WEDGE6);
    SMESH_TEST_EQ(par_ex->distributed()->n_nodes_global(), serial_ex->n_nodes());
    SMESH_TEST_EQ(par_ex->distributed()->n_elements_global(), serial_ex->n_elements());
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_ssquad_square_factory() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    auto            serial = Mesh::create_semistructured_quad_square(Communicator::self(), 2, nx, ny);
    auto            par    = Mesh::create_semistructured_quad_square(comm, 2, nx, ny);
    SMESH_TEST_ASSERT(serial != nullptr);
    SMESH_TEST_ASSERT(par != nullptr);
    SMESH_TEST_EQ(par->element_type(0), proteus_quad_type(2));
    SMESH_TEST_EQ(par->distributed()->n_nodes_global(), serial->n_nodes());
    SMESH_TEST_EQ(par->distributed()->n_elements_global(), serial->n_elements());
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_sfc_checkerboard_cube() {
    auto mesh = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
    return check_multiblock_sfc(*mesh, "morton3");
}

static int test_sfc_checkerboard_cube_hilbert() {
    auto mesh = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    return check_multiblock_sfc(*mesh, "hilbert3");
}

static int test_sfc_hex8_tet4_cube() {
    auto mesh = Mesh::create_hex8_tet4_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
    SMESH_TEST_EQ(mesh->element_type(0), HEX8);
    SMESH_TEST_EQ(mesh->element_type(1), TET4);
    return check_multiblock_sfc(*mesh, "morton3");
}

static int test_sfc_quad4_split_block() {
    auto quad = Mesh::create_quad4_square(Communicator::self(), 4, 3);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_EQ(quad->spatial_dimension(), 2);
    auto mesh = split_first_half(quad);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
    SMESH_TEST_EQ(mesh->spatial_dimension(), 2);
    return check_multiblock_sfc(*mesh, "morton3");
}

static int test_mpi_sfc_hex8_cube() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t n = std::max<ptrdiff_t>(2 * comm->size(), 4);
    auto mesh = Mesh::create_hex8_cube(comm, n, 4, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 1);
    return check_distributed_sfc(*mesh, "morton3");
#endif
}

static int test_mpi_sfc_quad4_square() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    auto mesh = Mesh::create_quad4_square(comm, nx, 3);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());
    SMESH_TEST_EQ(mesh->spatial_dimension(), 2);
    return check_distributed_sfc(*mesh, "morton3");
#endif
}

static int test_mpi_sfc_checkerboard_cube() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t n = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t even = n + (n & 1);
    auto mesh = Mesh::create_hex8_checkerboard_cube(comm, even, even, even);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
    return check_distributed_sfc(*mesh, "hilbert3");
#endif
}

static int test_mpi_sfc_hex8_tet4_cube() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t n = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t even = n + (n & 1);
    auto mesh = Mesh::create_hex8_tet4_cube(comm, even, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
    return check_distributed_sfc(*mesh, "morton3");
#endif
}

static int test_mpi_sfc_checkerboard_sideset_remap() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const ptrdiff_t n = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t even = n + (n & 1);
    auto mesh = Mesh::create_hex8_checkerboard_cube(comm, even, even, even);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto sidesets = Sideset::create_from_selector(
        mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_ASSERT(!sidesets.empty());
    for (size_t i = 0; i < sidesets.size(); ++i) {
        mesh->add_sideset("left", sidesets[i]);
    }
    SMESH_TEST_EQ(SFC("morton3").reorder(*mesh), SMESH_SUCCESS);
    auto recreated = Sideset::create_from_selector(
        mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_EQ(recreated.size(), sidesets.size());
    const auto registered = mesh->sidesets("left");
    SMESH_TEST_EQ(registered.size(), sidesets.size());
    for (size_t i = 0; i < registered.size(); ++i) {
        SMESH_TEST_EQ(registered[i]->size(), recreated[i]->size());
        SMESH_TEST_EQ(registered[i]->block_id(), recreated[i]->block_id());
        std::vector<std::pair<element_idx_t, i16>> a(static_cast<size_t>(registered[i]->size()));
        std::vector<std::pair<element_idx_t, i16>> b(static_cast<size_t>(recreated[i]->size()));
        for (ptrdiff_t s = 0; s < registered[i]->size(); ++s) {
            a[static_cast<size_t>(s)] = {registered[i]->parent()->data()[s],
                                         registered[i]->lfi()->data()[s]};
            b[static_cast<size_t>(s)] = {recreated[i]->parent()->data()[s],
                                         recreated[i]->lfi()->data()[s]};
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        SMESH_TEST_ASSERT(a == b);
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_quad4_jacobian_fff_analytic);
    SMESH_RUN_TEST(test_tri3_jacobian_fff_analytic);
    SMESH_RUN_TEST(test_sfc_unpadded_quad_square);
    SMESH_RUN_TEST(test_sfc_checkerboard_cube);
    SMESH_RUN_TEST(test_sfc_checkerboard_cube_hilbert);
    SMESH_RUN_TEST(test_sfc_hex8_tet4_cube);
    SMESH_RUN_TEST(test_sfc_quad4_split_block);
    SMESH_RUN_TEST(test_serial_unpadded_extrude);
    SMESH_RUN_TEST(test_circle_quad_square_boundary);
    SMESH_RUN_TEST(test_promote_quad4_to_quad9);
    SMESH_RUN_TEST(test_promote_trishell3_to_trishell6);
    SMESH_RUN_TEST(test_ssquad_square_factory);
    SMESH_RUN_TEST(test_serial_trishell3_extrude);
    SMESH_RUN_TEST(test_mpi_unpadded_extrude);
    SMESH_RUN_TEST(test_mpi_trishell3_extrude);
    SMESH_RUN_TEST(test_mpi_ssquad_square_factory);
    SMESH_RUN_TEST(test_mpi_sfc_hex8_cube);
    SMESH_RUN_TEST(test_mpi_sfc_quad4_square);
    SMESH_RUN_TEST(test_mpi_sfc_checkerboard_cube);
    SMESH_RUN_TEST(test_mpi_sfc_hex8_tet4_cube);
    SMESH_RUN_TEST(test_mpi_sfc_checkerboard_sideset_remap);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
