#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "smesh_edgeset.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_parametrization.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#endif

using namespace smesh;

namespace {

constexpr geom_t kTol = sizeof(geom_t) == 8 ? static_cast<geom_t>(1e-10)
                                            : static_cast<geom_t>(1e-5);

std::shared_ptr<Nodeset> nodeset_from_sideset(const std::shared_ptr<Mesh>    &mesh,
                                              const std::shared_ptr<Sideset> &ss) {
    auto ids = create_nodeset_from_sideset(mesh, ss);
    if (!ids) {
        return nullptr;
    }
    return Nodeset::create(mesh->comm(), ids);
}

std::shared_ptr<Nodeset> unique_block_nodes(const std::shared_ptr<Mesh> &mesh, block_idx_t block_id) {
    auto block = mesh->block(block_id);
    if (!block || !block->elements()) {
        return nullptr;
    }
    const int       nxe     = block->n_nodes_per_element();
    const ptrdiff_t ne      = block->n_elements();
    const ptrdiff_t n_nodes = mesh->n_nodes();
    idx_t **const   elems   = block->elements()->data();
    std::vector<char> mark(static_cast<size_t>(n_nodes), 0);
    ptrdiff_t         count = 0;
    for (int d = 0; d < nxe; ++d) {
        for (ptrdiff_t e = 0; e < ne; ++e) {
            const idx_t id = elems[d][e];
            if (id < 0 || static_cast<ptrdiff_t>(id) >= n_nodes) {
                return nullptr;
            }
            if (!mark[static_cast<size_t>(id)]) {
                mark[static_cast<size_t>(id)] = 1;
                ++count;
            }
        }
    }
    auto buf = create_host_buffer<idx_t>(static_cast<size_t>(count));
    ptrdiff_t k = 0;
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        if (mark[static_cast<size_t>(i)]) {
            buf->data()[k++] = static_cast<idx_t>(i);
        }
    }
    return Nodeset::create(mesh->comm(), buf);
}

SharedBuffer<geom_t *> copy_points(const std::shared_ptr<Mesh> &mesh) {
    const int       sdim = mesh->spatial_dimension();
    const ptrdiff_t n    = mesh->n_nodes();
    auto            out  = create_host_buffer<geom_t>(static_cast<size_t>(sdim), static_cast<size_t>(n));
    geom_t **const  s    = mesh->points()->data();
    geom_t **const  d    = out->data();
    for (int c = 0; c < sdim; ++c) {
        std::memcpy(d[c], s[c], static_cast<size_t>(n) * sizeof(geom_t));
    }
    return out;
}

int mark_nodeset(const Nodeset &ns, ptrdiff_t n_nodes, std::vector<char> *mark) {
    mark->assign(static_cast<size_t>(n_nodes), 0);
    const idx_t *ids = ns.size() > 0 ? ns.nodes()->data() : nullptr;
    for (ptrdiff_t i = 0; i < ns.size(); ++i) {
        const idx_t id = ids[i];
        if (id < 0 || static_cast<ptrdiff_t>(id) >= n_nodes) {
            return SMESH_TEST_FAILURE;
        }
        (*mark)[static_cast<size_t>(id)] = 1;
    }
    return SMESH_TEST_SUCCESS;
}

}  // namespace

static int test_sphere_skin_apply() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss = skin_sideset(mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    auto ns = nodeset_from_sideset(mesh, ss);
    SMESH_TEST_ASSERT(ns != nullptr);
    SMESH_TEST_ASSERT(ns->size() > 0);

    auto before = copy_points(mesh);
    std::vector<char> on_skin;
    SMESH_TEST_EQ(mark_nodeset(*ns, mesh->n_nodes(), &on_skin), SMESH_TEST_SUCCESS);

    const geom_t cx = static_cast<geom_t>(0.5);
    const geom_t cy = static_cast<geom_t>(0.5);
    const geom_t cz = static_cast<geom_t>(0.5);
    const geom_t r  = static_cast<geom_t>(0.5);
    auto         p  = SphereParametrization::create(ns, cx, cy, cz, r);
    SMESH_TEST_ASSERT(p != nullptr);
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);

    geom_t **const pts = mesh->points()->data();
    geom_t **const old = before->data();
    const ptrdiff_t n  = mesh->n_nodes();
    for (ptrdiff_t i = 0; i < n; ++i) {
        if (!on_skin[static_cast<size_t>(i)]) {
            SMESH_TEST_APPROXEQ(pts[0][i], old[0][i], kTol);
            SMESH_TEST_APPROXEQ(pts[1][i], old[1][i], kTol);
            SMESH_TEST_APPROXEQ(pts[2][i], old[2][i], kTol);
            continue;
        }
        const geom_t dx = pts[0][i] - cx;
        const geom_t dy = pts[1][i] - cy;
        const geom_t dz = pts[2][i] - cz;
        const geom_t rr = std::sqrt(dx * dx + dy * dy + dz * dz);
        SMESH_TEST_APPROXEQ(rr, r, kTol);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_sphere_center_fallback() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const ptrdiff_t n = mesh->n_nodes();
    geom_t **const  pts = mesh->points()->data();
    idx_t           center_id = -1;
    for (ptrdiff_t i = 0; i < n; ++i) {
        if (std::fabs(pts[0][i] - static_cast<geom_t>(0.5)) < kTol &&
            std::fabs(pts[1][i] - static_cast<geom_t>(0.5)) < kTol &&
            std::fabs(pts[2][i] - static_cast<geom_t>(0.5)) < kTol) {
            center_id = static_cast<idx_t>(i);
            break;
        }
    }
    SMESH_TEST_ASSERT(center_id >= 0);
    auto ids = create_host_buffer<idx_t>(1);
    ids->data()[0] = center_id;
    auto ns        = Nodeset::create(mesh->comm(), ids);
    auto p         = SphereParametrization::create(ns,
                                           static_cast<geom_t>(0.5),
                                           static_cast<geom_t>(0.5),
                                           static_cast<geom_t>(0.5),
                                           static_cast<geom_t>(0.5));
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);
    SMESH_TEST_APPROXEQ(pts[0][center_id], static_cast<geom_t>(1.0), kTol);
    SMESH_TEST_APPROXEQ(pts[1][center_id], static_cast<geom_t>(0.5), kTol);
    SMESH_TEST_APPROXEQ(pts[2][center_id], static_cast<geom_t>(0.5), kTol);
    return SMESH_TEST_SUCCESS;
}

static int test_circle_edgeset_apply() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 1, 1, 1);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const int n_edges = elem_num_edges(HEX8);
    SMESH_TEST_EQ(n_edges, 12);
    auto parent = create_host_buffer<element_idx_t>(static_cast<size_t>(n_edges));
    auto lei    = create_host_buffer<i16>(static_cast<size_t>(n_edges));
    for (int e = 0; e < n_edges; ++e) {
        parent->data()[e] = 0;
        lei->data()[e]    = static_cast<i16>(e);
    }
    auto es = Edgeset::create(mesh->comm(), parent, lei, 0);
    SMESH_TEST_ASSERT(es != nullptr);
    auto ns = create_nodeset_from_edgeset(mesh, es);
    SMESH_TEST_ASSERT(ns != nullptr);
    SMESH_TEST_ASSERT(ns->size() > 0);

    const geom_t cx = static_cast<geom_t>(0.5);
    const geom_t cy = static_cast<geom_t>(0.5);
    const geom_t cz = static_cast<geom_t>(0.5);
    const geom_t r  = static_cast<geom_t>(0.5);
    auto         p  = CircleParametrization::create(ns, cx, cy, cz, 0, 0, 1, r);
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);

    geom_t **const pts = mesh->points()->data();
    const idx_t   *ids = ns->nodes()->data();
    for (ptrdiff_t i = 0; i < ns->size(); ++i) {
        const idx_t  id = ids[i];
        const geom_t dx = pts[0][id] - cx;
        const geom_t dy = pts[1][id] - cy;
        const geom_t dz = pts[2][id] - cz;
        SMESH_TEST_APPROXEQ(dz, static_cast<geom_t>(0), kTol);
        const geom_t rr = std::sqrt(dx * dx + dy * dy);
        SMESH_TEST_APPROXEQ(rr, r, kTol);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_identity_block_nodes_noop() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 2, 1, 1);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ns = unique_block_nodes(mesh, 0);
    SMESH_TEST_ASSERT(ns != nullptr);
    SMESH_TEST_EQ(ns->size(), mesh->n_nodes());
    auto before = copy_points(mesh);
    auto p      = IdentityParametrization::create(ns);
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);
    geom_t **const pts = mesh->points()->data();
    geom_t **const old = before->data();
    const ptrdiff_t n  = mesh->n_nodes();
    for (ptrdiff_t i = 0; i < n; ++i) {
        SMESH_TEST_APPROXEQ(pts[0][i], old[0][i], kTol);
        SMESH_TEST_APPROXEQ(pts[1][i], old[1][i], kTol);
        SMESH_TEST_APPROXEQ(pts[2][i], old[2][i], kTol);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_polynomial_surface_apply() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 1);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ns = unique_block_nodes(mesh, 0);
    SMESH_TEST_ASSERT(ns != nullptr);

    // ζ = 0.1 (ξ² + η²). Monomial order: i outer, j inner, i+j <= 2.
    const geom_t coeffs[6] = {
            0, 0, static_cast<geom_t>(0.1), 0, 0, static_cast<geom_t>(0.1)};
    auto p = PolynomialSurfaceParametrization::create(
            ns, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, coeffs, 6);
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);

    geom_t **const pts = mesh->points()->data();
    const idx_t   *ids = ns->nodes()->data();
    for (ptrdiff_t i = 0; i < ns->size(); ++i) {
        const idx_t  id = ids[i];
        const geom_t x  = pts[0][id];
        const geom_t y  = pts[1][id];
        const geom_t z  = pts[2][id];
        const geom_t expect = static_cast<geom_t>(0.1) * (x * x + y * y);
        SMESH_TEST_APPROXEQ(z, expect, kTol);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_mesh_parametrization_registry_clone() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 1, 1, 1);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss = skin_sideset(mesh);
    auto ns = nodeset_from_sideset(mesh, ss);
    auto p  = SphereParametrization::create(
            ns, static_cast<geom_t>(0.5), static_cast<geom_t>(0.5), static_cast<geom_t>(0.5), static_cast<geom_t>(0.5));
    mesh->add_parametrization("sphere", p);
    SMESH_TEST_EQ(mesh->parametrizations().size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(mesh->parametrizations("sphere").size(), static_cast<size_t>(1));
    SMESH_TEST_ASSERT(mesh->parametrizations("sphere")[0].get() == p.get());

    auto cloned = mesh->clone();
    SMESH_TEST_ASSERT(cloned != nullptr);
    SMESH_TEST_EQ(cloned->parametrizations("sphere").size(), static_cast<size_t>(1));
    SMESH_TEST_ASSERT(cloned->parametrizations("sphere")[0].get() == p.get());
    SMESH_TEST_EQ(cloned->parametrizations("sphere")[0]->apply(*cloned), SMESH_SUCCESS);

    geom_t **const pts = cloned->points()->data();
    const idx_t   *ids = ns->nodes()->data();
    for (ptrdiff_t i = 0; i < ns->size(); ++i) {
        const idx_t  id = ids[i];
        const geom_t dx = pts[0][id] - static_cast<geom_t>(0.5);
        const geom_t dy = pts[1][id] - static_cast<geom_t>(0.5);
        const geom_t dz = pts[2][id] - static_cast<geom_t>(0.5);
        const geom_t rr = std::sqrt(dx * dx + dy * dy + dz * dz);
        SMESH_TEST_APPROXEQ(rr, static_cast<geom_t>(0.5), kTol);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_mpi_sphere_skin_gid_parity() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    const geom_t    cx = static_cast<geom_t>(0.5);
    const geom_t    cy = static_cast<geom_t>(0.5);
    const geom_t    cz = static_cast<geom_t>(0.5);
    const geom_t    r  = static_cast<geom_t>(0.5);

    i64 n_serial = 0;
    SharedBuffer<geom_t *> serial_pts;
    if (comm->rank() == 0) {
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
        SMESH_TEST_ASSERT(serial != nullptr);
        auto ss = skin_sideset(serial);
        SMESH_TEST_ASSERT(ss != nullptr);
        auto ns = nodeset_from_sideset(serial, ss);
        SMESH_TEST_ASSERT(ns != nullptr);
        auto p = SphereParametrization::create(ns, cx, cy, cz, r);
        SMESH_TEST_EQ(p->apply(*serial), SMESH_SUCCESS);
        n_serial   = static_cast<i64>(serial->n_nodes());
        serial_pts = copy_points(serial);
    }
    comm->broadcast(&n_serial, 1, 0);
    SMESH_TEST_ASSERT(n_serial > 0);
    if (comm->rank() != 0) {
        serial_pts = create_host_buffer<geom_t>(3, static_cast<size_t>(n_serial));
    }
    comm->broadcast(serial_pts->data()[0], static_cast<int>(n_serial), 0);
    comm->broadcast(serial_pts->data()[1], static_cast<int>(n_serial), 0);
    comm->broadcast(serial_pts->data()[2], static_cast<int>(n_serial), 0);

    auto mesh = Mesh::create_hex8_cube(comm, nx, ny, nz);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());
    auto ss = skin_sideset(mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    auto ns = nodeset_from_sideset(mesh, ss);
    SMESH_TEST_ASSERT(ns != nullptr);
    auto p = SphereParametrization::create(ns, cx, cy, cz, r);
    SMESH_TEST_EQ(p->apply(*mesh), SMESH_SUCCESS);

    geom_t **const        pts = mesh->points()->data();
    const idx_t          *ids = ns->size() > 0 ? ns->nodes()->data() : nullptr;
    const large_idx_t    *map = mesh->distributed()->node_mapping()->data();
    geom_t **const        ref = serial_pts->data();
    for (ptrdiff_t i = 0; i < ns->size(); ++i) {
        const idx_t id = ids[i];
        SMESH_TEST_ASSERT(id >= 0 && static_cast<ptrdiff_t>(id) < mesh->n_nodes());
        const large_idx_t gid = map[id];
        SMESH_TEST_ASSERT(gid >= 0 && gid < n_serial);
        SMESH_TEST_APPROXEQ(pts[0][id], ref[0][gid], kTol);
        SMESH_TEST_APPROXEQ(pts[1][id], ref[1][gid], kTol);
        SMESH_TEST_APPROXEQ(pts[2][id], ref[2][gid], kTol);
        const geom_t dx = pts[0][id] - cx;
        const geom_t dy = pts[1][id] - cy;
        const geom_t dz = pts[2][id] - cz;
        const geom_t rr = std::sqrt(dx * dx + dy * dy + dz * dz);
        SMESH_TEST_APPROXEQ(rr, r, kTol);
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_sphere_skin_apply);
    SMESH_RUN_TEST(test_sphere_center_fallback);
    SMESH_RUN_TEST(test_circle_edgeset_apply);
    SMESH_RUN_TEST(test_identity_block_nodes_noop);
    SMESH_RUN_TEST(test_polynomial_surface_apply);
    SMESH_RUN_TEST(test_mesh_parametrization_registry_clone);
    SMESH_RUN_TEST(test_mpi_sphere_skin_gid_parity);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
