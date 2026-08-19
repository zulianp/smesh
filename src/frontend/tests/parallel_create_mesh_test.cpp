#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <vector>

#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static ptrdiff_t even_at_least(const ptrdiff_t v) {
    const ptrdiff_t a = std::max<ptrdiff_t>(v, 2);
    return a + (a & 1);
}

static int check_unique_node_ids(const Mesh &mesh) {
    const ptrdiff_t n_nodes = mesh.n_nodes();
    SMESH_TEST_ASSERT(n_nodes > 0);
    for (size_t b = 0; b < mesh.n_blocks(); ++b) {
        auto            block = mesh.block(b);
        const ptrdiff_t n_e   = block->n_elements();
        const int       nxe   = block->n_nodes_per_element();
        idx_t         **elems = block->elements()->data();
        for (int v = 0; v < nxe; ++v) {
            for (ptrdiff_t e = 0; e < n_e; ++e) {
                SMESH_TEST_ASSERT(elems[v][e] >= 0);
                SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(elems[v][e]) < n_nodes);
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_create_size_one_serial() {
    auto comm = Communicator::self();
    SMESH_TEST_EQ(comm->size(), 1);

    auto hex = Mesh::create_hex8_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    SMESH_TEST_ASSERT(!hex->is_distributed());
    SMESH_TEST_EQ(hex->n_nodes(), static_cast<ptrdiff_t>(3 * 3 * 3));
    SMESH_TEST_EQ(hex->n_elements(), static_cast<ptrdiff_t>(8));
    SMESH_TEST_EQ(static_cast<int>(hex->n_blocks()), 1);
    SMESH_TEST_EQ(hex->element_type(0), HEX8);

    auto tet = Mesh::create_tet4_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_ASSERT(!tet->is_distributed());
    SMESH_TEST_EQ(tet->n_nodes(), static_cast<ptrdiff_t>(27 + 8));
    SMESH_TEST_EQ(tet->n_elements(), static_cast<ptrdiff_t>(12 * 8));
    SMESH_TEST_EQ(tet->element_type(0), TET4);

    auto quad = Mesh::create_quad4_square(comm, 4, 3);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_ASSERT(!quad->is_distributed());
    SMESH_TEST_EQ(quad->n_nodes(), static_cast<ptrdiff_t>(5 * 4));
    SMESH_TEST_EQ(quad->n_elements(), static_cast<ptrdiff_t>(12));
    SMESH_TEST_EQ(quad->spatial_dimension(), 2);
    SMESH_TEST_EQ(quad->element_type(0), QUAD4);

    auto tri = Mesh::create_tri3_square(comm, 4, 3);
    SMESH_TEST_ASSERT(tri != nullptr);
    SMESH_TEST_ASSERT(!tri->is_distributed());
    SMESH_TEST_EQ(tri->n_nodes(), static_cast<ptrdiff_t>(5 * 4));
    SMESH_TEST_EQ(tri->n_elements(), static_cast<ptrdiff_t>(24));
    SMESH_TEST_EQ(tri->spatial_dimension(), 2);
    SMESH_TEST_EQ(tri->element_type(0), TRI3);

    auto cb = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(cb != nullptr);
    SMESH_TEST_ASSERT(!cb->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(cb->n_blocks()), 2);
    SMESH_TEST_ASSERT(cb->block(0)->name() == "white");
    SMESH_TEST_ASSERT(cb->block(1)->name() == "black");
    SMESH_TEST_EQ(cb->element_type(0), HEX8);
    SMESH_TEST_EQ(cb->element_type(1), HEX8);
    SMESH_TEST_EQ(cb->n_nodes(), static_cast<ptrdiff_t>(27));

    auto mixed = Mesh::create_hex8_tet4_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_ASSERT(!mixed->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(mixed->n_blocks()), 2);
    SMESH_TEST_ASSERT(mixed->block(0)->name() == "hex");
    SMESH_TEST_ASSERT(mixed->block(1)->name() == "tet");
    SMESH_TEST_EQ(mixed->element_type(0), HEX8);
    SMESH_TEST_EQ(mixed->element_type(1), TET4);
    SMESH_TEST_EQ(mixed->n_nodes(), static_cast<ptrdiff_t>(27));
    SMESH_TEST_EQ(mixed->n_elements(0), static_cast<ptrdiff_t>(4));
    SMESH_TEST_EQ(mixed->n_elements(1), static_cast<ptrdiff_t>(24));
    SMESH_TEST_EQ(check_unique_node_ids(*mixed), SMESH_TEST_SUCCESS);

    auto bidomain = Mesh::create_hex8_bidomain_cube(comm, 4, 2, 2);
    SMESH_TEST_ASSERT(bidomain != nullptr);
    SMESH_TEST_ASSERT(!bidomain->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(bidomain->n_blocks()), 2);
    SMESH_TEST_ASSERT(bidomain->block(0)->name() == "left");
    SMESH_TEST_ASSERT(bidomain->block(1)->name() == "right");
    SMESH_TEST_EQ(bidomain->n_elements(0), static_cast<ptrdiff_t>(8));
    SMESH_TEST_EQ(bidomain->n_elements(1), static_cast<ptrdiff_t>(8));

    auto sshex = Mesh::create_semistructured_hex_cube(comm, 2, 2, 2, 2);
    SMESH_TEST_ASSERT(sshex != nullptr);
    SMESH_TEST_ASSERT(!sshex->is_distributed());
    SMESH_TEST_EQ(sshex->element_type(0), proteus_hex_type(2));

    auto ring = Mesh::create_quad4_ring(comm, 0.5, 1.0, 2, 8);
    SMESH_TEST_ASSERT(ring != nullptr);
    SMESH_TEST_ASSERT(!ring->is_distributed());
    SMESH_TEST_EQ(ring->element_type(0), QUAD4);
    SMESH_TEST_EQ(ring->n_elements(), static_cast<ptrdiff_t>(16));

    return SMESH_TEST_SUCCESS;
}

static int test_create_hex8_matches_file() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr));
    }
    comm->broadcast(&token, 1, 0);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = std::max<ptrdiff_t>(comm->size(), 3);
    const ptrdiff_t nz = std::max<ptrdiff_t>(comm->size(), 2);

    char path_buffer[256];
    std::snprintf(path_buffer, sizeof(path_buffer), "/tmp/smesh_parallel_create_hex8_%d_%d",
                  comm->size(), token);
    const Path mesh_path(path_buffer);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(mesh_path.to_string());
        create_directory(mesh_path);
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto from_file = Mesh::create_from_file(comm, mesh_path);
    auto created   = Mesh::create_hex8_cube(comm, nx, ny, nz);
    SMESH_TEST_ASSERT(from_file != nullptr);
    SMESH_TEST_ASSERT(created != nullptr);
    SMESH_TEST_ASSERT(created->is_distributed());
    SMESH_TEST_ASSERT(from_file->is_distributed());

    auto cd = created->distributed();
    auto fd = from_file->distributed();
    SMESH_TEST_EQ(cd->n_nodes_global(), fd->n_nodes_global());
    SMESH_TEST_EQ(cd->n_nodes_owned(), fd->n_nodes_owned());
    SMESH_TEST_EQ(cd->n_nodes_shared(), fd->n_nodes_shared());
    SMESH_TEST_EQ(cd->n_nodes_ghosts(), fd->n_nodes_ghosts());
    SMESH_TEST_EQ(cd->n_nodes_aura(), fd->n_nodes_aura());
    SMESH_TEST_EQ(cd->n_elements_global(), fd->n_elements_global());
    SMESH_TEST_EQ(cd->n_elements_owned(), fd->n_elements_owned());
    SMESH_TEST_EQ(cd->n_elements_shared(), fd->n_elements_shared());
    SMESH_TEST_EQ(cd->n_elements_ghosts(), fd->n_elements_ghosts());
    SMESH_TEST_EQ(created->n_nodes(), from_file->n_nodes());
    SMESH_TEST_EQ(static_cast<int>(created->n_blocks()), 1);
    SMESH_TEST_ASSERT(created->block(0)->name() == "default");
    SMESH_TEST_EQ(created->element_type(0), HEX8);

    comm->barrier();
    if (comm->rank() == 0) {
        std::filesystem::remove_all(mesh_path.to_string());
    }
    comm->barrier();
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_create_tet4_quad4_vs_serial() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    const ptrdiff_t nx = std::max<ptrdiff_t>(comm->size(), 2);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;

    auto tet_serial = Mesh::create_tet4_cube(Communicator::self(), nx, ny, nz);
    auto tet        = Mesh::create_tet4_cube(comm, nx, ny, nz);
    SMESH_TEST_ASSERT(tet_serial != nullptr);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_ASSERT(tet->is_distributed());
    SMESH_TEST_EQ(tet->distributed()->n_nodes_global(), tet_serial->n_nodes());
    SMESH_TEST_EQ(tet->distributed()->n_elements_global(), tet_serial->n_elements());
    SMESH_TEST_EQ(tet->element_type(0), TET4);

    const ptrdiff_t qnx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t qny = std::max<ptrdiff_t>(comm->size(), 3);
    auto quad_serial = Mesh::create_quad4_square(Communicator::self(), qnx, qny);
    auto quad        = Mesh::create_quad4_square(comm, qnx, qny);
    SMESH_TEST_ASSERT(quad_serial != nullptr);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_ASSERT(quad->is_distributed());
    SMESH_TEST_EQ(quad->distributed()->n_nodes_global(), quad_serial->n_nodes());
    SMESH_TEST_EQ(quad->distributed()->n_elements_global(), quad_serial->n_elements());
    SMESH_TEST_EQ(quad->spatial_dimension(), 2);
    SMESH_TEST_EQ(quad->element_type(0), QUAD4);

    auto tri_serial = Mesh::create_tri3_square(Communicator::self(), qnx, qny);
    auto tri        = Mesh::create_tri3_square(comm, qnx, qny);
    SMESH_TEST_ASSERT(tri_serial != nullptr);
    SMESH_TEST_ASSERT(tri != nullptr);
    SMESH_TEST_ASSERT(tri->is_distributed());
    SMESH_TEST_EQ(tri->distributed()->n_nodes_global(), tri_serial->n_nodes());
    SMESH_TEST_EQ(tri->distributed()->n_elements_global(), tri_serial->n_elements());
    SMESH_TEST_EQ(tri->spatial_dimension(), 2);
    SMESH_TEST_EQ(tri->element_type(0), TRI3);

    return SMESH_TEST_SUCCESS;
#endif
}

static int test_create_checkerboard_and_hex_tet() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    const ptrdiff_t n = even_at_least(std::max<ptrdiff_t>(comm->size(), 2));

    auto cb_serial = Mesh::create_hex8_checkerboard_cube(Communicator::self(), n, n, n);
    auto cb        = Mesh::create_hex8_checkerboard_cube(comm, n, n, n);
    SMESH_TEST_ASSERT(cb_serial != nullptr);
    SMESH_TEST_ASSERT(cb != nullptr);
    SMESH_TEST_ASSERT(cb->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(cb->n_blocks()), 2);
    SMESH_TEST_ASSERT(cb->block(0)->name() == "white");
    SMESH_TEST_ASSERT(cb->block(1)->name() == "black");
    SMESH_TEST_EQ(cb->element_type(0), HEX8);
    SMESH_TEST_EQ(cb->element_type(1), HEX8);
    SMESH_TEST_EQ(cb->distributed()->n_nodes_global(), cb_serial->n_nodes());
    SMESH_TEST_EQ(cb->distributed()->n_elements_global(), cb_serial->n_elements());

    auto mixed_serial = Mesh::create_hex8_tet4_cube(Communicator::self(), n, n, n);
    auto mixed        = Mesh::create_hex8_tet4_cube(comm, n, n, n);
    SMESH_TEST_ASSERT(mixed_serial != nullptr);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_ASSERT(mixed->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(mixed->n_blocks()), 2);
    SMESH_TEST_ASSERT(mixed->block(0)->name() == "hex");
    SMESH_TEST_ASSERT(mixed->block(1)->name() == "tet");
    SMESH_TEST_EQ(mixed->element_type(0), HEX8);
    SMESH_TEST_EQ(mixed->element_type(1), TET4);
    SMESH_TEST_EQ(mixed->distributed()->n_nodes_global(), mixed_serial->n_nodes());
    SMESH_TEST_EQ(mixed->distributed()->n_elements_global(), mixed_serial->n_elements());
    SMESH_TEST_EQ(check_unique_node_ids(*mixed), SMESH_TEST_SUCCESS);

    return SMESH_TEST_SUCCESS;
#endif
}

static int test_create_remaining_a17_generators() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    const ptrdiff_t nx = even_at_least(std::max<ptrdiff_t>(comm->size(), 4));
    const ptrdiff_t ny = 4;
    const ptrdiff_t nz = 2;

    auto bidomain_serial = Mesh::create_hex8_bidomain_cube(Communicator::self(), nx, ny, nz);
    auto bidomain        = Mesh::create_hex8_bidomain_cube(comm, nx, ny, nz);
    SMESH_TEST_ASSERT(bidomain_serial != nullptr);
    SMESH_TEST_ASSERT(bidomain != nullptr);
    SMESH_TEST_ASSERT(bidomain->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(bidomain->n_blocks()), 2);
    SMESH_TEST_ASSERT(bidomain->block(0)->name() == "left");
    SMESH_TEST_ASSERT(bidomain->block(1)->name() == "right");
    SMESH_TEST_EQ(bidomain->element_type(0), HEX8);
    SMESH_TEST_EQ(bidomain->element_type(1), HEX8);
    SMESH_TEST_EQ(bidomain->distributed()->n_nodes_global(), bidomain_serial->n_nodes());
    SMESH_TEST_EQ(bidomain->distributed()->n_elements_global(), bidomain_serial->n_elements());

    auto ss_serial = Mesh::create_semistructured_hex_cube(Communicator::self(), 2, 2, 2, 2);
    auto ss        = Mesh::create_semistructured_hex_cube(comm, 2, 2, 2, 2);
    SMESH_TEST_ASSERT(ss_serial != nullptr);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), ss_serial->n_nodes());
    SMESH_TEST_EQ(ss->distributed()->n_elements_global(), ss_serial->n_elements());
    SMESH_TEST_EQ(ss->element_type(0), proteus_hex_type(2));

    auto ring_serial = Mesh::create_quad4_ring(Communicator::self(), 0.5, 1.0, 3, std::max<ptrdiff_t>(8, 2 * comm->size()));
    auto ring        = Mesh::create_quad4_ring(comm, 0.5, 1.0, 3, std::max<ptrdiff_t>(8, 2 * comm->size()));
    SMESH_TEST_ASSERT(ring_serial != nullptr);
    SMESH_TEST_ASSERT(ring != nullptr);
    SMESH_TEST_ASSERT(ring->is_distributed());
    SMESH_TEST_EQ(ring->distributed()->n_nodes_global(), ring_serial->n_nodes());
    SMESH_TEST_EQ(ring->distributed()->n_elements_global(), ring_serial->n_elements());
    SMESH_TEST_EQ(ring->spatial_dimension(), 3);
    SMESH_TEST_EQ(ring->element_type(0), QUAD4);

    auto hump_hex_serial = Mesh::create_wall_mounted_hump(Communicator::self(), HEX8, nx, ny, nz);
    auto hump_hex        = Mesh::create_wall_mounted_hump(comm, HEX8, nx, ny, nz);
    SMESH_TEST_ASSERT(hump_hex_serial != nullptr);
    SMESH_TEST_ASSERT(hump_hex != nullptr);
    SMESH_TEST_ASSERT(hump_hex->is_distributed());
    SMESH_TEST_ASSERT(hump_hex->block(0)->name() == "fluid");
    SMESH_TEST_EQ(hump_hex->distributed()->n_nodes_global(), hump_hex_serial->n_nodes());
    SMESH_TEST_EQ(hump_hex->distributed()->n_elements_global(), hump_hex_serial->n_elements());

    auto hump_tet_serial = Mesh::create_wall_mounted_hump(Communicator::self(), TET4, nx, ny, nz);
    auto hump_tet        = Mesh::create_wall_mounted_hump(comm, TET4, nx, ny, nz);
    SMESH_TEST_ASSERT(hump_tet_serial != nullptr);
    SMESH_TEST_ASSERT(hump_tet != nullptr);
    SMESH_TEST_ASSERT(hump_tet->is_distributed());
    SMESH_TEST_ASSERT(hump_tet->block(0)->name() == "fluid");
    SMESH_TEST_EQ(hump_tet->distributed()->n_nodes_global(), hump_tet_serial->n_nodes());
    SMESH_TEST_EQ(hump_tet->distributed()->n_elements_global(), hump_tet_serial->n_elements());

    auto sphere_hex_serial = Mesh::create_half_sphere(Communicator::self(), HEX8, 1.0, nx, ny, nz);
    auto sphere_hex        = Mesh::create_half_sphere(comm, HEX8, 1.0, nx, ny, nz);
    SMESH_TEST_ASSERT(sphere_hex_serial != nullptr);
    SMESH_TEST_ASSERT(sphere_hex != nullptr);
    SMESH_TEST_ASSERT(sphere_hex->is_distributed());
    SMESH_TEST_EQ(sphere_hex->distributed()->n_nodes_global(), sphere_hex_serial->n_nodes());
    SMESH_TEST_EQ(sphere_hex->distributed()->n_elements_global(), sphere_hex_serial->n_elements());

    auto sphere_tet_serial = Mesh::create_half_sphere(Communicator::self(), TET4, 1.0, nx, ny, nz);
    auto sphere_tet        = Mesh::create_half_sphere(comm, TET4, 1.0, nx, ny, nz);
    SMESH_TEST_ASSERT(sphere_tet_serial != nullptr);
    SMESH_TEST_ASSERT(sphere_tet != nullptr);
    SMESH_TEST_ASSERT(sphere_tet->is_distributed());
    SMESH_TEST_EQ(sphere_tet->distributed()->n_nodes_global(), sphere_tet_serial->n_nodes());
    SMESH_TEST_EQ(sphere_tet->distributed()->n_elements_global(), sphere_tet_serial->n_elements());

    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_create_size_one_serial);
    SMESH_RUN_TEST(test_create_hex8_matches_file);
    SMESH_RUN_TEST(test_create_tet4_quad4_vs_serial);
    SMESH_RUN_TEST(test_create_checkerboard_and_hex_tet);
    SMESH_RUN_TEST(test_create_remaining_a17_generators);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
