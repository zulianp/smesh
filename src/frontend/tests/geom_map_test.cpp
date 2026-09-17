#include "smesh_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_test.hpp"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using namespace smesh;

static Path unique_tmp(const char *tag) {
    const auto token = static_cast<long long>(
            std::chrono::steady_clock::now().time_since_epoch().count());
    char buf[256];
    std::snprintf(buf, sizeof(buf), "/tmp/smesh_geom_map_%s_%lld", tag, token);
    return Path(buf);
}

static int test_factory_geom_maps() {
    auto comm = Communicator::self();
    auto hex  = Mesh::create_hex8_cube(comm, 2, 2, 2);
    auto quad = Mesh::create_quad4_square(comm, 2, 2);
    auto tet  = Mesh::create_tet4_cube(comm, 1, 1, 1);
    auto tri  = Mesh::create_tri3_square(comm, 2, 2);
    auto ss   = Mesh::create_semistructured_hex_cube(comm, 2, 1, 1, 1);
    auto lsh  = Mesh::create_hex8_lshape(comm, 4, 2, 2, 4, 2, 1, 1, 1);
    const std::vector<geom_t>    nxb = {0, 1, 2};
    const std::vector<geom_t>    nrb = {2, 2, 1};
    const std::vector<ptrdiff_t> nna = {2, 2};
    auto noz = Mesh::create_hex8_nozzle(comm, nxb, nrb, nna, -1, 3, 2, 1, 1);
    auto mixed = Mesh::create_hex8_tet4_cube(comm, 2, 2, 2);
    auto board = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto ring = Mesh::create_quad4_ring(comm, 1, 2, 2, 8);
    auto sph_hex = Mesh::create_hex8_half_sphere(comm, 1, 2, 2, 2);
    auto sph_tet = Mesh::create_tet4_half_sphere(comm, 1, 2, 2, 2);
    auto hump_hex = Mesh::create_wall_mounted_hump(comm, HEX8, 8, 4, 2);
    auto hump_tet = Mesh::create_wall_mounted_hump(comm, TET4, 4, 4, 2);
    auto hex27 = Mesh::create_cube(comm, HEX27, 1, 1, 1);
    auto tri6 = Mesh::create_square(comm, TRI6, 2, 2);
    auto refc = Mesh::create_hex8_reference_cube();
    auto hexdom = Mesh::create_hex_dominant_serial(comm);
    auto cyl = Mesh::create_hex_dominant_cylinder(comm, 1, 1, 1, 8, 2, 0);
    auto bidom = Mesh::create_hex8_bidomain_cube(comm, 2, 2, 2);

    SMESH_TEST_ASSERT(hex != nullptr);
    SMESH_TEST_EQ(hex->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_EQ(quad->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(tri != nullptr);
    SMESH_TEST_EQ(tri->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_EQ(ss->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(lsh != nullptr);
    SMESH_TEST_EQ(lsh->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(noz != nullptr);
    SMESH_TEST_EQ(noz->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_EQ(mixed->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(mixed->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(mixed->geom_map(1), AFFINE);
    SMESH_TEST_ASSERT(board != nullptr);
    SMESH_TEST_EQ(board->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(board->geom_map(1), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(ring != nullptr);
    SMESH_TEST_EQ(ring->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_ASSERT(sph_hex != nullptr);
    SMESH_TEST_EQ(sph_hex->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_ASSERT(sph_tet != nullptr);
    SMESH_TEST_EQ(sph_tet->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(hump_hex != nullptr);
    SMESH_TEST_EQ(hump_hex->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_ASSERT(hump_tet != nullptr);
    SMESH_TEST_EQ(hump_tet->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(hex27 != nullptr);
    SMESH_TEST_EQ(hex27->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(tri6 != nullptr);
    SMESH_TEST_EQ(tri6->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(refc != nullptr);
    SMESH_TEST_EQ(refc->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(hexdom != nullptr);
    SMESH_TEST_EQ(hexdom->n_blocks(), static_cast<size_t>(4));
    SMESH_TEST_EQ(hexdom->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(hexdom->geom_map(1), ISOPARAMETRIC);
    SMESH_TEST_EQ(hexdom->geom_map(2), AFFINE);
    SMESH_TEST_EQ(hexdom->geom_map(3), AFFINE);
    SMESH_TEST_ASSERT(cyl != nullptr);
    SMESH_TEST_EQ(cyl->geom_map(0), ISOPARAMETRIC);
    if (cyl->n_blocks() > 1) {
        SMESH_TEST_EQ(cyl->geom_map(1), AFFINE);
    }
    SMESH_TEST_ASSERT(bidom != nullptr);
    SMESH_TEST_EQ(bidom->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(bidom->geom_map(1), AXIS_ALIGNED);
    return SMESH_TEST_SUCCESS;
}

static int test_axis_aligned_rejected_on_tet() {
    auto tet = Mesh::create_tet4_cube(Communicator::self(), 1, 1, 1);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->geom_map(0), AFFINE);
    SMESH_TEST_EQ(tet->set_geom_map(0, AXIS_ALIGNED), SMESH_FAILURE);
    SMESH_TEST_EQ(tet->geom_map(0), AFFINE);
    SMESH_TEST_EQ(tet->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(tet->geom_map(0), ISOPARAMETRIC);
    return SMESH_TEST_SUCCESS;
}

static int test_clone_split_convert_promote() {
    auto hex = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);

    auto cloned = hex->clone();
    SMESH_TEST_ASSERT(cloned != nullptr);
    SMESH_TEST_EQ(cloned->geom_map(0), AXIS_ALIGNED);

    auto tets = convert_to(TET4, hex);
    SMESH_TEST_ASSERT(tets != nullptr);
    SMESH_TEST_EQ(tets->geom_map(0), AFFINE);

    auto ss = to_semistructured(2, hex, false, false);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_EQ(ss->geom_map(0), AXIS_ALIGNED);

    auto quad = Mesh::create_quad4_square(Communicator::self(), 2, 2);
    auto q9   = promote_to(QUAD9, quad);
    SMESH_TEST_ASSERT(q9 != nullptr);
    SMESH_TEST_EQ(q9->geom_map(0), AXIS_ALIGNED);

    auto parents = create_host_buffer<element_idx_t>(1);
    parents->data()[0] = 0;
    SMESH_TEST_EQ(cloned->split_block(parents, "part0", 0), SMESH_SUCCESS);
    SMESH_TEST_EQ(cloned->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(cloned->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(cloned->geom_map(1), AXIS_ALIGNED);
    return SMESH_TEST_SUCCESS;
}

static int test_write_read_roundtrip() {
    auto hex = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    const Path path = unique_tmp("hex");
    std::filesystem::remove_all(path.to_string());
    SMESH_TEST_EQ(hex->write(path), SMESH_SUCCESS);

    auto loaded = Mesh::create_from_file(Communicator::self(), path);
    SMESH_TEST_ASSERT(loaded != nullptr);
    SMESH_TEST_EQ(loaded->geom_map(0), AXIS_ALIGNED);
    std::filesystem::remove_all(path.to_string());

    auto mixed = Mesh::create_hex8_tet4_cube(Communicator::self(), 2, 2, 2);
    const Path mpath = unique_tmp("mixed");
    std::filesystem::remove_all(mpath.to_string());
    SMESH_TEST_EQ(mixed->write(mpath), SMESH_SUCCESS);
    auto mloaded = Mesh::create_from_file(Communicator::self(), mpath);
    SMESH_TEST_ASSERT(mloaded != nullptr);
    SMESH_TEST_EQ(mloaded->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(mloaded->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(mloaded->geom_map(1), AFFINE);
    std::filesystem::remove_all(mpath.to_string());
    return SMESH_TEST_SUCCESS;
}

static int test_detect_geom_maps() {
    auto comm = Communicator::self();
    auto hex  = Mesh::create_hex8_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    SMESH_TEST_EQ(hex->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_EQ(hex->detect_geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(hex->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_EQ(hex->detect_and_set_geom_map(0), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex->geom_map(0), AXIS_ALIGNED);

    auto tet = Mesh::create_tet4_cube(comm, 1, 1, 1);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(tet->detect_geom_map(0), AFFINE);

    auto quad = Mesh::create_quad4_square(comm, 2, 2);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_EQ(quad->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(quad->detect_geom_map(0), AXIS_ALIGNED);

    auto ss = Mesh::create_semistructured_hex_cube(comm, 2, 1, 1, 1);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_EQ(ss->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(ss->detect_geom_map(0), AXIS_ALIGNED);

    auto hex27 = Mesh::create_cube(comm, HEX27, 1, 1, 1);
    SMESH_TEST_ASSERT(hex27 != nullptr);
    SMESH_TEST_EQ(hex27->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex27->detect_geom_map(0), AXIS_ALIGNED);

    auto q9 = promote_to(QUAD9, Mesh::create_quad4_square(comm, 2, 2));
    SMESH_TEST_ASSERT(q9 != nullptr);
    SMESH_TEST_EQ(q9->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(q9->detect_geom_map(0), AXIS_ALIGNED);

    auto t10 = promote_to(TET10, Mesh::create_tet4_cube(comm, 1, 1, 1));
    SMESH_TEST_ASSERT(t10 != nullptr);
    SMESH_TEST_EQ(t10->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(t10->detect_geom_map(0), AFFINE);

    const std::vector<geom_t>    nxb = {0, 1, 2};
    const std::vector<geom_t>    nrb = {2, 2, 1};
    const std::vector<ptrdiff_t> nna = {2, 2};
    auto noz = Mesh::create_hex8_nozzle(comm, nxb, nrb, nna, -1, 3, 2, 1, 1);
    SMESH_TEST_ASSERT(noz != nullptr);
    SMESH_TEST_EQ(noz->detect_geom_map(0), ISOPARAMETRIC);

    auto sph = Mesh::create_hex8_half_sphere(comm, 1, 2, 2, 2);
    SMESH_TEST_ASSERT(sph != nullptr);
    SMESH_TEST_EQ(sph->detect_geom_map(0), ISOPARAMETRIC);

    auto mixed = Mesh::create_hex8_tet4_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_EQ(mixed->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(mixed->set_geom_map(1, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(mixed->detect_and_set_geom_maps(), SMESH_SUCCESS);
    SMESH_TEST_EQ(mixed->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(mixed->geom_map(1), AFFINE);
    return SMESH_TEST_SUCCESS;
}

static int test_detect_sheared_and_warped_hex() {
    auto hex = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    auto            p  = hex->points()->data();
    const ptrdiff_t nn = hex->n_nodes();
    for (ptrdiff_t i = 0; i < nn; ++i) {
        p[1][i] += geom_t(0.25) * p[0][i];
    }
    SMESH_TEST_EQ(hex->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(hex->detect_geom_map(0), AFFINE);
    SMESH_TEST_EQ(hex->detect_and_set_geom_map(0), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex->geom_map(0), AFFINE);

    p[0][0] += geom_t(0.2);
    SMESH_TEST_EQ(hex->detect_geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_EQ(hex->detect_and_set_geom_map(0), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex->geom_map(0), ISOPARAMETRIC);
    return SMESH_TEST_SUCCESS;
}

#ifdef SMESH_ENABLE_MPI
static int test_mpi_factory_geom_maps() {
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    auto hex  = Mesh::create_hex8_cube(comm, 4, 2, 2);
    auto quad = Mesh::create_quad4_square(comm, 4, 2);
    auto tet  = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto lsh  = Mesh::create_hex8_lshape(comm, 4, 2, 2, 4, 2, 1, 1, 1);
    auto ring = Mesh::create_quad4_ring(comm, 1, 2, 2, 8);
    auto sph_hex = Mesh::create_hex8_half_sphere(comm, 1, 4, 2, 2);
    auto sph_tet = Mesh::create_tet4_half_sphere(comm, 1, 2, 2, 2);
    auto hump = Mesh::create_wall_mounted_hump(comm, HEX8, 8, 4, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    SMESH_TEST_EQ(hex->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(quad != nullptr);
    SMESH_TEST_EQ(quad->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(lsh != nullptr);
    SMESH_TEST_EQ(lsh->geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_ASSERT(ring != nullptr);
    SMESH_TEST_EQ(ring->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_ASSERT(sph_hex != nullptr);
    SMESH_TEST_EQ(sph_hex->geom_map(0), ISOPARAMETRIC);
    SMESH_TEST_ASSERT(sph_tet != nullptr);
    SMESH_TEST_EQ(sph_tet->geom_map(0), AFFINE);
    SMESH_TEST_ASSERT(hump != nullptr);
    SMESH_TEST_EQ(hump->geom_map(0), ISOPARAMETRIC);
    return SMESH_TEST_SUCCESS;
}

static int test_mpi_detect_geom_maps() {
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    auto hex = Mesh::create_hex8_cube(comm, 4, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    SMESH_TEST_EQ(hex->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex->detect_geom_map(0), AXIS_ALIGNED);
    SMESH_TEST_EQ(hex->detect_and_set_geom_map(0), SMESH_SUCCESS);
    SMESH_TEST_EQ(hex->geom_map(0), AXIS_ALIGNED);

    auto tet = Mesh::create_tet4_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->set_geom_map(0, ISOPARAMETRIC), SMESH_SUCCESS);
    SMESH_TEST_EQ(tet->detect_geom_map(0), AFFINE);

    const std::vector<geom_t>    nxb = {0, 1, 2};
    const std::vector<geom_t>    nrb = {2, 2, 1};
    const std::vector<ptrdiff_t> nna = {2, 2};
    auto noz = Mesh::create_hex8_nozzle(comm, nxb, nrb, nna, -1, 3, 2, 1, 1);
    SMESH_TEST_ASSERT(noz != nullptr);
    SMESH_TEST_EQ(noz->detect_geom_map(0), ISOPARAMETRIC);
    return SMESH_TEST_SUCCESS;
}
#endif

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_factory_geom_maps);
    SMESH_RUN_TEST(test_axis_aligned_rejected_on_tet);
    SMESH_RUN_TEST(test_clone_split_convert_promote);
    SMESH_RUN_TEST(test_write_read_roundtrip);
    SMESH_RUN_TEST(test_detect_geom_maps);
    SMESH_RUN_TEST(test_detect_sheared_and_warped_hex);
#ifdef SMESH_ENABLE_MPI
    SMESH_RUN_TEST(test_mpi_factory_geom_maps);
    SMESH_RUN_TEST(test_mpi_detect_geom_maps);
#endif
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
