#include "smesh_adjacency.hpp"
#include "smesh_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_path.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
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

static void write_bin(const Path &path, const void *data, size_t nbytes) {
    FILE *fp = std::fopen(path.c_str(), "wb");
    if (fp) {
        std::fwrite(data, 1, nbytes, fp);
        std::fclose(fp);
    }
}

static int test_stale_coord_files_do_not_scramble() {
    auto mesh = Mesh::create_tet4_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const Path path = unique_tmp("stale_coords");
    std::filesystem::remove_all(path.to_string());
    SMESH_TEST_EQ(mesh->write(path), SMESH_SUCCESS);

    const ptrdiff_t nn = mesh->n_nodes();
    auto *orig = mesh->points()->data();
    std::vector<double> x64(static_cast<size_t>(nn));
    std::vector<float> x32(static_cast<size_t>(nn));
    for (ptrdiff_t i = 0; i < nn; ++i) {
        x64[static_cast<size_t>(i)] = static_cast<double>(orig[0][i]);
        x32[static_cast<size_t>(i)] = static_cast<float>(orig[0][i]);
    }
    if ((path / "x.float32").exists()) {
        write_bin(path / "x.float64", x64.data(), x64.size() * sizeof(double));
    } else {
        write_bin(path / "x.float32", x32.data(), x32.size() * sizeof(float));
    }

    auto loaded = Mesh::create_from_file(Communicator::self(), path);
    SMESH_TEST_ASSERT(loaded != nullptr);
    SMESH_TEST_EQ(loaded->n_nodes(), nn);
    auto *got = loaded->points()->data();
    for (ptrdiff_t i = 0; i < nn; ++i) {
        SMESH_TEST_ASSERT(std::abs(got[0][i] - orig[0][i]) < geom_t(1e-6));
        SMESH_TEST_ASSERT(std::abs(got[1][i] - orig[1][i]) < geom_t(1e-6));
        SMESH_TEST_ASSERT(std::abs(got[2][i] - orig[2][i]) < geom_t(1e-6));
    }

    std::vector<double> y64(static_cast<size_t>(nn));
    std::vector<double> z64(static_cast<size_t>(nn));
    for (ptrdiff_t i = 0; i < nn; ++i) {
        x64[static_cast<size_t>(i)] = static_cast<double>(orig[0][i]);
        y64[static_cast<size_t>(i)] = static_cast<double>(orig[1][i]);
        z64[static_cast<size_t>(i)] = static_cast<double>(orig[2][i]);
    }
    write_bin(path / "x.float64", x64.data(), x64.size() * sizeof(double));
    write_bin(path / "y.float64", y64.data(), y64.size() * sizeof(double));
    write_bin(path / "z.float64", z64.data(), z64.size() * sizeof(double));

    std::vector<float> junk(static_cast<size_t>(nn), 999.f);
    write_bin(path / "x.float32", junk.data(), junk.size() * sizeof(float));
    write_bin(path / "y.float32", junk.data(), junk.size() * sizeof(float));
    write_bin(path / "z.float32", junk.data(), junk.size() * sizeof(float));

    {
        const Path meta = path / "meta.yaml";
        std::ifstream in(meta.c_str());
        std::string yaml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
        for (size_t at = 0; (at = yaml.find(".float32", at)) != std::string::npos;) {
            yaml.replace(at, 8, ".float64");
            at += 8;
        }
        std::ofstream out(meta.c_str());
        out << yaml;
    }

    auto loaded64 = Mesh::create_from_file(Communicator::self(), path);
    SMESH_TEST_ASSERT(loaded64 != nullptr);
    auto *got64 = loaded64->points()->data();
    for (ptrdiff_t i = 0; i < nn; ++i) {
        SMESH_TEST_ASSERT(std::abs(got64[0][i] - orig[0][i]) < geom_t(1e-5));
        SMESH_TEST_ASSERT(std::abs(got64[1][i] - orig[1][i]) < geom_t(1e-5));
        SMESH_TEST_ASSERT(std::abs(got64[2][i] - orig[2][i]) < geom_t(1e-5));
    }

    std::filesystem::remove_all(path.to_string());
    return SMESH_TEST_SUCCESS;
}

static int test_named_single_block_write() {
    auto mesh = Mesh::create_tet4_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    mesh->block(0)->set_name("FLUID");
    const Path path = unique_tmp("named_block");
    std::filesystem::remove_all(path.to_string());
    SMESH_TEST_EQ(mesh->write(path), SMESH_SUCCESS);
    SMESH_TEST_ASSERT(std::filesystem::is_directory((path / "blocks" / "FLUID").to_string()));
    auto loaded = Mesh::create_from_file(Communicator::self(), path);
    SMESH_TEST_ASSERT(loaded != nullptr);
    SMESH_TEST_EQ(loaded->n_blocks(), static_cast<size_t>(1));
    SMESH_TEST_ASSERT(loaded->block(0)->name() == "FLUID");
    std::filesystem::remove_all(path.to_string());
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

static geom_t tet_volume(const geom_t *const *p, const idx_t a, const idx_t b, const idx_t c,
                         const idx_t d) {
    const geom_t abx = p[0][b] - p[0][a];
    const geom_t aby = p[1][b] - p[1][a];
    const geom_t abz = p[2][b] - p[2][a];
    const geom_t acx = p[0][c] - p[0][a];
    const geom_t acy = p[1][c] - p[1][a];
    const geom_t acz = p[2][c] - p[2][a];
    const geom_t adx = p[0][d] - p[0][a];
    const geom_t ady = p[1][d] - p[1][a];
    const geom_t adz = p[2][d] - p[2][a];
    const geom_t cx  = aby * acz - abz * acy;
    const geom_t cy  = abz * acx - abx * acz;
    const geom_t cz  = abx * acy - aby * acx;
    return (cx * adx + cy * ady + cz * adz) / static_cast<geom_t>(6);
}

static void sort6(geom_t *v) {
    for (int i = 1; i < 6; ++i) {
        const geom_t x = v[i];
        int          j = i;
        while (j > 0 && v[j - 1] > x) {
            v[j] = v[j - 1];
            --j;
        }
        v[j] = x;
    }
}

static int test_hex8_quality_tet4() {
    auto comm = Communicator::self();
    auto hex  = Mesh::create_hex8_cube(comm, 1, 1, 1);
    SMESH_TEST_ASSERT(hex != nullptr);
    auto tet = convert_to(TET4, hex);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->n_nodes(), static_cast<ptrdiff_t>(15));
    SMESH_TEST_EQ(tet->n_elements(), static_cast<ptrdiff_t>(24));
    SMESH_TEST_EQ(tet->geom_map(0), AFFINE);
    SMESH_TEST_EQ(tet->n_nodes(), hex->n_nodes() + 7);

    auto         p     = tet->points()->data();
    auto         el    = tet->elements(0)->data();
    geom_t       vol   = 0;
    geom_t       edges0[6];
    const geom_t tol   = static_cast<geom_t>(1e-9);
    for (ptrdiff_t e = 0; e < tet->n_elements(); ++e) {
        const idx_t n[4] = {el[0][e], el[1][e], el[2][e], el[3][e]};
        const geom_t ve  = tet_volume(p, n[0], n[1], n[2], n[3]);
        SMESH_TEST_ASSERT(ve > 0);
        vol += ve;
        geom_t edges[6];
        int    k = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                const geom_t dx = p[0][n[i]] - p[0][n[j]];
                const geom_t dy = p[1][n[i]] - p[1][n[j]];
                const geom_t dz = p[2][n[i]] - p[2][n[j]];
                edges[k++]      = dx * dx + dy * dy + dz * dz;
            }
        }
        sort6(edges);
        if (e == 0) {
            for (int i = 0; i < 6; ++i) {
                edges0[i] = edges[i];
            }
        } else {
            for (int i = 0; i < 6; ++i) {
                SMESH_TEST_ASSERT(std::fabs(edges[i] - edges0[i]) < tol);
            }
        }
    }
    const geom_t vol_tol = static_cast<geom_t>(1e-5);
    SMESH_TEST_ASSERT(std::fabs(vol - static_cast<geom_t>(1)) < vol_tol);

    auto pair = Mesh::create_hex8_cube(comm, 2, 1, 1, 0, 0, 0, 2, 1, 1);
    auto pair_tet = convert_to(TET4, pair);
    SMESH_TEST_ASSERT(pair_tet != nullptr);
    SMESH_TEST_EQ(pair_tet->n_nodes(), static_cast<ptrdiff_t>(25));
    SMESH_TEST_EQ(pair_tet->n_elements(), static_cast<ptrdiff_t>(48));
    auto   pp      = pair_tet->points()->data();
    auto   pe      = pair_tet->elements(0)->data();
    int    n_face  = 0;
    geom_t pair_vol = 0;
    idx_t  tris[8][3];
    for (ptrdiff_t i = 0; i < pair_tet->n_nodes(); ++i) {
        if (std::fabs(pp[0][i] - static_cast<geom_t>(1)) < tol &&
            std::fabs(pp[1][i] - static_cast<geom_t>(0.5)) < tol &&
            std::fabs(pp[2][i] - static_cast<geom_t>(0.5)) < tol) {
            ++n_face;
        }
    }
    SMESH_TEST_EQ(n_face, 1);
    int ntri = 0;
    for (ptrdiff_t e = 0; e < pair_tet->n_elements(); ++e) {
        const idx_t n[4] = {pe[0][e], pe[1][e], pe[2][e], pe[3][e]};
        const geom_t ve  = tet_volume(pp, n[0], n[1], n[2], n[3]);
        SMESH_TEST_ASSERT(ve > 0);
        pair_vol += ve;
        int on_face = 0;
        for (int k = 0; k < 3; ++k) {
            if (std::fabs(pp[0][n[k]] - static_cast<geom_t>(1)) < tol) {
                ++on_face;
            }
        }
        if (on_face == 3) {
            SMESH_TEST_ASSERT(ntri < 8);
            tris[ntri][0] = n[0];
            tris[ntri][1] = n[1];
            tris[ntri][2] = n[2];
            ++ntri;
        }
    }
    SMESH_TEST_EQ(ntri, 8);
    SMESH_TEST_ASSERT(std::fabs(pair_vol - static_cast<geom_t>(2)) < static_cast<geom_t>(1e-5));
    for (int i = 0; i < ntri; ++i) {
        int mates = 0;
        for (int j = 0; j < ntri; ++j) {
            if (i == j) {
                continue;
            }
            for (int r = 0; r < 3; ++r) {
                if (tris[j][0] == tris[i][r] && tris[j][1] == tris[i][(r + 2) % 3] &&
                    tris[j][2] == tris[i][(r + 1) % 3]) {
                    ++mates;
                }
            }
        }
        SMESH_TEST_EQ(mates, 1);
    }

    auto board = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto cube  = Mesh::create_hex8_cube(comm, 2, 2, 2);
    auto board_tet = convert_to(TET4, board);
    auto cube_tet  = convert_to(TET4, cube);
    SMESH_TEST_ASSERT(board_tet != nullptr);
    SMESH_TEST_ASSERT(cube_tet != nullptr);
    SMESH_TEST_EQ(board_tet->n_nodes(), cube_tet->n_nodes());
    SMESH_TEST_EQ(board_tet->n_elements(), cube_tet->n_elements());
    SMESH_TEST_EQ(cube_tet->n_nodes(), static_cast<ptrdiff_t>(71));
    SMESH_TEST_EQ(cube_tet->n_elements(), static_cast<ptrdiff_t>(192));
    return SMESH_TEST_SUCCESS;
}

static int sides_lie_on_plane(const Mesh &mesh, const Sideset &ss, const geom_t nx, const geom_t ny,
                              const geom_t nz, const geom_t distance) {
    auto block = mesh.block(ss.block_id());
    SMESH_TEST_ASSERT(block != nullptr);
    LocalSideTable lst;
    SMESH_TEST_EQ(lst.fill(block->element_type()), SMESH_SUCCESS);
    const int n_sides = elem_num_sides(block->element_type());
    auto el = block->elements()->data();
    auto p = mesh.points()->data();
    const int dim = mesh.spatial_dimension();
    const element_idx_t *parent = ss.parent()->data();
    const i16 *lfi = ss.lfi()->data();
    const geom_t tol = static_cast<geom_t>(1e-5);
    for (ptrdiff_t i = 0; i < ss.size(); ++i) {
        const int side = lfi[i];
        SMESH_TEST_ASSERT(side >= 0 && side < n_sides);
        SMESH_TEST_ASSERT(parent[i] >= 0 && static_cast<ptrdiff_t>(parent[i]) < block->n_elements());
        const int n_nodes = lst.nnxs_side[side];
        for (int k = 0; k < n_nodes; ++k) {
            const idx_t node = el[lst(side, k)][parent[i]];
            const geom_t x = dim > 0 ? p[0][node] : static_cast<geom_t>(0);
            const geom_t y = dim > 1 ? p[1][node] : static_cast<geom_t>(0);
            const geom_t z = dim > 2 ? p[2][node] : static_cast<geom_t>(0);
            SMESH_TEST_ASSERT(std::fabs(nx * x + ny * y + nz * z - distance) < tol);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_convert_keeps_sets() {
    auto comm = Communicator::self();
    auto hex = Mesh::create_hex8_cube(comm, 1, 1, 1);
    SMESH_TEST_ASSERT(hex != nullptr);
    auto face = Sideset::create_from_plane(hex, 0, 1, 0, 0);
    SMESH_TEST_EQ(face.size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(face[0]->size(), static_cast<ptrdiff_t>(1));
    hex->add_sideset("bottom", face[0]);
    auto ns_buf = create_nodeset_from_sideset(hex, face[0]);
    SMESH_TEST_ASSERT(ns_buf != nullptr);
    SMESH_TEST_EQ(ns_buf->size(), static_cast<size_t>(4));
    hex->add_nodeset("bottom_nodes", Nodeset::create(comm, ns_buf));

    auto tet = convert_to(TET4, hex);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->block(0)->name(), std::string("default"));
    SMESH_TEST_EQ(tet->sidesets("bottom").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(tet->sidesets("bottom")[0]->size(), face[0]->size() * 4);
    SMESH_TEST_EQ(tet->sidesets("bottom")[0]->block_id(), face[0]->block_id());
    SMESH_TEST_EQ(sides_lie_on_plane(*tet, *tet->sidesets("bottom")[0], 0, 1, 0, 0), SMESH_TEST_SUCCESS);

    auto kept = tet->nodesets("bottom_nodes");
    SMESH_TEST_EQ(kept.size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(static_cast<size_t>(kept[0]->size()), ns_buf->size());
    auto hp = hex->points()->data();
    auto tp = tet->points()->data();
    const idx_t *coarse_ids = ns_buf->data();
    const idx_t *fine_ids = kept[0]->nodes()->data();
    const geom_t tol = static_cast<geom_t>(1e-6);
    for (ptrdiff_t i = 0; i < kept[0]->size(); ++i) {
        SMESH_TEST_EQ(fine_ids[i], coarse_ids[i]);
        for (int d = 0; d < 3; ++d) {
            SMESH_TEST_ASSERT(std::fabs(tp[d][fine_ids[i]] - hp[d][coarse_ids[i]]) < tol);
        }
    }

    auto pair = Mesh::create_hex8_bidomain_cube(comm, 2, 1, 1);
    SMESH_TEST_ASSERT(pair != nullptr);
    auto right = Sideset::create_from_plane(pair, 1, 0, 0, 1, static_cast<geom_t>(1e-6), {"right"});
    SMESH_TEST_EQ(right.size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(right[0]->block_id(), static_cast<block_idx_t>(1));
    SMESH_TEST_EQ(right[0]->size(), static_cast<ptrdiff_t>(1));
    pair->add_sideset("xmax", right[0]);
    auto pair_tet = convert_to(TET4, pair);
    SMESH_TEST_ASSERT(pair_tet != nullptr);
    SMESH_TEST_EQ(pair_tet->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(pair_tet->block(0)->name(), std::string("left"));
    SMESH_TEST_EQ(pair_tet->block(1)->name(), std::string("right"));
    SMESH_TEST_EQ(pair_tet->sidesets("xmax").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(pair_tet->sidesets("xmax")[0]->block_id(), static_cast<block_idx_t>(1));
    SMESH_TEST_EQ(pair_tet->sidesets("xmax")[0]->size(), right[0]->size() * 4);
    SMESH_TEST_EQ(sides_lie_on_plane(*pair_tet, *pair_tet->sidesets("xmax")[0], 1, 0, 0, 1), SMESH_TEST_SUCCESS);

    auto quad = Mesh::create_quad4_square(comm, 1, 1);
    SMESH_TEST_ASSERT(quad != nullptr);
    auto edge = Sideset::create_from_plane(quad, 1, 0, 0, 0);
    SMESH_TEST_EQ(edge.size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(edge[0]->size(), static_cast<ptrdiff_t>(1));
    quad->add_sideset("xmin", edge[0]);
    auto tri = convert_to(TRI3, quad);
    SMESH_TEST_ASSERT(tri != nullptr);
    SMESH_TEST_EQ(tri->block(0)->name(), std::string("default"));
    SMESH_TEST_EQ(tri->element_type(0), TRI3);
    SMESH_TEST_EQ(tri->sidesets("xmin").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(tri->sidesets("xmin")[0]->size(), edge[0]->size());
    SMESH_TEST_EQ(sides_lie_on_plane(*tri, *tri->sidesets("xmin")[0], 1, 0, 0, 0), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_factory_geom_maps);
    SMESH_RUN_TEST(test_axis_aligned_rejected_on_tet);
    SMESH_RUN_TEST(test_clone_split_convert_promote);
    SMESH_RUN_TEST(test_hex8_quality_tet4);
    SMESH_RUN_TEST(test_convert_keeps_sets);
    SMESH_RUN_TEST(test_write_read_roundtrip);
    SMESH_RUN_TEST(test_stale_coord_files_do_not_scramble);
    SMESH_RUN_TEST(test_named_single_block_write);
    SMESH_RUN_TEST(test_detect_geom_maps);
    SMESH_RUN_TEST(test_detect_sheared_and_warped_hex);
#ifdef SMESH_ENABLE_MPI
    SMESH_RUN_TEST(test_mpi_factory_geom_maps);
    SMESH_RUN_TEST(test_mpi_detect_geom_maps);
#endif
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
