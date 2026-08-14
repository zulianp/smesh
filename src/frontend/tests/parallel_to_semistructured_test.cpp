#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <vector>

#include "smesh_conversion.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static std::shared_ptr<Mesh> split_first_half(const std::shared_ptr<Mesh> &mesh) {
    auto            out     = mesh->clone();
    const ptrdiff_t n       = out->n_elements(0);
    const ptrdiff_t n_split = n / 2;
    auto            parents = create_host_buffer<element_idx_t>(static_cast<size_t>(n_split));
    for (ptrdiff_t i = 0; i < n_split; ++i) {
        parents->data()[i] = static_cast<element_idx_t>(i);
    }
    if (out->split_block(parents, "part0", 0) != SMESH_SUCCESS) {
        return nullptr;
    }
    return out;
}

static std::shared_ptr<Mesh> create_hex8_tet4_serial(const ptrdiff_t nx, const ptrdiff_t ny, const ptrdiff_t nz) {
    auto            cube       = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
    const ptrdiff_t n_hex_all  = cube->n_elements();
    const ptrdiff_t n_hex_keep = n_hex_all / 2;
    const ptrdiff_t n_hex_conv = n_hex_all - n_hex_keep;
    auto            hex_src    = cube->elements(0)->data();
    auto            hex_keep   = create_host_buffer<idx_t>(8, static_cast<size_t>(n_hex_keep));
    for (int d = 0; d < 8; ++d) {
        std::memcpy(hex_keep->data()[d], hex_src[d], static_cast<size_t>(n_hex_keep) * sizeof(idx_t));
    }
    idx_t *hex_tail[8];
    for (int d = 0; d < 8; ++d) {
        hex_tail[d] = hex_src[d] + n_hex_keep;
    }
    auto tet_buf = create_host_buffer<idx_t>(4, static_cast<size_t>(n_hex_conv * 6));
    mesh_hex8_to_6x_tet4<idx_t>(n_hex_conv, hex_tail, tet_buf->data());
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(std::make_shared<Mesh::Block>("hex", HEX8, hex_keep));
    blocks.push_back(std::make_shared<Mesh::Block>("tet", TET4, tet_buf));
    return std::make_shared<Mesh>(Communicator::self(), blocks, cube->points());
}

static Path make_tmp_path(const char *prefix, const int token) {
    auto          comm = Communicator::world();
    char          buf[256];
    std::snprintf(buf, sizeof(buf), "/tmp/%s_%d_%d", prefix, comm->size(), token);
    return Path(buf);
}

static int check_owned_gids_unique(const Mesh &ss) {
#ifndef SMESH_ENABLE_MPI
    SMESH_UNUSED(ss);
    return SMESH_TEST_SUCCESS;
#else
    SMESH_TEST_ASSERT(ss.is_distributed());
    auto              dist    = ss.distributed();
    const ptrdiff_t   n_owned = dist->n_nodes_owned();
    const large_idx_t *map    = dist->node_mapping()->data();
    std::vector<large_idx_t> local(static_cast<size_t>(n_owned));
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
        local[static_cast<size_t>(i)] = map[i];
    }
    std::sort(local.begin(), local.end());
    SMESH_TEST_ASSERT(std::unique(local.begin(), local.end()) == local.end());

    ptrdiff_t n_owned_sum = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&n_owned, &n_owned_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, ss.comm()->get()));
    SMESH_TEST_EQ(n_owned_sum, dist->n_nodes_global());
    return SMESH_TEST_SUCCESS;
#endif
}

static int check_block_layout_copied(const Mesh &coarse, const Mesh &ss) {
    SMESH_TEST_EQ(static_cast<int>(ss.n_blocks()), static_cast<int>(coarse.n_blocks()));
    for (size_t b = 0; b < coarse.n_blocks(); ++b) {
        auto cb = coarse.block(b);
        auto sb = ss.block(b);
        SMESH_TEST_ASSERT(cb->name() == sb->name());
        SMESH_TEST_EQ(sb->n_elements(), cb->n_elements());
        SMESH_TEST_EQ(sb->n_elements_owned(), cb->n_elements_owned());
        SMESH_TEST_EQ(sb->n_elements_shared(), cb->n_elements_shared());
        SMESH_TEST_EQ(sb->n_elements_ghosts(), cb->n_elements_ghosts());
        SMESH_TEST_EQ(sb->n_elements(), sb->n_elements_owned() + sb->n_elements_ghosts());
    }
    SMESH_TEST_EQ(ss.distributed()->n_elements_global(), coarse.distributed()->n_elements_global());
    return SMESH_TEST_SUCCESS;
}

static int test_mpi_hex8_cube_to_ss() {
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
    const Path path = make_tmp_path("smesh_mpi_ss_hex8", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    const int       L  = 2;
    const ptrdiff_t serial_nnodes    = (static_cast<ptrdiff_t>(L) * nx + 1) *
                                    (static_cast<ptrdiff_t>(L) * ny + 1) *
                                    (static_cast<ptrdiff_t>(L) * nz + 1);
    const ptrdiff_t serial_nelements = nx * ny * nz;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(serial->n_elements(), serial_nelements);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());

    SMESH_TEST_ASSERT(to_semistructured(2, mesh, true, false) == nullptr);

    auto ss = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 1);
    SMESH_TEST_EQ(ss->element_type(0), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), sshex8_nxe(2));
    SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(ss->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(ss->n_nodes(), ss->distributed()->n_nodes_local());
    SMESH_TEST_EQ(ss->n_nodes(),
                  ss->distributed()->n_nodes_owned() + ss->distributed()->n_nodes_ghosts() +
                          ss->distributed()->n_nodes_aura());
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_checkerboard_to_ss() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 17;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_ss_checkerboard", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    const int       L  = 2;
    const ptrdiff_t serial_nnodes = (static_cast<ptrdiff_t>(L) * nx + 1) *
                                    (static_cast<ptrdiff_t>(L) * ny + 1) *
                                    (static_cast<ptrdiff_t>(L) * nz + 1);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_checkerboard_cube(Communicator::self(), nx, ny, nz);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(static_cast<int>(serial->n_blocks()), 2);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 2);
    SMESH_TEST_ASSERT(ss->block(0)->name() == "white");
    SMESH_TEST_ASSERT(ss->block(1)->name() == "black");
    SMESH_TEST_EQ(ss->element_type(0), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(ss->element_type(1), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_tet4_to_ss() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 31;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_ss_tet4", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    ptrdiff_t       serial_nnodes = 0;
    ptrdiff_t       n_blocks      = 0;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto single = Mesh::create_tet4_cube(Communicator::self(), nx, 2, 2);
        auto multi  = split_first_half(single);
        SMESH_TEST_ASSERT(multi != nullptr);
        n_blocks          = static_cast<ptrdiff_t>(multi->n_blocks());
        auto ss_serial    = to_semistructured(2, multi);
        SMESH_TEST_ASSERT(ss_serial != nullptr);
        serial_nnodes = ss_serial->n_nodes();
        SMESH_TEST_ASSERT(multi->write(path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes, 1, 0);
    comm->broadcast(&n_blocks, 1, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), static_cast<int>(n_blocks));
    SMESH_TEST_EQ(ss->element_type(0), semistructured_type(TET4, 2));
    SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), sstet4_nxe(2));
    SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static std::shared_ptr<Mesh> repeat_block_elements(const std::shared_ptr<Mesh> &mesh, const ptrdiff_t copies) {
    if (!mesh || copies < 1) {
        return nullptr;
    }
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        auto            src = mesh->block(b);
        const int       nxe = src->n_nodes_per_element();
        const ptrdiff_t n0  = src->n_elements();
        const ptrdiff_t n1  = n0 * copies;
        auto            dst = create_host_buffer<idx_t>(static_cast<size_t>(nxe), static_cast<size_t>(n1));
        for (int d = 0; d < nxe; ++d) {
            for (ptrdiff_t c = 0; c < copies; ++c) {
                std::memcpy(dst->data()[d] + c * n0, src->elements()->data()[d], static_cast<size_t>(n0) * sizeof(idx_t));
            }
        }
        blocks.push_back(std::make_shared<Mesh::Block>(src->name(), src->element_type(), dst));
    }
    return std::make_shared<Mesh>(mesh->comm(), blocks, mesh->points());
}

static int test_mpi_mixed_rejected() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 47;
    }
    comm->broadcast(&token, 1, 0);
    const Path mixed_path = make_tmp_path("smesh_mpi_ss_mixed", token);
    const Path hexdom_path = make_tmp_path("smesh_mpi_ss_hexdom", token + 1);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(mixed_path.to_string());
        std::filesystem::remove_all(hexdom_path.to_string());
        auto mixed = create_hex8_tet4_serial(4, 2, 2);
        SMESH_TEST_ASSERT(mixed != nullptr);
        SMESH_TEST_ASSERT(mixed->write(mixed_path) == SMESH_SUCCESS);
        auto hexdom = Mesh::create_hex_dominant_serial(Communicator::self());
        SMESH_TEST_ASSERT(hexdom != nullptr);
        hexdom = repeat_block_elements(hexdom, std::max<ptrdiff_t>(2 * comm->size(), 4));
        SMESH_TEST_ASSERT(hexdom != nullptr);
        SMESH_TEST_ASSERT(hexdom->write(hexdom_path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mixed = Mesh::create_from_file(comm, mixed_path);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_ASSERT(to_semistructured(2, mixed) == nullptr);

    auto hexdom = Mesh::create_from_file(comm, hexdom_path);
    SMESH_TEST_ASSERT(hexdom != nullptr);
    SMESH_TEST_ASSERT(to_semistructured(2, hexdom) == nullptr);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(mixed_path.to_string());
        std::filesystem::remove_all(hexdom_path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_mpi_hex8_cube_to_ss);
    SMESH_RUN_TEST(test_mpi_checkerboard_to_ss);
    SMESH_RUN_TEST(test_mpi_tet4_to_ss);
    SMESH_RUN_TEST(test_mpi_mixed_rejected);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

