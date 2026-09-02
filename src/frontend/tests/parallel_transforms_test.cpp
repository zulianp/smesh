#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <utility>
#include <vector>

#include "smesh_distributed_base.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static Path make_tmp_path(const char *prefix, const int token) {
    auto comm = Communicator::world();
    char buf[256];
    std::snprintf(buf, sizeof(buf), "/tmp/%s_%d_%d", prefix, comm->size(), token);
    return Path(buf);
}

static std::shared_ptr<Mesh> create_quad4_square_3d(const ptrdiff_t nx, const ptrdiff_t ny) {
    auto q2 = Mesh::create_quad4_square(Communicator::self(), nx, ny);
    auto p3 = create_host_buffer<geom_t>(3, static_cast<size_t>(q2->n_nodes()));
    std::memcpy(p3->data()[0], q2->points()->data()[0], static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
    std::memcpy(p3->data()[1], q2->points()->data()[1], static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
    std::memset(p3->data()[2], 0, static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
    return std::make_shared<Mesh>(Communicator::self(), q2->blocks(), p3);
}

static int check_owned_gids_unique(const Mesh &mesh) {
#ifndef SMESH_ENABLE_MPI
    SMESH_UNUSED(mesh);
    return SMESH_TEST_SUCCESS;
#else
    SMESH_TEST_ASSERT(mesh.is_distributed());
    auto               dist    = mesh.distributed();
    const ptrdiff_t    n_owned = dist->n_nodes_owned();
    const large_idx_t *map     = dist->node_mapping()->data();
    std::vector<large_idx_t> local(static_cast<size_t>(n_owned));
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
        local[static_cast<size_t>(i)] = map[i];
    }
    std::sort(local.begin(), local.end());
    SMESH_TEST_ASSERT(std::unique(local.begin(), local.end()) == local.end());

    ptrdiff_t n_owned_sum = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&n_owned, &n_owned_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, mesh.comm()->get()));
    SMESH_TEST_EQ(n_owned_sum, dist->n_nodes_global());
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_hex8_refine() {
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
    const Path path = make_tmp_path("smesh_mpi_xf_hex_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
        SMESH_TEST_ASSERT(serial != nullptr);
        auto refined = refine(serial, 1);
        SMESH_TEST_ASSERT(refined != nullptr);
        serial_nnodes    = refined->n_nodes();
        serial_nelements = refined->n_elements();
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes, 1, 0);
    comm->broadcast(&serial_nelements, 1, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());

    auto refined = refine(mesh, 1);
    SMESH_TEST_ASSERT(refined != nullptr);
    SMESH_TEST_ASSERT(refined->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), static_cast<int>(mesh->n_blocks()));
    SMESH_TEST_EQ(refined->element_type(0), HEX8);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 8);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_plane(mesh, 1, 0, 0, 0.0);
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), coarse_ss[0]->size() * 4);
    SMESH_TEST_EQ(mapped->block_id(), coarse_ss[0]->block_id());
    for (ptrdiff_t i = 0; i < mapped->size(); ++i) {
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] >= 0);
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] < refined->n_elements(mapped->block_id()));
    }

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_quad_extrude() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 1;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_quad_extrude", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_quad4_square_3d(nx, ny);
        SMESH_TEST_ASSERT(serial != nullptr);
        auto extruded = extrude(serial, 1.0, 2);
        SMESH_TEST_ASSERT(extruded != nullptr);
        serial_nnodes    = extruded->n_nodes();
        serial_nelements = extruded->n_elements();
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes, 1, 0);
    comm->broadcast(&serial_nelements, 1, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());

    auto extruded = extrude(mesh, 1.0, 2);
    SMESH_TEST_ASSERT(extruded != nullptr);
    SMESH_TEST_ASSERT(extruded->is_distributed());
    SMESH_TEST_EQ(extruded->element_type(0), HEX8);
    SMESH_TEST_EQ(extruded->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(extruded->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(check_owned_gids_unique(*extruded), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_clone_convert() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 2;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_clone_convert", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, 2, 2);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());

    auto cloned = mesh->clone();
    SMESH_TEST_ASSERT(cloned != nullptr);
    SMESH_TEST_ASSERT(cloned->is_distributed());
    SMESH_TEST_EQ(cloned->distributed()->n_nodes_global(), mesh->distributed()->n_nodes_global());
    SMESH_TEST_EQ(cloned->distributed()->n_elements_global(), mesh->distributed()->n_elements_global());
    SMESH_TEST_EQ(cloned->n_nodes(), mesh->n_nodes());
    SMESH_TEST_EQ(cloned->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned());
    SMESH_TEST_EQ(check_owned_gids_unique(*cloned), SMESH_TEST_SUCCESS);

    auto tet = convert_to(TET4, mesh);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_ASSERT(tet->is_distributed());
    SMESH_TEST_EQ(tet->element_type(0), TET4);
    SMESH_TEST_EQ(tet->distributed()->n_nodes_global(), mesh->distributed()->n_nodes_global());
    SMESH_TEST_EQ(tet->distributed()->n_elements_global(), mesh->distributed()->n_elements_global() * 6);
    SMESH_TEST_EQ(tet->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 6);
    SMESH_TEST_EQ(check_owned_gids_unique(*tet), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_ss_derefine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 3;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_derefine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, 2, 2);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const ptrdiff_t n_coarse = mesh->distributed()->n_nodes_global();

    auto ss = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());

    auto coarse = derefine(ss, 1);
    SMESH_TEST_ASSERT(coarse != nullptr);
    SMESH_TEST_ASSERT(coarse->is_distributed());
    SMESH_TEST_EQ(coarse->distributed()->n_nodes_global(), n_coarse);
    SMESH_TEST_EQ(coarse->distributed()->n_elements_global(), mesh->distributed()->n_elements_global());
    SMESH_TEST_EQ(check_owned_gids_unique(*coarse), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_mesh_sideset_io() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 11;
    }
    comm->broadcast(&token, 1, 0);
    const Path path  = make_tmp_path("smesh_mpi_mesh_ss", token);
    const Path path2 = make_tmp_path("smesh_mpi_mesh_ss_w", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        std::filesystem::remove_all(path2.to_string());
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, 2, 2);
        SMESH_TEST_ASSERT(serial != nullptr);
        auto left = Sideset::create_from_selector(
                serial, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
        SMESH_TEST_EQ(left.size(), static_cast<size_t>(1));
        serial->add_sideset("left", left[0]);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto loaded = mesh->sidesets("left");
    SMESH_TEST_EQ(loaded.size(), static_cast<size_t>(1));

    auto recreated = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_EQ(recreated.size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(loaded[0]->block_id(), recreated[0]->block_id());
    SMESH_TEST_EQ(loaded[0]->size(), recreated[0]->size());

    std::vector<std::pair<element_idx_t, i16>> keys_loaded(static_cast<size_t>(loaded[0]->size()));
    std::vector<std::pair<element_idx_t, i16>> keys_recreated(
            static_cast<size_t>(recreated[0]->size()));
    for (ptrdiff_t i = 0; i < loaded[0]->size(); ++i) {
        keys_loaded[static_cast<size_t>(i)] =
                std::make_pair(loaded[0]->parent()->data()[i], loaded[0]->lfi()->data()[i]);
        keys_recreated[static_cast<size_t>(i)] =
                std::make_pair(recreated[0]->parent()->data()[i], recreated[0]->lfi()->data()[i]);
    }
    std::sort(keys_loaded.begin(), keys_loaded.end());
    std::sort(keys_recreated.begin(), keys_recreated.end());
    for (size_t i = 0; i < keys_loaded.size(); ++i) {
        SMESH_TEST_EQ(keys_loaded[i].first, keys_recreated[i].first);
        SMESH_TEST_EQ(keys_loaded[i].second, keys_recreated[i].second);
    }

    SMESH_TEST_ASSERT(mesh->write(path2) == SMESH_SUCCESS);
    auto mesh2 = Mesh::create_from_file(comm, path2);
    SMESH_TEST_ASSERT(mesh2 != nullptr);
    auto loaded2 = mesh2->sidesets("left");
    SMESH_TEST_EQ(loaded2.size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(loaded2[0]->size(), loaded[0]->size());
    SMESH_TEST_EQ(comm->sum(static_cast<i64>(loaded[0]->size())),
                  comm->sum(static_cast<i64>(recreated[0]->size())));

    auto cloned = mesh->clone();
    SMESH_TEST_ASSERT(cloned != nullptr);
    SMESH_TEST_EQ(cloned->sidesets("left").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(cloned->sidesets("left")[0]->size(), loaded[0]->size());

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        std::filesystem::remove_all(path2.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_mpi_hex8_refine);
    SMESH_RUN_TEST(test_mpi_quad_extrude);
    SMESH_RUN_TEST(test_mpi_clone_convert);
    SMESH_RUN_TEST(test_mpi_ss_derefine);
    SMESH_RUN_TEST(test_mpi_mesh_sideset_io);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

