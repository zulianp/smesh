#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "smesh_distributed_base.hpp"
#include "smesh_edgeset.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
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

static std::shared_ptr<Mesh> create_tri3_square_3d(const ptrdiff_t nx, const ptrdiff_t ny) {
    auto t2 = Mesh::create_tri3_square(Communicator::self(), nx, ny);
    auto p3 = create_host_buffer<geom_t>(3, static_cast<size_t>(t2->n_nodes()));
    std::memcpy(p3->data()[0], t2->points()->data()[0], static_cast<size_t>(t2->n_nodes()) * sizeof(geom_t));
    std::memcpy(p3->data()[1], t2->points()->data()[1], static_cast<size_t>(t2->n_nodes()) * sizeof(geom_t));
    std::memset(p3->data()[2], 0, static_cast<size_t>(t2->n_nodes()) * sizeof(geom_t));
    return std::make_shared<Mesh>(Communicator::self(), t2->blocks(), p3);
}

static std::shared_ptr<Mesh> create_quad4_square_3d(const ptrdiff_t nx, const ptrdiff_t ny) {
    auto q2 = Mesh::create_quad4_square(Communicator::self(), nx, ny);
    auto p3 = create_host_buffer<geom_t>(3, static_cast<size_t>(q2->n_nodes()));
    std::memcpy(p3->data()[0], q2->points()->data()[0], static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
    std::memcpy(p3->data()[1], q2->points()->data()[1], static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
    std::memset(p3->data()[2], 0, static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
    return std::make_shared<Mesh>(Communicator::self(), q2->blocks(), p3);
}

static std::shared_ptr<Mesh> create_edge2_line_3d(const ptrdiff_t n_seg, const enum ElemType et = EDGE2) {
    auto elems = create_host_buffer<idx_t>(2, static_cast<size_t>(n_seg));
    auto pts   = create_host_buffer<geom_t>(3, static_cast<size_t>(n_seg + 1));
    for (ptrdiff_t i = 0; i < n_seg; ++i) {
        elems->data()[0][i] = static_cast<idx_t>(i);
        elems->data()[1][i] = static_cast<idx_t>(i + 1);
    }
    for (ptrdiff_t i = 0; i <= n_seg; ++i) {
        pts->data()[0][i] = static_cast<geom_t>(i);
        pts->data()[1][i] = 0;
        pts->data()[2][i] = 0;
    }
    return std::make_shared<Mesh>(Communicator::self(), et, elems, pts);
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
    SMESH_TEST_ASSERT(refined->sidesets().empty());

    auto coarse_ss = Sideset::create_from_plane(mesh, 1, 0, 0, 0.0);
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    mesh->add_sideset("left", coarse_ss[0]);
    auto ns_buf = create_nodeset_from_sideset(mesh, coarse_ss[0]);
    SMESH_TEST_ASSERT(ns_buf != nullptr);
    mesh->add_nodeset("left_nodes", Nodeset::create(mesh->comm(), ns_buf));
    auto refined_reg = refine(mesh, 1);
    SMESH_TEST_ASSERT(refined_reg != nullptr);
    SMESH_TEST_EQ(refined_reg->sidesets("left").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(refined_reg->sidesets("left")[0]->size(), coarse_ss[0]->size() * 4);
    SMESH_TEST_EQ(refined_reg->nodesets("left_nodes").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(refined_reg->nodesets("left_nodes")[0]->size(),
                 static_cast<ptrdiff_t>(ns_buf->size()));
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined_reg);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), refined_reg->sidesets("left")[0]->size());
    SMESH_TEST_EQ(mapped->block_id(), coarse_ss[0]->block_id());
    for (ptrdiff_t i = 0; i < mapped->size(); ++i) {
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] >= 0);
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] < refined_reg->n_elements(mapped->block_id()));
    }

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_serial_trishell3_refine() {
    auto tri = Mesh::create_tri3_square(Communicator::self(), 2, 2);
    SMESH_TEST_ASSERT(tri != nullptr);
    const ptrdiff_t n_e = tri->n_elements();
    tri->set_element_type(0, TRISHELL3);
    SMESH_TEST_EQ(tri->element_type(0), TRISHELL3);
    auto refined = refine(tri, 1);
    SMESH_TEST_ASSERT(refined != nullptr);
    SMESH_TEST_EQ(refined->element_type(0), TRISHELL3);
    SMESH_TEST_EQ(refined->n_elements(), n_e * 4);
    return SMESH_TEST_SUCCESS;
}

static int test_mpi_tet4_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 10;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_tet_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_tet4_cube(Communicator::self(), nx, ny, nz);
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
    SMESH_TEST_EQ(refined->element_type(0), TET4);
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

static int test_mpi_tri3_refine() {
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
    const Path path = make_tmp_path("smesh_mpi_xf_tri_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_tri3_square_3d(nx, ny);
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
    SMESH_TEST_EQ(refined->element_type(0), TRI3);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 4);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 4);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), coarse_ss[0]->size() * 2);
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

static int test_mpi_trishell3_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 14;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_trishell_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_tri3_square_3d(nx, ny);
        SMESH_TEST_ASSERT(serial != nullptr);
        serial->set_element_type(0, TRISHELL3);
        auto refined = refine(serial, 1);
        SMESH_TEST_ASSERT(refined != nullptr);
        SMESH_TEST_EQ(refined->element_type(0), TRISHELL3);
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
    SMESH_TEST_EQ(mesh->element_type(0), TRISHELL3);

    auto refined = refine(mesh, 1);
    SMESH_TEST_ASSERT(refined != nullptr);
    SMESH_TEST_ASSERT(refined->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), static_cast<int>(mesh->n_blocks()));
    SMESH_TEST_EQ(refined->element_type(0), TRISHELL3);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 4);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 4);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), coarse_ss[0]->size() * 2);
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

static int test_mpi_edge2_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 15;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_edge_refine", token);

    const ptrdiff_t n_seg            = std::max<ptrdiff_t>(4 * comm->size(), 8);
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_edge2_line_3d(n_seg);
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
    SMESH_TEST_EQ(refined->element_type(0), EDGE2);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 2);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 2);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), coarse_ss[0]->size());
    SMESH_TEST_EQ(mapped->block_id(), coarse_ss[0]->block_id());
    for (ptrdiff_t i = 0; i < mapped->size(); ++i) {
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] >= 0);
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] < refined->n_elements(mapped->block_id()));
        SMESH_TEST_EQ(mapped->lfi()->data()[i], coarse_ss[0]->lfi()->data()[i]);
    }

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_quad4_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 12;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_quad_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_quad4_square_3d(nx, ny);
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
    SMESH_TEST_EQ(refined->element_type(0), QUAD4);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 4);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 4);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), coarse_ss[0]->size() * 2);
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

static int test_mpi_wedge6_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 13;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_wedge_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = extrude(create_tri3_square_3d(nx, ny), 1.0, 1);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(serial->element_type(0), WEDGE6);
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
    SMESH_TEST_EQ(refined->element_type(0), WEDGE6);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 8);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
    SMESH_TEST_ASSERT(!coarse_ss.empty());
    auto mapped = map_sideset_through_refine(mesh, coarse_ss[0], refined);
    SMESH_TEST_ASSERT(mapped != nullptr);
    SMESH_TEST_EQ(mapped->size(), coarse_ss[0]->size() * 4);
    SMESH_TEST_EQ(mapped->block_id(), coarse_ss[0]->block_id());
    const enum ElemType parent_side = side_type(WEDGE6, coarse_ss[0]->lfi()->data()[0]);
    for (ptrdiff_t i = 0; i < mapped->size(); ++i) {
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] >= 0);
        SMESH_TEST_ASSERT(mapped->parent()->data()[i] < refined->n_elements(mapped->block_id()));
        SMESH_TEST_EQ(side_type(WEDGE6, mapped->lfi()->data()[i]), parent_side);
    }

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_hex_wedge_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 19;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_hex_wedge_refine", token);

    ptrdiff_t serial_nnodes    = 0;
    ptrdiff_t serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex_dominant_cylinder(Communicator::self(), 1, 1, 1, 8, 2, 0);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(static_cast<int>(serial->n_blocks()), 2);
        SMESH_TEST_EQ(serial->element_type(0), HEX8);
        SMESH_TEST_EQ(serial->element_type(1), WEDGE6);
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
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
    SMESH_TEST_EQ(mesh->element_type(0), HEX8);
    SMESH_TEST_EQ(mesh->element_type(1), WEDGE6);

    auto refined = refine(mesh, 1);
    SMESH_TEST_ASSERT(refined != nullptr);
    SMESH_TEST_ASSERT(refined->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), 2);
    SMESH_TEST_EQ(refined->element_type(0), HEX8);
    SMESH_TEST_EQ(refined->element_type(1), WEDGE6);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(1)->n_elements(), mesh->block(1)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 8);
    SMESH_TEST_EQ(refined->block(1)->n_elements_owned(), mesh->block(1)->n_elements_owned() * 8);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    auto coarse_ss = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
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

static std::shared_ptr<Mesh> create_pyramid_pairs_serial(const ptrdiff_t pairs) {
    if (pairs < 1) {
        return nullptr;
    }

    // Connected strip along x so MPI partitions share nodes (disconnected pairs
    // yield n_ghosts == 0 and abort in mesh_create_parallel).
    const ptrdiff_t n_line  = pairs + 1;
    const ptrdiff_t n_nodes = 2 * n_line + 2 * pairs;
    auto            points  = create_host_buffer<geom_t>(3, static_cast<size_t>(n_nodes));
    auto            elems   = create_host_buffer<idx_t>(5, static_cast<size_t>(2 * pairs));
    auto *const     x       = points->data()[0];
    auto *const     y       = points->data()[1];
    auto *const     z       = points->data()[2];
    for (ptrdiff_t i = 0; i < n_line; ++i) {
        x[i]          = static_cast<geom_t>(i);
        y[i]          = 0;
        z[i]          = 0;
        x[n_line + i] = static_cast<geom_t>(i);
        y[n_line + i] = 1;
        z[n_line + i] = 0;
    }
    const ptrdiff_t apex0 = 2 * n_line;
    for (ptrdiff_t p = 0; p < pairs; ++p) {
        x[apex0 + p]           = static_cast<geom_t>(p) + geom_t(0.5);
        y[apex0 + p]           = geom_t(0.5);
        z[apex0 + p]           = 1;
        x[apex0 + pairs + p]   = static_cast<geom_t>(p) + geom_t(0.5);
        y[apex0 + pairs + p]   = geom_t(0.5);
        z[apex0 + pairs + p]   = -1;
    }
    for (ptrdiff_t p = 0; p < pairs; ++p) {
        const idx_t n00 = static_cast<idx_t>(p);
        const idx_t n10 = static_cast<idx_t>(p + 1);
        const idx_t n01 = static_cast<idx_t>(n_line + p);
        const idx_t n11 = static_cast<idx_t>(n_line + p + 1);
        const idx_t ap  = static_cast<idx_t>(apex0 + p);
        const idx_t am  = static_cast<idx_t>(apex0 + pairs + p);
        const idx_t p0[5] = {n00, n10, n11, n01, ap};
        const idx_t p1[5] = {n00, n01, n11, n10, am};
        const idx_t e0    = static_cast<idx_t>(2 * p);
        for (int d = 0; d < 5; ++d) {
            elems->data()[d][e0 + 0] = p0[d];
            elems->data()[d][e0 + 1] = p1[d];
        }
    }

    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(std::make_shared<Mesh::Block>("pyramid", PYRAMID5, elems));
    return std::make_shared<Mesh>(Communicator::self(), blocks, points);
}

static int test_mpi_pyramid_refine() {
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
    const Path pyr_path = make_tmp_path("smesh_mpi_xf_pyr_refine", token);

    const ptrdiff_t n_pairs = std::max<ptrdiff_t>(2 * comm->size(), 4);
    ptrdiff_t       serial_nnodes    = 0;
    ptrdiff_t       serial_nelements = 0;
    ptrdiff_t       serial_npyr      = 0;
    ptrdiff_t       serial_ntet      = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(pyr_path.to_string());
        auto pyr = create_pyramid_pairs_serial(n_pairs);
        SMESH_TEST_ASSERT(pyr != nullptr);
        SMESH_TEST_EQ(pyr->element_type(0), PYRAMID5);
        auto ref_serial = refine(pyr, 1);
        SMESH_TEST_ASSERT(ref_serial != nullptr);
        SMESH_TEST_EQ(static_cast<int>(ref_serial->n_blocks()), 2);
        SMESH_TEST_EQ(ref_serial->element_type(0), PYRAMID5);
        SMESH_TEST_EQ(ref_serial->element_type(1), TET4);
        SMESH_TEST_EQ(ref_serial->n_elements(0), pyr->n_elements() * (ptrdiff_t)sspyramid_n_pyr(2));
        SMESH_TEST_EQ(ref_serial->n_elements(1), pyr->n_elements() * (ptrdiff_t)sspyramid_n_tet(2));
        serial_nnodes    = ref_serial->n_nodes();
        serial_nelements = ref_serial->n_elements();
        serial_npyr      = ref_serial->n_elements(0);
        serial_ntet      = ref_serial->n_elements(1);
        auto pyr_ss = to_semistructured(2, pyr);
        SMESH_TEST_ASSERT(pyr_ss != nullptr);
        SMESH_TEST_EQ(pyr_ss->element_type(0), semistructured_type(PYRAMID5, 2));
        SMESH_TEST_EQ(pyr_ss->n_nodes(), serial_nnodes);
        SMESH_TEST_ASSERT(pyr->write(pyr_path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes, 1, 0);
    comm->broadcast(&serial_nelements, 1, 0);
    comm->broadcast(&serial_npyr, 1, 0);
    comm->broadcast(&serial_ntet, 1, 0);
    comm->barrier();

    auto pyr = Mesh::create_from_file(comm, pyr_path);
    SMESH_TEST_ASSERT(pyr != nullptr);
    SMESH_TEST_ASSERT(pyr->is_distributed());
    SMESH_TEST_EQ(pyr->element_type(0), PYRAMID5);
    auto refined = refine(pyr, 1);
    SMESH_TEST_ASSERT(refined != nullptr);
    SMESH_TEST_ASSERT(refined->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), 2);
    SMESH_TEST_EQ(refined->element_type(0), PYRAMID5);
    SMESH_TEST_EQ(refined->element_type(1), TET4);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), pyr->block(0)->n_elements() * (ptrdiff_t)sspyramid_n_pyr(2));
    SMESH_TEST_EQ(refined->block(1)->n_elements(), pyr->block(0)->n_elements() * (ptrdiff_t)sspyramid_n_tet(2));
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(),
                 pyr->block(0)->n_elements_owned() * (ptrdiff_t)sspyramid_n_pyr(2));
    SMESH_TEST_EQ(refined->block(1)->n_elements_owned(),
                 pyr->block(0)->n_elements_owned() * (ptrdiff_t)sspyramid_n_tet(2));
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    ptrdiff_t n_pyr_owned = refined->block(0)->n_elements_owned();
    ptrdiff_t n_tet_owned = refined->block(1)->n_elements_owned();
    ptrdiff_t n_pyr_sum   = 0;
    ptrdiff_t n_tet_sum   = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&n_pyr_owned, &n_pyr_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, comm->get()));
    SMESH_MPI_CATCH(MPI_Allreduce(&n_tet_owned, &n_tet_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, comm->get()));
    SMESH_TEST_EQ(n_pyr_sum, serial_npyr);
    SMESH_TEST_EQ(n_tet_sum, serial_ntet);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(pyr_path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static std::shared_ptr<Mesh> repeat_mesh_components(const std::shared_ptr<Mesh> &mesh, const ptrdiff_t copies) {
    if (!mesh || copies < 1) {
        return nullptr;
    }
    const int       dim     = mesh->spatial_dimension();
    const ptrdiff_t n_nodes = mesh->n_nodes();
    const geom_t    dx      = geom_t(2);
    auto            points  = create_host_buffer<geom_t>(static_cast<size_t>(dim), static_cast<size_t>(n_nodes * copies));
    for (ptrdiff_t c = 0; c < copies; ++c) {
        for (int d = 0; d < dim; ++d) {
            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                points->data()[d][c * n_nodes + i] =
                        mesh->points()->data()[d][i] + (d == 0 ? dx * static_cast<geom_t>(c) : geom_t(0));
            }
        }
    }
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        auto            src = mesh->block(static_cast<block_idx_t>(b));
        const int       nxe = src->n_nodes_per_element();
        const ptrdiff_t n0  = src->n_elements();
        auto            dst = create_host_buffer<idx_t>(static_cast<size_t>(nxe), static_cast<size_t>(n0 * copies));
        for (int d = 0; d < nxe; ++d) {
            for (ptrdiff_t c = 0; c < copies; ++c) {
                const idx_t off = static_cast<idx_t>(c * n_nodes);
                for (ptrdiff_t e = 0; e < n0; ++e) {
                    dst->data()[d][c * n0 + e] = src->elements()->data()[d][e] + off;
                }
            }
        }
        blocks.push_back(std::make_shared<Mesh::Block>(src->name(), src->element_type(), dst));
    }
    return std::make_shared<Mesh>(mesh->comm(), blocks, points);
}

static int test_mpi_hex_dominant_refine() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 37;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_xf_hexdom_refine", token);

    ptrdiff_t serial_nnodes    = 0;
    ptrdiff_t serial_nelements = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = repeat_mesh_components(Mesh::create_hex_dominant_serial(Communicator::self()),
                                             std::max<ptrdiff_t>(2 * comm->size(), 4));
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(static_cast<int>(serial->n_blocks()), 4);
        SMESH_TEST_EQ(serial->element_type(0), HEX8);
        SMESH_TEST_EQ(serial->element_type(1), PYRAMID5);
        SMESH_TEST_EQ(serial->element_type(2), TET4);
        SMESH_TEST_EQ(serial->element_type(3), WEDGE6);
        auto refined = refine(serial, 1);
        SMESH_TEST_ASSERT(refined != nullptr);
        SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), 5);
        SMESH_TEST_EQ(refined->element_type(0), HEX8);
        SMESH_TEST_EQ(refined->element_type(1), PYRAMID5);
        SMESH_TEST_EQ(refined->element_type(2), TET4);
        SMESH_TEST_EQ(refined->element_type(3), TET4);
        SMESH_TEST_EQ(refined->element_type(4), WEDGE6);
        SMESH_TEST_EQ(refined->n_elements(0), serial->n_elements(0) * 8);
        SMESH_TEST_EQ(refined->n_elements(1), serial->n_elements(1) * (ptrdiff_t)sspyramid_n_pyr(2));
        SMESH_TEST_EQ(refined->n_elements(2), serial->n_elements(1) * (ptrdiff_t)sspyramid_n_tet(2));
        SMESH_TEST_EQ(refined->n_elements(3), serial->n_elements(2) * 8);
        SMESH_TEST_EQ(refined->n_elements(4), serial->n_elements(3) * 8);
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
    SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 4);

    auto refined = refine(mesh, 1);
    SMESH_TEST_ASSERT(refined != nullptr);
    SMESH_TEST_ASSERT(refined->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), 5);
    SMESH_TEST_EQ(refined->element_type(0), HEX8);
    SMESH_TEST_EQ(refined->element_type(1), PYRAMID5);
    SMESH_TEST_EQ(refined->element_type(2), TET4);
    SMESH_TEST_EQ(refined->element_type(3), TET4);
    SMESH_TEST_EQ(refined->element_type(4), WEDGE6);
    SMESH_TEST_EQ(refined->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(refined->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(refined->block(0)->n_elements(), mesh->block(0)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(1)->n_elements(), mesh->block(1)->n_elements() * (ptrdiff_t)sspyramid_n_pyr(2));
    SMESH_TEST_EQ(refined->block(2)->n_elements(), mesh->block(1)->n_elements() * (ptrdiff_t)sspyramid_n_tet(2));
    SMESH_TEST_EQ(refined->block(3)->n_elements(), mesh->block(2)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(4)->n_elements(), mesh->block(3)->n_elements() * 8);
    SMESH_TEST_EQ(refined->block(0)->n_elements_owned(), mesh->block(0)->n_elements_owned() * 8);
    SMESH_TEST_EQ(refined->block(1)->n_elements_owned(),
                 mesh->block(1)->n_elements_owned() * (ptrdiff_t)sspyramid_n_pyr(2));
    SMESH_TEST_EQ(refined->block(2)->n_elements_owned(),
                 mesh->block(1)->n_elements_owned() * (ptrdiff_t)sspyramid_n_tet(2));
    SMESH_TEST_EQ(refined->block(3)->n_elements_owned(), mesh->block(2)->n_elements_owned() * 8);
    SMESH_TEST_EQ(refined->block(4)->n_elements_owned(), mesh->block(3)->n_elements_owned() * 8);
    SMESH_TEST_EQ(check_owned_gids_unique(*refined), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_higher_order_ss_refine_rejected() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 41;
    }
    comm->broadcast(&token, 1, 0);
    const Path hex_path   = make_tmp_path("smesh_mpi_xf_ss_refine", token);
    const Path tet10_path = make_tmp_path("smesh_mpi_xf_tet10_refine", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    if (comm->rank() == 0) {
        std::filesystem::remove_all(hex_path.to_string());
        std::filesystem::remove_all(tet10_path.to_string());
        auto hex = Mesh::create_hex8_cube(Communicator::self(), nx, 2, 2);
        SMESH_TEST_ASSERT(hex != nullptr);
        auto ss = to_semistructured(2, hex);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_ASSERT(refine(ss, 1) == nullptr);
        SMESH_TEST_ASSERT(hex->write(hex_path) == SMESH_SUCCESS);

        auto tet   = Mesh::create_tet4_cube(Communicator::self(), nx, 2, 2);
        auto tet10 = promote_to(TET10, tet);
        SMESH_TEST_ASSERT(tet10 != nullptr);
        SMESH_TEST_EQ(tet10->element_type(0), TET10);
        SMESH_TEST_ASSERT(refine(tet10, 1) == nullptr);
        SMESH_TEST_ASSERT(tet10->write(tet10_path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto hex = Mesh::create_from_file(comm, hex_path);
    SMESH_TEST_ASSERT(hex != nullptr);
    auto ss = to_semistructured(2, hex);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_ASSERT(is_semistructured_type(ss->element_type(0)));
    SMESH_TEST_ASSERT(refine(ss, 1) == nullptr);

    auto tet10 = Mesh::create_from_file(comm, tet10_path);
    SMESH_TEST_ASSERT(tet10 != nullptr);
    SMESH_TEST_ASSERT(tet10->is_distributed());
    SMESH_TEST_EQ(tet10->element_type(0), TET10);
    SMESH_TEST_ASSERT(refine(tet10, 1) == nullptr);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(hex_path.to_string());
        std::filesystem::remove_all(tet10_path.to_string());
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

static int test_mpi_write_with_xdmf() {
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
    const Path path  = make_tmp_path("smesh_mpi_write_xdmf", token);
    const Path pathx = make_tmp_path("smesh_mpi_write_xdmf_out", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        std::filesystem::remove_all(pathx.to_string());
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
    SMESH_TEST_ASSERT(mesh->write_with_xdmf(pathx) == SMESH_SUCCESS);

    if (comm->rank() == 0) {
        SMESH_TEST_ASSERT((pathx / "mesh.xdmf").exists());
        const Path aos =
                pathx / (std::string("connectivity.") + std::string(TypeToString<idx_t>::value()));
        SMESH_TEST_ASSERT(aos.exists());
        FILE *fp = fopen(aos.c_str(), "rb");
        SMESH_TEST_ASSERT(fp != nullptr);
        fseek(fp, 0, SEEK_END);
        const long nbytes = ftell(fp);
        fclose(fp);
        SMESH_TEST_EQ(nbytes,
                      (long)mesh->distributed()->n_elements_global() * 8L * (long)sizeof(idx_t));
        const Path surf = pathx / "sidesets" / "left" /
                          (std::string("surface.") + std::string(TypeToString<idx_t>::value()));
        SMESH_TEST_ASSERT(surf.exists());
        std::ifstream xdmf((pathx / "mesh.xdmf").to_string());
        std::string xml((std::istreambuf_iterator<char>(xdmf)), std::istreambuf_iterator<char>());
        SMESH_TEST_ASSERT(xml.find("Hexahedron") != std::string::npos);
        SMESH_TEST_ASSERT(xml.find("Quadrilateral") != std::string::npos);
        std::filesystem::remove_all(path.to_string());
        std::filesystem::remove_all(pathx.to_string());
    }
    comm->barrier();
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_mesh_edgeset_nodeset_io() {
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
    const Path path = make_tmp_path("smesh_mpi_mesh_en", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, 2, 2);
        SMESH_TEST_ASSERT(serial != nullptr);
        const ptrdiff_t n_e = serial->n_elements(0);
        auto parent = create_host_buffer<element_idx_t>((size_t)n_e);
        auto lei    = create_host_buffer<i16>((size_t)n_e);
        for (ptrdiff_t i = 0; i < n_e; ++i) {
            parent->data()[i] = (element_idx_t)i;
            lei->data()[i]    = 0;
        }
        serial->add_edgeset("e0", Edgeset::create(serial->comm(), parent, lei, 0));
        auto left = Sideset::create_from_selector(
                serial, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
        auto ns_buf = create_nodeset_from_sideset(serial, left[0]);
        serial->add_nodeset("left_nodes", Nodeset::create(serial->comm(), ns_buf));
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto loaded_e = mesh->edgesets("e0");
    SMESH_TEST_EQ(loaded_e.size(), static_cast<size_t>(1));
    SMESH_TEST_ASSERT(comm->sum(static_cast<i64>(loaded_e[0]->size())) > 0);

    auto loaded_n = mesh->nodesets("left_nodes");
    SMESH_TEST_EQ(loaded_n.size(), static_cast<size_t>(1));
    SMESH_TEST_ASSERT(comm->sum(static_cast<i64>(loaded_n[0]->size())) > 0);

    auto cloned = mesh->clone();
    SMESH_TEST_EQ(cloned->edgesets("e0").size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(cloned->nodesets("left_nodes").size(), static_cast<size_t>(1));

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_serial_trishell3_refine);
    SMESH_RUN_TEST(test_mpi_hex8_refine);
    SMESH_RUN_TEST(test_mpi_tet4_refine);
    SMESH_RUN_TEST(test_mpi_tri3_refine);
    SMESH_RUN_TEST(test_mpi_trishell3_refine);
    SMESH_RUN_TEST(test_mpi_edge2_refine);
    SMESH_RUN_TEST(test_mpi_quad4_refine);
    SMESH_RUN_TEST(test_mpi_wedge6_refine);
    SMESH_RUN_TEST(test_mpi_hex_wedge_refine);
    SMESH_RUN_TEST(test_mpi_pyramid_refine);
    SMESH_RUN_TEST(test_mpi_hex_dominant_refine);
    SMESH_RUN_TEST(test_mpi_higher_order_ss_refine_rejected);
    SMESH_RUN_TEST(test_mpi_quad_extrude);
    SMESH_RUN_TEST(test_mpi_clone_convert);
    SMESH_RUN_TEST(test_mpi_ss_derefine);
    SMESH_RUN_TEST(test_mpi_mesh_sideset_io);
    SMESH_RUN_TEST(test_mpi_mesh_edgeset_nodeset_io);
    SMESH_RUN_TEST(test_mpi_write_with_xdmf);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

