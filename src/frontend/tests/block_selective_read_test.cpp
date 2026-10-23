// Mesh::read(path, block_names): only the named blocks, only the nodes they reference, each node
// mapped to its id in the file; sets restricted to what was read. Serial, and with more than one
// rank the distributed read against the serial one.

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <vector>

#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

using namespace smesh;

namespace {
    // A two-block (white/black) checkerboard with a sideset on black and a nodeset, written to a
    // fresh folder on rank 0; returns the folder.
    Path write_checkerboard(const std::shared_ptr<Communicator> &comm, const char *tag) {
        int token = 0;
        if (comm->rank() == 0) {
            token = static_cast<int>(std::time(nullptr));
        }
        comm->broadcast(&token, 1, 0);
        char buffer[256];
        std::snprintf(buffer, sizeof(buffer), "/tmp/smesh_block_selective_%s_%d_%d", tag, comm->size(), token);
        const Path path(buffer);

        if (comm->rank() == 0) {
            std::filesystem::remove_all(path.to_string());
            auto mesh = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 4, 4, 2, 0, 0, 0, 1, 1, 1);

            auto parent = create_host_buffer<element_idx_t>(2);
            auto lfi    = create_host_buffer<i16>(2);
            parent->data()[0] = 0;
            parent->data()[1] = 1;
            lfi->data()[0]    = 0;
            lfi->data()[1]    = 3;
            mesh->add_sideset("black_faces", Sideset::create(Communicator::self(), parent, lfi, 1));

            auto nodes = create_host_buffer<idx_t>(3);
            nodes->data()[0] = 0;
            nodes->data()[1] = static_cast<idx_t>(mesh->n_nodes() / 2);
            nodes->data()[2] = static_cast<idx_t>(mesh->n_nodes() - 1);
            mesh->add_nodeset("probe", Nodeset::create(Communicator::self(), nodes));

            mesh->write(path);
        }
        comm->barrier();
        return path;
    }

    // The sorted file ids of the nodes `block` of `full` references.
    std::vector<idx_t> referenced_nodes(const Mesh &full, const std::string &block) {
        auto            b    = full.find_block(block);
        idx_t **const   e    = b->elements()->data();
        std::vector<idx_t> ids;
        for (int v = 0; v < b->n_nodes_per_element(); ++v) {
            ids.insert(ids.end(), e[v], e[v] + b->n_elements());
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        return ids;
    }
}  // namespace

int test_serial_block_selection() {
    auto       self = Communicator::self();
    const Path path = write_checkerboard(self, "serial");

    Mesh full(self);
    SMESH_TEST_ASSERT(full.read(path) == SMESH_SUCCESS);
    SMESH_TEST_EQ(static_cast<int>(full.n_blocks()), 2);
    SMESH_TEST_ASSERT(full.node_mapping() == nullptr);

    for (const char *name : {"white", "black"}) {
        Mesh part(self);
        SMESH_TEST_ASSERT(part.read(path, {name}) == SMESH_SUCCESS);
        SMESH_TEST_EQ(static_cast<int>(part.n_blocks()), 1);
        SMESH_TEST_ASSERT(part.block(0)->name() == name);
        SMESH_TEST_EQ(part.n_elements(), full.find_block(name)->n_elements());

        // Only the referenced nodes, in file order, each with its file id.
        const auto ids = referenced_nodes(full, name);
        SMESH_TEST_EQ(part.n_nodes(), static_cast<ptrdiff_t>(ids.size()));
        SMESH_TEST_ASSERT(part.n_nodes() < full.n_nodes());
        auto map = part.node_mapping();
        SMESH_TEST_ASSERT(map != nullptr);
        for (ptrdiff_t i = 0; i < part.n_nodes(); ++i) {
            SMESH_TEST_EQ(map->data()[i], ids[i]);
            for (int d = 0; d < 3; ++d) {
                SMESH_TEST_ASSERT(part.points()->data()[d][i] == full.points()->data()[d][map->data()[i]]);
            }
        }
        idx_t **const pe = part.block(0)->elements()->data();
        idx_t **const fe = full.find_block(name)->elements()->data();
        for (int v = 0; v < 8; ++v) {
            for (ptrdiff_t e = 0; e < part.n_elements(); ++e) {
                SMESH_TEST_EQ(map->data()[pe[v][e]], fe[v][e]);
            }
        }

        // The black sideset is kept, as block 0, only with black; the nodeset keeps the probe
        // nodes the block references, renumbered.
        const bool black = std::string(name) == "black";
        SMESH_TEST_EQ(static_cast<int>(part.sidesets("black_faces").size()), black ? 1 : 0);
        if (black) {
            auto ss = part.sidesets("black_faces")[0];
            SMESH_TEST_EQ(static_cast<int>(ss->block_id()), 0);
            SMESH_TEST_EQ(ss->size(), static_cast<ptrdiff_t>(2));
        }
        auto            probe      = full.nodesets("probe")[0];
        ptrdiff_t       n_expected = 0;
        for (ptrdiff_t i = 0; i < probe->size(); ++i) {
            n_expected += std::binary_search(ids.begin(), ids.end(), probe->nodes()->data()[i]);
        }
        auto kept = part.nodesets("probe")[0];
        SMESH_TEST_EQ(kept->size(), n_expected);
        for (ptrdiff_t i = 0; i < kept->size(); ++i) {
            const idx_t file_id = map->data()[kept->nodes()->data()[i]];
            SMESH_TEST_ASSERT(std::find(probe->nodes()->data(), probe->nodes()->data() + probe->size(), file_id) !=
                              probe->nodes()->data() + probe->size());
        }
    }

    // Both blocks is the whole mesh, untouched.
    Mesh both(self);
    SMESH_TEST_ASSERT(both.read(path, {"black", "white"}) == SMESH_SUCCESS);
    SMESH_TEST_EQ(both.n_nodes(), full.n_nodes());
    SMESH_TEST_ASSERT(both.node_mapping() == nullptr);

    std::filesystem::remove_all(path.to_string());
    return SMESH_TEST_SUCCESS;
}

int test_distributed_block_selection() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }
    const Path path = write_checkerboard(comm, "mpi");

    // The serial selection is the reference.
    Mesh serial(Communicator::self());
    SMESH_TEST_ASSERT(serial.read(path, {"black"}) == SMESH_SUCCESS);
    const auto serial_map = serial.node_mapping();

    Mesh mesh(comm);
    SMESH_TEST_ASSERT(mesh.read(path, {"black"}) == SMESH_SUCCESS);
    SMESH_TEST_EQ(static_cast<int>(mesh.n_blocks()), 1);
    auto dist = mesh.distributed();
    SMESH_TEST_ASSERT(dist != nullptr);

    // Owned nodes over all ranks are the serial selection's nodes, by file id.
    const ptrdiff_t n_owned = dist->n_nodes_owned();
    SMESH_TEST_EQ(comm->sum(n_owned), serial.n_nodes());
    std::vector<long long> owned_ids(n_owned);
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
        owned_ids[i] = static_cast<long long>(dist->node_mapping()->data()[i]);
    }
    std::vector<int> counts(comm->size()), displs(comm->size() + 1, 0);
    const int        n_owned_int = static_cast<int>(n_owned);
    MPI_Allgather(&n_owned_int, 1, MPI_INT, counts.data(), 1, MPI_INT, comm->get());
    for (int r = 0; r < comm->size(); ++r) {
        displs[r + 1] = displs[r] + counts[r];
    }
    std::vector<long long> all_ids(displs[comm->size()]);
    MPI_Allgatherv(owned_ids.data(), n_owned_int, MPI_LONG_LONG, all_ids.data(), counts.data(), displs.data(),
                   MPI_LONG_LONG, comm->get());
    std::sort(all_ids.begin(), all_ids.end());
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(all_ids.size()), serial.n_nodes());
    for (ptrdiff_t i = 0; i < serial.n_nodes(); ++i) {
        SMESH_TEST_EQ(all_ids[i], static_cast<long long>(serial_map->data()[i]));
    }

    // Every local node is a node of a local element: none of white's is left behind.
    std::vector<char> used(mesh.n_nodes(), 0);
    auto              block = mesh.block(0);
    for (int v = 0; v < block->n_nodes_per_element(); ++v) {
        for (ptrdiff_t e = 0; e < block->n_elements(); ++e) {
            used[block->elements()->data()[v][e]] = 1;
        }
    }
    SMESH_TEST_ASSERT(std::all_of(used.begin(), used.end(), [](char u) { return u == 1; }));
    SMESH_TEST_EQ(comm->sum(block->n_elements_owned()), serial.n_elements());

    comm->barrier();
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    if (Communicator::world()->rank() == 0) {
        SMESH_RUN_TEST(test_serial_block_selection);
    }
    SMESH_RUN_TEST(test_distributed_block_selection);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
