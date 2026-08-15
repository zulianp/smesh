#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <vector>

#include "smesh_distributed_base.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_test.hpp"

using namespace smesh;

using NodeKey = std::array<long long, 3>;
using FaceKey = std::vector<NodeKey>;

static Path make_tmp_path(const char *prefix, const int token) {
    auto comm = Communicator::world();
    char buf[256];
    std::snprintf(buf, sizeof(buf), "/tmp/%s_%d_%d", prefix, comm->size(), token);
    return Path(buf);
}

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

static NodeKey quantize_node(const geom_t *const *pts, const idx_t node) {
    return {std::llround(static_cast<double>(pts[0][node]) * 1e9),
            std::llround(static_cast<double>(pts[1][node]) * 1e9),
            std::llround(static_cast<double>(pts[2][node]) * 1e9)};
}

static int collect_face_keys(const std::shared_ptr<Mesh>                 &mesh,
                             const std::vector<std::shared_ptr<Sideset>> &sidesets,
                             std::vector<FaceKey>                        *faces,
                             enum ElemType                               *type_out,
                             int                                         *nnxs_out) {
    faces->clear();
    if (type_out) {
        *type_out = INVALID;
    }
    if (nnxs_out) {
        *nnxs_out = 0;
    }
    auto pts = mesh->points()->data();
    for (const auto &ss : sidesets) {
        if (!ss) {
            continue;
        }
        const ptrdiff_t n_e = mesh->n_elements(ss->block_id());
        for (ptrdiff_t i = 0; i < ss->size(); ++i) {
            SMESH_TEST_ASSERT(ss->parent()->data()[i] >= 0);
            SMESH_TEST_ASSERT(ss->parent()->data()[i] < n_e);
        }
        auto [st, surface] = create_surface_from_sideset(mesh, ss);
        SMESH_TEST_ASSERT(surface != nullptr);
        if (type_out && *type_out == INVALID) {
            *type_out = st;
        }
        const int nnxs = static_cast<int>(surface->extent(0));
        if (nnxs_out && *nnxs_out == 0) {
            *nnxs_out = nnxs;
        }
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(surface->extent(1)), ss->size());
        for (ptrdiff_t e = 0; e < static_cast<ptrdiff_t>(surface->extent(1)); ++e) {
            FaceKey key(static_cast<size_t>(nnxs));
            for (int ln = 0; ln < nnxs; ++ln) {
                const idx_t node = surface->data()[ln][e];
                SMESH_TEST_ASSERT(node >= 0);
                SMESH_TEST_ASSERT(node < mesh->n_nodes());
                key[static_cast<size_t>(ln)] = quantize_node(pts, node);
            }
            std::sort(key.begin(), key.end());
            faces->push_back(std::move(key));
        }
    }
    return SMESH_TEST_SUCCESS;
}

#ifdef SMESH_ENABLE_MPI
static int match_faces_mpi(const std::vector<FaceKey> &serial_faces, const std::vector<FaceKey> &local_faces,
                           MPI_Comm comm) {
    std::vector<int> claimed(serial_faces.size(), 0);
    for (const auto &face : local_faces) {
        bool found = false;
        for (size_t j = 0; j < serial_faces.size(); ++j) {
            if (claimed[j]) {
                continue;
            }
            if (face == serial_faces[j]) {
                claimed[j] = 1;
                found      = true;
                break;
            }
        }
        SMESH_TEST_ASSERT(found);
    }
    if (!claimed.empty()) {
        MPI_Allreduce(MPI_IN_PLACE, claimed.data(), static_cast<int>(claimed.size()), MPI_INT, MPI_SUM, comm);
        for (size_t j = 0; j < claimed.size(); ++j) {
            SMESH_TEST_EQ(claimed[j], 1);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static void bcast_faces(std::vector<FaceKey> *faces, int *nnxs, const int root) {
    auto comm    = Communicator::world();
    int  n_faces = static_cast<int>(faces->size());
    comm->broadcast(&n_faces, 1, root);
    comm->broadcast(nnxs, 1, root);
    faces->resize(static_cast<size_t>(n_faces));
    if (n_faces == 0 || *nnxs == 0) {
        return;
    }
    std::vector<i64> flat(static_cast<size_t>(n_faces) * static_cast<size_t>(*nnxs) * 3);
    if (comm->rank() == root) {
        size_t w = 0;
        for (int f = 0; f < n_faces; ++f) {
            for (int ln = 0; ln < *nnxs; ++ln) {
                flat[w++] = (*faces)[static_cast<size_t>(f)][static_cast<size_t>(ln)][0];
                flat[w++] = (*faces)[static_cast<size_t>(f)][static_cast<size_t>(ln)][1];
                flat[w++] = (*faces)[static_cast<size_t>(f)][static_cast<size_t>(ln)][2];
            }
        }
    }
    comm->broadcast(flat.data(), static_cast<int>(flat.size()), root);
    if (comm->rank() != root) {
        size_t r = 0;
        for (int f = 0; f < n_faces; ++f) {
            (*faces)[static_cast<size_t>(f)].assign(static_cast<size_t>(*nnxs), NodeKey{});
            for (int ln = 0; ln < *nnxs; ++ln) {
                (*faces)[static_cast<size_t>(f)][static_cast<size_t>(ln)][0] = flat[r++];
                (*faces)[static_cast<size_t>(f)][static_cast<size_t>(ln)][1] = flat[r++];
                (*faces)[static_cast<size_t>(f)][static_cast<size_t>(ln)][2] = flat[r++];
            }
        }
    }
}
#endif

static int test_mpi_hex8_ss_sideset_surface() {
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
    const Path path = make_tmp_path("smesh_mpi_ss_sideset_hex", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;

    std::vector<FaceKey> serial_faces;
    enum ElemType        serial_st = INVALID;
    int                  nnxs      = 0;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = Mesh::create_hex8_checkerboard_cube(Communicator::self(), nx, ny, nz);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(static_cast<int>(serial->n_blocks()), 2);
        auto sidesets = Sideset::create_from_plane(serial, 1, 0, 0, 0.0);
        auto ss       = to_semistructured(2, serial);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_EQ(collect_face_keys(ss, sidesets, &serial_faces, &serial_st, &nnxs), SMESH_TEST_SUCCESS);
        SMESH_TEST_EQ(serial_st, PROTEUS_QUADSHELL9);
        SMESH_TEST_EQ(nnxs, 9);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    {
        int st_i = static_cast<int>(serial_st);
        comm->broadcast(&st_i, 1, 0);
        serial_st = static_cast<enum ElemType>(st_i);
    }
    bcast_faces(&serial_faces, &nnxs, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());
    auto local_ss = Sideset::create_from_plane(mesh, 1, 0, 0, 0.0);
    auto ss       = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 2);

    enum ElemType local_st   = INVALID;
    int           local_nnxs = 0;
    std::vector<FaceKey> local_faces;
    SMESH_TEST_EQ(collect_face_keys(ss, local_ss, &local_faces, &local_st, &local_nnxs), SMESH_TEST_SUCCESS);
    if (!local_faces.empty()) {
        SMESH_TEST_EQ(local_st, PROTEUS_QUADSHELL9);
        SMESH_TEST_EQ(local_nnxs, 9);
    }
    SMESH_TEST_EQ(match_faces_mpi(serial_faces, local_faces, comm->get()), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_tet4_ss_sideset_surface() {
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
    const Path path = make_tmp_path("smesh_mpi_ss_sideset_tet", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);

    std::vector<FaceKey> serial_faces;
    enum ElemType        serial_st = INVALID;
    int                  nnxs      = 0;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto single = Mesh::create_tet4_cube(Communicator::self(), nx, 2, 2);
        auto serial = split_first_half(single);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(static_cast<int>(serial->n_blocks()), 2);
        auto sidesets = Sideset::create_from_plane(serial, 1, 0, 0, 0.0);
        auto ss       = to_semistructured(2, serial);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_EQ(collect_face_keys(ss, sidesets, &serial_faces, &serial_st, &nnxs), SMESH_TEST_SUCCESS);
        SMESH_TEST_EQ(serial_st, TRISHELL6);
        SMESH_TEST_EQ(nnxs, sstet4_n_tri(2));
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    {
        int st_i = static_cast<int>(serial_st);
        comm->broadcast(&st_i, 1, 0);
        serial_st = static_cast<enum ElemType>(st_i);
    }
    bcast_faces(&serial_faces, &nnxs, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto local_ss = Sideset::create_from_plane(mesh, 1, 0, 0, 0.0);
    auto ss       = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 2);

    enum ElemType local_st   = INVALID;
    int           local_nnxs = 0;
    std::vector<FaceKey> local_faces;
    SMESH_TEST_EQ(collect_face_keys(ss, local_ss, &local_faces, &local_st, &local_nnxs), SMESH_TEST_SUCCESS);
    if (!local_faces.empty()) {
        SMESH_TEST_EQ(local_st, TRISHELL6);
        SMESH_TEST_EQ(local_nnxs, 6);
    }
    SMESH_TEST_EQ(match_faces_mpi(serial_faces, local_faces, comm->get()), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_mpi_hex8_ss_sideset_surface);
    SMESH_RUN_TEST(test_mpi_tet4_ss_sideset_surface);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}



