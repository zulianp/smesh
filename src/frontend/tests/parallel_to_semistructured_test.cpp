#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <map>
#include <vector>

#include "smesh_conversion.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_ssquad4.hpp"
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

static int check_owned_gid_prefix(const Mesh &ss, const ptrdiff_t n_prefix) {
#ifndef SMESH_ENABLE_MPI
    SMESH_UNUSED(ss);
    SMESH_UNUSED(n_prefix);
    return SMESH_TEST_SUCCESS;
#else
    SMESH_TEST_ASSERT(ss.is_distributed());
    auto                dist    = ss.distributed();
    const ptrdiff_t     n_owned = dist->n_nodes_owned();
    const large_idx_t  *map     = dist->node_mapping()->data();
    ptrdiff_t           local   = 0;
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
        if (map[i] < (large_idx_t)n_prefix) {
            local++;
        }
    }
    ptrdiff_t sum = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&local, &sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, ss.comm()->get()));
    SMESH_TEST_EQ(sum, n_prefix);
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

    auto ss_h = to_semistructured(2, mesh, true, false);
    SMESH_TEST_ASSERT(ss_h != nullptr);
    SMESH_TEST_ASSERT(ss_h->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss_h->n_blocks()), 1);
    SMESH_TEST_EQ(ss_h->element_type(0), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), ss->distributed()->n_nodes_global());
    SMESH_TEST_EQ(ss_h->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, mesh->distributed()->n_nodes_global()), SMESH_TEST_SUCCESS);

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

static int test_mpi_mixed_hex_tet() {
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
    const Path mixed_path  = make_tmp_path("smesh_mpi_ss_mixed", token);
    const Path hexdom_path = make_tmp_path("smesh_mpi_ss_hexdom", token + 1);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;

    ptrdiff_t serial_nnodes_l2 = 0;
    ptrdiff_t serial_nnodes_l4 = 0;
    ptrdiff_t serial_ncoarse   = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(mixed_path.to_string());
        std::filesystem::remove_all(hexdom_path.to_string());
        auto mixed = create_hex8_tet4_serial(nx, ny, nz);
        SMESH_TEST_ASSERT(mixed != nullptr);
        serial_ncoarse = mixed->n_nodes();
        auto ss2       = to_semistructured(2, mixed);
        auto ss4       = to_semistructured(4, mixed);
        SMESH_TEST_ASSERT(ss2 != nullptr);
        SMESH_TEST_ASSERT(ss4 != nullptr);
        serial_nnodes_l2 = ss2->n_nodes();
        serial_nnodes_l4 = ss4->n_nodes();
        SMESH_TEST_ASSERT(mixed->write(mixed_path) == SMESH_SUCCESS);
        auto hexdom = Mesh::create_hex_dominant_serial(Communicator::self());
        SMESH_TEST_ASSERT(hexdom != nullptr);
        hexdom = repeat_block_elements(hexdom, std::max<ptrdiff_t>(2 * comm->size(), 4));
        SMESH_TEST_ASSERT(hexdom != nullptr);
        SMESH_TEST_ASSERT(hexdom->write(hexdom_path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes_l2, 1, 0);
    comm->broadcast(&serial_nnodes_l4, 1, 0);
    comm->broadcast(&serial_ncoarse, 1, 0);
    comm->barrier();

    auto mixed = Mesh::create_from_file(comm, mixed_path);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_ASSERT(mixed->is_distributed());

    for (int L : {2, 4}) {
        auto ss = to_semistructured(L, mixed);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_ASSERT(ss->is_distributed());
        SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 2);
        SMESH_TEST_ASSERT(ss->block(0)->name() == "hex");
        SMESH_TEST_ASSERT(ss->block(1)->name() == "tet");
        SMESH_TEST_EQ(ss->element_type(0), semistructured_type(HEX8, L));
        SMESH_TEST_EQ(ss->element_type(1), semistructured_type(TET4, L));
        SMESH_TEST_EQ(check_block_layout_copied(*mixed, *ss), SMESH_TEST_SUCCESS);
        SMESH_TEST_EQ(check_owned_gids_unique(*ss), SMESH_TEST_SUCCESS);
        const ptrdiff_t expect = (L == 2) ? serial_nnodes_l2 : serial_nnodes_l4;
        SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), expect);

        const large_idx_t *cnmap = mixed->distributed()->node_mapping()->data();
        const large_idx_t *snmap = ss->distributed()->node_mapping()->data();
        auto               hex_c = mixed->elements(0)->data();
        auto               tet_c = mixed->elements(1)->data();
        auto               hex_s = ss->elements(0)->data();
        auto               tet_s = ss->elements(1)->data();
        const ptrdiff_t    n_hex_owned = mixed->block(0)->n_elements_owned();
        const ptrdiff_t    n_tet_owned = mixed->block(1)->n_elements_owned();
        const int          n_int       = L > 1 ? (L - 1) : 0;

        std::map<std::pair<large_idx_t, large_idx_t>, std::vector<large_idx_t>> hex_gids;
        std::map<std::pair<large_idx_t, large_idx_t>, std::vector<large_idx_t>> tet_gids;
        static const int hex_edges[12][2] = {{0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3},
                                             {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};
        static const int hex_xyz[8][3]    = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        static const int tet_xyz[4][3]    = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (ptrdiff_t e = 0; e < n_hex_owned; ++e) {
            for (int k = 0; k < 12; ++k) {
                const idx_t       a  = hex_c[hex_edges[k][0]][e];
                const idx_t       b  = hex_c[hex_edges[k][1]][e];
                const large_idx_t ga = cnmap[a];
                const large_idx_t gb = cnmap[b];
                const large_idx_t lo = std::min(ga, gb);
                const large_idx_t hi = std::max(ga, gb);
                const int         c0 = (ga == lo) ? hex_edges[k][0] : hex_edges[k][1];
                const int         c1 = (ga == lo) ? hex_edges[k][1] : hex_edges[k][0];
                const int         P1[3] = {hex_xyz[c0][0] * L, hex_xyz[c0][1] * L, hex_xyz[c0][2] * L};
                const int         P2[3] = {hex_xyz[c1][0] * L, hex_xyz[c1][1] * L, hex_xyz[c1][2] * L};
                std::vector<large_idx_t> ids(static_cast<size_t>(n_int));
                for (int t = 1; t < L; ++t) {
                    const int xi = (P1[0] * (L - t) + P2[0] * t) / L;
                    const int yi = (P1[1] * (L - t) + P2[1] * t) / L;
                    const int zi = (P1[2] * (L - t) + P2[2] * t) / L;
                    ids[static_cast<size_t>(t - 1)] = snmap[hex_s[sshex8_lidx(L, xi, yi, zi)][e]];
                }
                hex_gids[std::make_pair(lo, hi)] = ids;
            }
        }
        for (ptrdiff_t e = 0; e < n_tet_owned; ++e) {
            for (int i = 0; i < 4; ++i) {
                for (int j = i + 1; j < 4; ++j) {
                    const large_idx_t ga = cnmap[tet_c[i][e]];
                    const large_idx_t gb = cnmap[tet_c[j][e]];
                    const large_idx_t lo = std::min(ga, gb);
                    const large_idx_t hi = std::max(ga, gb);
                    const int         c0 = (ga == lo) ? i : j;
                    const int         c1 = (ga == lo) ? j : i;
                    const int         P1[3] = {tet_xyz[c0][0] * L, tet_xyz[c0][1] * L, tet_xyz[c0][2] * L};
                    const int         P2[3] = {tet_xyz[c1][0] * L, tet_xyz[c1][1] * L, tet_xyz[c1][2] * L};
                    std::vector<large_idx_t> ids(static_cast<size_t>(n_int));
                    for (int t = 1; t < L; ++t) {
                        const int xi = (P1[0] * (L - t) + P2[0] * t) / L;
                        const int yi = (P1[1] * (L - t) + P2[1] * t) / L;
                        const int zi = (P1[2] * (L - t) + P2[2] * t) / L;
                        ids[static_cast<size_t>(t - 1)] = snmap[tet_s[sstet4_lidx(L, xi, yi, zi)][e]];
                    }
                    tet_gids[std::make_pair(lo, hi)] = ids;
                }
            }
        }

        const int rec = 2 + n_int;
        auto pack_map = [&](const std::map<std::pair<large_idx_t, large_idx_t>, std::vector<large_idx_t>> &m) {
            std::vector<large_idx_t> buf;
            buf.reserve(m.size() * static_cast<size_t>(rec));
            for (const auto &kv : m) {
                buf.push_back(kv.first.first);
                buf.push_back(kv.first.second);
                for (int t = 0; t < n_int; ++t) {
                    buf.push_back(kv.second[static_cast<size_t>(t)]);
                }
            }
            return buf;
        };
        auto hex_buf = pack_map(hex_gids);
        auto tet_buf = pack_map(tet_gids);
        int  n_hex   = static_cast<int>(hex_buf.size());
        int  n_tet   = static_cast<int>(tet_buf.size());
        std::vector<int> hex_counts(static_cast<size_t>(comm->size()));
        std::vector<int> tet_counts(static_cast<size_t>(comm->size()));
        SMESH_MPI_CATCH(MPI_Gather(&n_hex, 1, MPI_INT, hex_counts.data(), 1, MPI_INT, 0, comm->get()));
        SMESH_MPI_CATCH(MPI_Gather(&n_tet, 1, MPI_INT, tet_counts.data(), 1, MPI_INT, 0, comm->get()));
        std::vector<int> hex_displs(static_cast<size_t>(comm->size()));
        std::vector<int> tet_displs(static_cast<size_t>(comm->size()));
        int              hex_total = 0;
        int              tet_total = 0;
        if (comm->rank() == 0) {
            for (int r = 0; r < comm->size(); ++r) {
                hex_displs[static_cast<size_t>(r)] = hex_total;
                tet_displs[static_cast<size_t>(r)] = tet_total;
                hex_total += hex_counts[static_cast<size_t>(r)];
                tet_total += tet_counts[static_cast<size_t>(r)];
            }
        }
        std::vector<large_idx_t> hex_all(static_cast<size_t>(hex_total));
        std::vector<large_idx_t> tet_all(static_cast<size_t>(tet_total));
        SMESH_MPI_CATCH(MPI_Gatherv(hex_buf.data(),
                                    n_hex,
                                    mpi_type<large_idx_t>(),
                                    hex_all.data(),
                                    hex_counts.data(),
                                    hex_displs.data(),
                                    mpi_type<large_idx_t>(),
                                    0,
                                    comm->get()));
        SMESH_MPI_CATCH(MPI_Gatherv(tet_buf.data(),
                                    n_tet,
                                    mpi_type<large_idx_t>(),
                                    tet_all.data(),
                                    tet_counts.data(),
                                    tet_displs.data(),
                                    mpi_type<large_idx_t>(),
                                    0,
                                    comm->get()));
        if (comm->rank() == 0) {
            std::map<std::pair<large_idx_t, large_idx_t>, std::vector<large_idx_t>> hex_all_map;
            std::map<std::pair<large_idx_t, large_idx_t>, std::vector<large_idx_t>> tet_all_map;
            for (int i = 0; i < hex_total; i += rec) {
                std::vector<large_idx_t> ids(static_cast<size_t>(n_int));
                for (int t = 0; t < n_int; ++t) {
                    ids[static_cast<size_t>(t)] = hex_all[static_cast<size_t>(i + 2 + t)];
                }
                hex_all_map[std::make_pair(hex_all[static_cast<size_t>(i)], hex_all[static_cast<size_t>(i + 1)])] = ids;
            }
            for (int i = 0; i < tet_total; i += rec) {
                std::vector<large_idx_t> ids(static_cast<size_t>(n_int));
                for (int t = 0; t < n_int; ++t) {
                    ids[static_cast<size_t>(t)] = tet_all[static_cast<size_t>(i + 2 + t)];
                }
                tet_all_map[std::make_pair(tet_all[static_cast<size_t>(i)], tet_all[static_cast<size_t>(i + 1)])] = ids;
            }
            int n_shared = 0;
            for (const auto &kv : hex_all_map) {
                auto it = tet_all_map.find(kv.first);
                if (it == tet_all_map.end()) {
                    continue;
                }
                n_shared++;
                SMESH_TEST_EQ(kv.second.size(), it->second.size());
                for (size_t t = 0; t < kv.second.size(); ++t) {
                    SMESH_TEST_EQ(kv.second[t], it->second[t]);
                }
            }
            SMESH_TEST_ASSERT(n_shared > 0);
        }
    }

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

static int test_mpi_mixed_hier_l4() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 91;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_ss_mixed_hier4", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    const int       L  = 4;

    ptrdiff_t serial_nnodes = 0;
    ptrdiff_t n_coarse      = 0;
    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto mixed = create_hex8_tet4_serial(nx, ny, nz);
        SMESH_TEST_ASSERT(mixed != nullptr);
        n_coarse      = mixed->n_nodes();
        auto ss       = to_semistructured(L, mixed);
        SMESH_TEST_ASSERT(ss != nullptr);
        serial_nnodes = ss->n_nodes();
        SMESH_TEST_ASSERT(mixed->write(path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes, 1, 0);
    comm->broadcast(&n_coarse, 1, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss   = to_semistructured(L, mesh);
    auto ss_h = to_semistructured(L, mesh, true, false);
    auto ss_l2 = to_semistructured(2, mesh, true, false);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss_h != nullptr);
    SMESH_TEST_ASSERT(ss_l2 != nullptr);
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), ss->distributed()->n_nodes_global());
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, n_coarse), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, ss_l2->distributed()->n_nodes_global()), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static std::shared_ptr<Mesh> create_quad4_square_3d(const ptrdiff_t nx, const ptrdiff_t ny) {
    auto q2 = Mesh::create_quad4_square(Communicator::self(), nx, ny);
    auto p2 = q2->points()->data();
    auto p3 = create_host_buffer<geom_t>(3, q2->n_nodes());
    for (ptrdiff_t i = 0; i < q2->n_nodes(); ++i) {
        p3->data()[0][i] = p2[0][i];
        p3->data()[1][i] = p2[1][i];
        p3->data()[2][i] = 0;
    }
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    for (auto &b : q2->blocks()) {
        blocks.push_back(std::make_shared<Mesh::Block>(b->name(), b->element_type(), b->elements()));
    }
    return std::make_shared<Mesh>(Communicator::self(), blocks, p3);
}

static int test_mpi_quad4_to_ss() {
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
    const Path path = make_tmp_path("smesh_mpi_ss_quad4", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 4;
    const int       L  = 2;
    const ptrdiff_t serial_nnodes    = (static_cast<ptrdiff_t>(L) * nx + 1) * (static_cast<ptrdiff_t>(L) * ny + 1);
    const ptrdiff_t serial_nelements = nx * ny;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_quad4_square_3d(nx, ny);
        SMESH_TEST_ASSERT(serial != nullptr);
        SMESH_TEST_EQ(serial->n_elements(), serial_nelements);
        auto ss_serial = to_semistructured(L, serial);
        SMESH_TEST_ASSERT(ss_serial != nullptr);
        SMESH_TEST_EQ(ss_serial->n_nodes(), serial_nnodes);
        SMESH_TEST_ASSERT(serial->write(path) == SMESH_SUCCESS);
    }
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->is_distributed());

    auto ss = to_semistructured(L, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_ASSERT(ss->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 1);
    SMESH_TEST_EQ(ss->element_type(0), semistructured_type(QUAD4, L));
    SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), ssquad4_nxe(L));
    SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(ss->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss), SMESH_TEST_SUCCESS);

    auto ss4 = to_semistructured(4, mesh);
    SMESH_TEST_ASSERT(ss4 != nullptr);
    const ptrdiff_t serial_nnodes4 = (static_cast<ptrdiff_t>(4) * nx + 1) * (static_cast<ptrdiff_t>(4) * ny + 1);
    SMESH_TEST_EQ(ss4->distributed()->n_nodes_global(), serial_nnodes4);

    auto ss_h = to_semistructured(L, mesh, true, false);
    SMESH_TEST_ASSERT(ss_h != nullptr);
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, mesh->distributed()->n_nodes_global()), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_quad4_split_to_ss() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 43;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_ss_quad4_split", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 4;
    const int       L  = 2;
    ptrdiff_t       serial_nnodes = 0;
    const int       n_blocks      = 2;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto serial = create_quad4_square_3d(nx, ny);
        auto multi  = split_first_half(serial);
        SMESH_TEST_ASSERT(multi != nullptr);
        SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), n_blocks);
        auto ss_serial = to_semistructured(L, multi);
        SMESH_TEST_ASSERT(ss_serial != nullptr);
        serial_nnodes = ss_serial->n_nodes();
        SMESH_TEST_ASSERT(multi->write(path) == SMESH_SUCCESS);
    }
    comm->broadcast(&serial_nnodes, 1, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss = to_semistructured(L, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), n_blocks);
    SMESH_TEST_EQ(ss->element_type(0), semistructured_type(QUAD4, L));
    SMESH_TEST_EQ(ss->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss), SMESH_TEST_SUCCESS);

    auto ss_h = to_semistructured(4, mesh, true, false);
    SMESH_TEST_ASSERT(ss_h != nullptr);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, mesh->distributed()->n_nodes_global()), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss_h), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_hex8_hier_l4() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 61;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_ss_hex8_hier4", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    const ptrdiff_t ny = 2;
    const ptrdiff_t nz = 2;
    const int       L  = 4;
    const ptrdiff_t serial_nnodes = (static_cast<ptrdiff_t>(L) * nx + 1) *
                                    (static_cast<ptrdiff_t>(L) * ny + 1) *
                                    (static_cast<ptrdiff_t>(L) * nz + 1);
    const ptrdiff_t n_nodes_l2 = (2 * nx + 1) * (2 * ny + 1) * (2 * nz + 1);
    const ptrdiff_t n_coarse   = (nx + 1) * (ny + 1) * (nz + 1);
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

    auto ss = to_semistructured(L, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    auto ss_h = to_semistructured(L, mesh, true, false);
    SMESH_TEST_ASSERT(ss_h != nullptr);
    SMESH_TEST_ASSERT(ss_h->is_distributed());
    SMESH_TEST_EQ(ss_h->element_type(0), semistructured_type(HEX8, L));
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), serial_nnodes);
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), ss->distributed()->n_nodes_global());
    SMESH_TEST_EQ(ss_h->distributed()->n_elements_global(), serial_nelements);
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, n_coarse), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, n_nodes_l2), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_tet4_hier_l4() {
#ifndef SMESH_ENABLE_MPI
    return SMESH_TEST_SUCCESS;
#else
    auto comm = Communicator::world();
    if (comm->size() < 2) {
        return SMESH_TEST_SUCCESS;
    }

    int token = 0;
    if (comm->rank() == 0) {
        token = static_cast<int>(std::time(nullptr)) + 79;
    }
    comm->broadcast(&token, 1, 0);
    const Path path = make_tmp_path("smesh_mpi_ss_tet4_hier4", token);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
    ptrdiff_t       n_blocks = 0;

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
        auto single = Mesh::create_tet4_cube(Communicator::self(), nx, 2, 2);
        auto multi  = split_first_half(single);
        SMESH_TEST_ASSERT(multi != nullptr);
        n_blocks = static_cast<ptrdiff_t>(multi->n_blocks());
        SMESH_TEST_ASSERT(multi->write(path) == SMESH_SUCCESS);
    }
    comm->broadcast(&n_blocks, 1, 0);
    comm->barrier();

    auto mesh = Mesh::create_from_file(comm, path);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ss_l2 = to_semistructured(2, mesh);
    SMESH_TEST_ASSERT(ss_l2 != nullptr);
    auto ss = to_semistructured(4, mesh);
    SMESH_TEST_ASSERT(ss != nullptr);
    auto ss_h = to_semistructured(4, mesh, true, false);
    SMESH_TEST_ASSERT(ss_h != nullptr);
    SMESH_TEST_ASSERT(ss_h->is_distributed());
    SMESH_TEST_EQ(static_cast<int>(ss_h->n_blocks()), static_cast<int>(n_blocks));
    SMESH_TEST_EQ(ss_h->element_type(0), semistructured_type(TET4, 4));
    SMESH_TEST_EQ(ss_h->distributed()->n_nodes_global(), ss->distributed()->n_nodes_global());
    SMESH_TEST_EQ(check_block_layout_copied(*mesh, *ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gids_unique(*ss_h), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, mesh->distributed()->n_nodes_global()), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(check_owned_gid_prefix(*ss_h, ss_l2->distributed()->n_nodes_global()), SMESH_TEST_SUCCESS);

    if (comm->rank() == 0) {
        std::filesystem::remove_all(path.to_string());
    }
    return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_mpi_hex8_cube_to_ss);
    SMESH_RUN_TEST(test_mpi_checkerboard_to_ss);
    SMESH_RUN_TEST(test_mpi_tet4_to_ss);
    SMESH_RUN_TEST(test_mpi_quad4_to_ss);
    SMESH_RUN_TEST(test_mpi_quad4_split_to_ss);
    SMESH_RUN_TEST(test_mpi_hex8_hier_l4);
    SMESH_RUN_TEST(test_mpi_tet4_hier_l4);
    SMESH_RUN_TEST(test_mpi_mixed_hex_tet);
    SMESH_RUN_TEST(test_mpi_mixed_hier_l4);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

