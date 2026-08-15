#include <algorithm>
#include <cstring>
#include <cstdint>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

#include "smesh_conversion.hpp"
#include "smesh_mesh.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static std::shared_ptr<Mesh> create_hex8_tet4_serial(const ptrdiff_t nx,
                                                     const ptrdiff_t ny,
                                                     const ptrdiff_t nz) {
    auto            cube      = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
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

static std::shared_ptr<Mesh> split_first_half(const std::shared_ptr<Mesh> &mesh) {
    auto              out     = mesh->clone();
    const ptrdiff_t   n       = out->n_elements(0);
    const ptrdiff_t   n_split = n / 2;
    auto              parents = create_host_buffer<element_idx_t>(static_cast<size_t>(n_split));
    for (ptrdiff_t i = 0; i < n_split; ++i) {
        parents->data()[i] = static_cast<element_idx_t>(i);
    }
    if (out->split_block(parents, "part0", 0) != SMESH_SUCCESS) {
        return nullptr;
    }
    return out;
}

static int node_map_is_permutation(const SharedBuffer<idx_t> &node_map, const ptrdiff_t n_nodes) {
    SMESH_TEST_ASSERT(node_map != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(node_map->size()), n_nodes);
    std::vector<unsigned char> seen(static_cast<size_t>(n_nodes), 0);
    auto                       map = node_map->data();
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const idx_t dst = map[i];
        SMESH_TEST_ASSERT(dst >= 0);
        SMESH_TEST_ASSERT(dst < n_nodes);
        SMESH_TEST_ASSERT(seen[static_cast<size_t>(dst)] == 0);
        seen[static_cast<size_t>(dst)] = 1;
    }
    return SMESH_TEST_SUCCESS;
}

static int test_checkerboard_to_semistructured() {
    auto comm  = Communicator::self();
    auto multi = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_ASSERT(single != nullptr);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_ASSERT(ss_multi->block(0)->name() == "white");
    SMESH_TEST_ASSERT(ss_multi->block(1)->name() == "black");
    SMESH_TEST_EQ(ss_multi->element_type(0), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(ss_multi->element_type(1), semistructured_type(HEX8, 2));

    auto hex_multi = sshex_to_hex8(ss_multi);
    SMESH_TEST_ASSERT(hex_multi != nullptr);
    SMESH_TEST_EQ(static_cast<int>(hex_multi->n_blocks()), 2);
    const ptrdiff_t txe = sshex8_txe(2);
    SMESH_TEST_EQ(hex_multi->n_elements(0), multi->n_elements(0) * txe);
    SMESH_TEST_EQ(hex_multi->n_elements(1), multi->n_elements(1) * txe);
    SMESH_TEST_EQ(hex_multi->n_elements(), ss_single->n_elements() * txe);

    return SMESH_TEST_SUCCESS;
}

static int test_checkerboard_to_semistructured_hierarchical() {
    auto comm  = Communicator::self();
    auto multi = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);

    auto ss_multi  = to_semistructured(2, multi, true, false);
    auto ss_single = to_semistructured(2, single, true, false);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    return SMESH_TEST_SUCCESS;
}

static int test_checkerboard_derefine() {
    auto comm  = Communicator::self();
    auto multi = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    auto d_multi   = derefine(ss_multi, 1);
    auto d_single  = derefine(ss_single, 1);
    SMESH_TEST_ASSERT(d_multi != nullptr);
    SMESH_TEST_ASSERT(d_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(d_multi->n_blocks()), 2);
    SMESH_TEST_EQ(d_multi->n_nodes(), d_single->n_nodes());
    SMESH_TEST_EQ(d_multi->element_type(0), semistructured_type(HEX8, 1));
    return SMESH_TEST_SUCCESS;
}

static int hex_edge_ss_ids(const int                    L,
                           const idx_t *const *         coarse,
                           const idx_t *const *         ss,
                           const ptrdiff_t              n_e,
                           const idx_t                  lo,
                           const idx_t                  hi,
                           std::vector<idx_t>          &ids) {
    static const int hex_edges[12][2] = {{0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3},
                                         {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};
    ids.clear();
    for (ptrdiff_t e = 0; e < n_e; ++e) {
        for (int k = 0; k < 12; ++k) {
            const idx_t a = coarse[hex_edges[k][0]][e];
            const idx_t b = coarse[hex_edges[k][1]][e];
            const idx_t u = std::min(a, b);
            const idx_t v = std::max(a, b);
            if (u != lo || v != hi) {
                continue;
            }
            const int c0 = (a == lo) ? hex_edges[k][0] : hex_edges[k][1];
            const int c1 = (a == lo) ? hex_edges[k][1] : hex_edges[k][0];
            static const int xyz[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                          {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
            const int P1[3] = {xyz[c0][0] * L, xyz[c0][1] * L, xyz[c0][2] * L};
            const int P2[3] = {xyz[c1][0] * L, xyz[c1][1] * L, xyz[c1][2] * L};
            ids.resize(static_cast<size_t>(L > 1 ? (L - 1) : 0));
            for (int t = 1; t < L; ++t) {
                const int xi = (P1[0] * (L - t) + P2[0] * t) / L;
                const int yi = (P1[1] * (L - t) + P2[1] * t) / L;
                const int zi = (P1[2] * (L - t) + P2[2] * t) / L;
                ids[static_cast<size_t>(t - 1)] = ss[sshex8_lidx(L, xi, yi, zi)][e];
            }
            return SMESH_TEST_SUCCESS;
        }
    }
    return SMESH_TEST_FAILURE;
}

static int tet_edge_ss_ids(const int                    L,
                           const idx_t *const *         coarse,
                           const idx_t *const *         ss,
                           const ptrdiff_t              n_e,
                           const idx_t                  lo,
                           const idx_t                  hi,
                           std::vector<idx_t>          &ids) {
    ids.clear();
    for (ptrdiff_t e = 0; e < n_e; ++e) {
        int c0 = -1;
        int c1 = -1;
        for (int d = 0; d < 4; ++d) {
            if (coarse[d][e] == lo) {
                c0 = d;
            }
            if (coarse[d][e] == hi) {
                c1 = d;
            }
        }
        if (c0 < 0 || c1 < 0) {
            continue;
        }
        static const int xyz[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const int P1[3] = {xyz[c0][0] * L, xyz[c0][1] * L, xyz[c0][2] * L};
        const int P2[3] = {xyz[c1][0] * L, xyz[c1][1] * L, xyz[c1][2] * L};
        ids.resize(static_cast<size_t>(L > 1 ? (L - 1) : 0));
        for (int t = 1; t < L; ++t) {
            const int xi = (P1[0] * (L - t) + P2[0] * t) / L;
            const int yi = (P1[1] * (L - t) + P2[1] * t) / L;
            const int zi = (P1[2] * (L - t) + P2[2] * t) / L;
            ids[static_cast<size_t>(t - 1)] = ss[sstet4_lidx(L, xi, yi, zi)][e];
        }
        return SMESH_TEST_SUCCESS;
    }
    return SMESH_TEST_FAILURE;
}

static void collect_hex_edges(const std::shared_ptr<Mesh> &mesh, std::set<std::pair<idx_t, idx_t>> &edges) {
    static const int hex_edges[12][2] = {{0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3},
                                         {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};
    auto             coarse           = mesh->elements(0)->data();
    const ptrdiff_t  n_e              = mesh->n_elements(0);
    for (ptrdiff_t e = 0; e < n_e; ++e) {
        for (int k = 0; k < 12; ++k) {
            const idx_t a = coarse[hex_edges[k][0]][e];
            const idx_t b = coarse[hex_edges[k][1]][e];
            edges.emplace(std::min(a, b), std::max(a, b));
        }
    }
}

static void collect_tet_edges(const std::shared_ptr<Mesh> &mesh, std::set<std::pair<idx_t, idx_t>> &edges) {
    auto            coarse = mesh->elements(1)->data();
    const ptrdiff_t n_e    = mesh->n_elements(1);
    for (ptrdiff_t e = 0; e < n_e; ++e) {
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                const idx_t a = coarse[i][e];
                const idx_t b = coarse[j][e];
                edges.emplace(std::min(a, b), std::max(a, b));
            }
        }
    }
}

static int test_mixed_to_semistructured() {
    auto mixed = create_hex8_tet4_serial(2, 2, 2);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_EQ(static_cast<int>(mixed->n_blocks()), 2);
    SMESH_TEST_EQ(mixed->element_type(0), HEX8);
    SMESH_TEST_EQ(mixed->element_type(1), TET4);

    SMESH_TEST_ASSERT(to_semistructured(2, mixed, false, true) == nullptr);

    std::set<std::pair<idx_t, idx_t>> hex_edges;
    std::set<std::pair<idx_t, idx_t>> tet_edges;
    collect_hex_edges(mixed, hex_edges);
    collect_tet_edges(mixed, tet_edges);
    std::set<std::pair<idx_t, idx_t>> shared;
    for (const auto &e : hex_edges) {
        if (tet_edges.count(e)) {
            shared.insert(e);
        }
    }
    SMESH_TEST_ASSERT(!shared.empty());

    std::vector<std::shared_ptr<Mesh::Block>> hex_blocks = {mixed->block(0)};
    std::vector<std::shared_ptr<Mesh::Block>> tet_blocks = {mixed->block(1)};
    auto hex_only = std::make_shared<Mesh>(mixed->comm(), hex_blocks, mixed->points());
    auto tet_only = std::make_shared<Mesh>(mixed->comm(), tet_blocks, mixed->points());

    const ptrdiff_t n_corners = mixed->n_nodes();
    auto            hex_c     = mixed->elements(0)->data();
    auto            tet_c     = mixed->elements(1)->data();

    for (int L : {1, 2, 4}) {
        auto ss = to_semistructured(L, mixed);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 2);
        SMESH_TEST_ASSERT(ss->block(0)->name() == "hex");
        SMESH_TEST_ASSERT(ss->block(1)->name() == "tet");
        SMESH_TEST_EQ(ss->element_type(0), semistructured_type(HEX8, L));
        SMESH_TEST_EQ(ss->element_type(1), semistructured_type(TET4, L));
        SMESH_TEST_EQ(ss->n_elements(0), mixed->n_elements(0));
        SMESH_TEST_EQ(ss->n_elements(1), mixed->n_elements(1));
        SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), sshex8_nxe(L));
        SMESH_TEST_EQ(static_cast<int>(ss->block(1)->n_nodes_per_element()), sstet4_nxe(L));

        auto hex_ss = ss->elements(0)->data();
        auto tet_ss = ss->elements(1)->data();
        const int hex_cidx[8] = {sshex8_lidx(L, 0, 0, 0),
                                 sshex8_lidx(L, L, 0, 0),
                                 sshex8_lidx(L, L, L, 0),
                                 sshex8_lidx(L, 0, L, 0),
                                 sshex8_lidx(L, 0, 0, L),
                                 sshex8_lidx(L, L, 0, L),
                                 sshex8_lidx(L, L, L, L),
                                 sshex8_lidx(L, 0, L, L)};
        const int tet_cidx[4] = {sstet4_lidx(L, 0, 0, 0),
                                 sstet4_lidx(L, L, 0, 0),
                                 sstet4_lidx(L, 0, L, 0),
                                 sstet4_lidx(L, 0, 0, L)};
        for (ptrdiff_t e = 0; e < mixed->n_elements(0); ++e) {
            for (int d = 0; d < 8; ++d) {
                SMESH_TEST_EQ(hex_ss[hex_cidx[d]][e], hex_c[d][e]);
            }
        }
        for (ptrdiff_t e = 0; e < mixed->n_elements(1); ++e) {
            for (int d = 0; d < 4; ++d) {
                SMESH_TEST_EQ(tet_ss[tet_cidx[d]][e], tet_c[d][e]);
            }
        }

        if (L == 1) {
            SMESH_TEST_EQ(ss->n_nodes(), n_corners);
            continue;
        }

        for (const auto &edge : shared) {
            std::vector<idx_t> hex_ids;
            std::vector<idx_t> tet_ids;
            SMESH_TEST_EQ(hex_edge_ss_ids(L, hex_c, hex_ss, mixed->n_elements(0), edge.first, edge.second, hex_ids),
                          SMESH_TEST_SUCCESS);
            SMESH_TEST_EQ(tet_edge_ss_ids(L, tet_c, tet_ss, mixed->n_elements(1), edge.first, edge.second, tet_ids),
                          SMESH_TEST_SUCCESS);
            SMESH_TEST_EQ(hex_ids.size(), tet_ids.size());
            for (size_t i = 0; i < hex_ids.size(); ++i) {
                SMESH_TEST_EQ(hex_ids[i], tet_ids[i]);
            }
        }

        auto ss_hex = to_semistructured(L, hex_only);
        auto ss_tet = to_semistructured(L, tet_only);
        SMESH_TEST_ASSERT(ss_hex != nullptr);
        SMESH_TEST_ASSERT(ss_tet != nullptr);
        const ptrdiff_t expected =
                ss_hex->n_nodes() + ss_tet->n_nodes() - n_corners -
                static_cast<ptrdiff_t>(shared.size()) * static_cast<ptrdiff_t>(L - 1);
        SMESH_TEST_EQ(ss->n_nodes(), expected);
    }

    auto ss2 = to_semistructured(2, mixed);
    SMESH_TEST_ASSERT(ss2 != nullptr);
    auto ss2h = to_semistructured(2, mixed, true, false);
    SMESH_TEST_ASSERT(ss2h != nullptr);
    SMESH_TEST_EQ(ss2h->n_blocks(), ss2->n_blocks());
    SMESH_TEST_EQ(ss2h->n_nodes(), ss2->n_nodes());
    SMESH_TEST_EQ(ss2h->element_type(0), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(ss2h->element_type(1), semistructured_type(TET4, 2));
    {
        auto hex_ss = ss2h->elements(0)->data();
        auto tet_ss = ss2h->elements(1)->data();
        auto hex_c  = mixed->elements(0)->data();
        auto tet_c  = mixed->elements(1)->data();
        const int hex_cidx[8] = {sshex8_lidx(2, 0, 0, 0),
                                 sshex8_lidx(2, 2, 0, 0),
                                 sshex8_lidx(2, 2, 2, 0),
                                 sshex8_lidx(2, 0, 2, 0),
                                 sshex8_lidx(2, 0, 0, 2),
                                 sshex8_lidx(2, 2, 0, 2),
                                 sshex8_lidx(2, 2, 2, 2),
                                 sshex8_lidx(2, 0, 2, 2)};
        const int tet_cidx[4] = {sstet4_lidx(2, 0, 0, 0),
                                 sstet4_lidx(2, 2, 0, 0),
                                 sstet4_lidx(2, 0, 2, 0),
                                 sstet4_lidx(2, 0, 0, 2)};
        for (ptrdiff_t e = 0; e < mixed->n_elements(0); ++e) {
            for (int d = 0; d < 8; ++d) {
                SMESH_TEST_EQ(hex_ss[hex_cidx[d]][e], hex_c[d][e]);
            }
        }
        for (ptrdiff_t e = 0; e < mixed->n_elements(1); ++e) {
            for (int d = 0; d < 4; ++d) {
                SMESH_TEST_EQ(tet_ss[tet_cidx[d]][e], tet_c[d][e]);
            }
        }
        for (const auto &edge : shared) {
            std::vector<idx_t> hex_ids;
            std::vector<idx_t> tet_ids;
            SMESH_TEST_EQ(hex_edge_ss_ids(2, hex_c, hex_ss, mixed->n_elements(0), edge.first, edge.second, hex_ids),
                          SMESH_TEST_SUCCESS);
            SMESH_TEST_EQ(tet_edge_ss_ids(2, tet_c, tet_ss, mixed->n_elements(1), edge.first, edge.second, tet_ids),
                          SMESH_TEST_SUCCESS);
            SMESH_TEST_EQ(hex_ids.size(), tet_ids.size());
            for (size_t i = 0; i < hex_ids.size(); ++i) {
                SMESH_TEST_EQ(hex_ids[i], tet_ids[i]);
            }
        }
    }

    auto d = derefine(ss2, 1);
    SMESH_TEST_ASSERT(d != nullptr);
    SMESH_TEST_EQ(static_cast<int>(d->n_blocks()), 2);
    SMESH_TEST_ASSERT(d->block(0)->name() == "hex");
    SMESH_TEST_ASSERT(d->block(1)->name() == "tet");
    SMESH_TEST_EQ(d->element_type(0), semistructured_type(HEX8, 1));
    SMESH_TEST_EQ(d->element_type(1), semistructured_type(TET4, 1));
    SMESH_TEST_EQ(d->n_nodes(), mixed->n_nodes());
    SMESH_TEST_ASSERT(sshex_to_hex8(ss2) == nullptr);

    return SMESH_TEST_SUCCESS;
}

static int test_tet4_to_semistructured() {
    auto comm  = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);
    SMESH_TEST_ASSERT(single != nullptr);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

    for (int L : {1, 2, 3, 4}) {
        auto ss_multi  = to_semistructured(L, multi);
        auto ss_single = to_semistructured(L, single);
        SMESH_TEST_ASSERT(ss_multi != nullptr);
        SMESH_TEST_ASSERT(ss_single != nullptr);
        SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
        SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
        SMESH_TEST_EQ(ss_multi->element_type(0), semistructured_type(TET4, L));
        SMESH_TEST_EQ(ss_multi->element_type(1), semistructured_type(TET4, L));
        SMESH_TEST_EQ(static_cast<int>(ss_multi->block(0)->n_nodes_per_element()), sstet4_nxe(L));
        SMESH_TEST_EQ(ss_multi->n_elements(0), multi->n_elements(0));
        SMESH_TEST_EQ(ss_multi->n_elements(1), multi->n_elements(1));
        SMESH_TEST_EQ(ss_multi->n_elements(), single->n_elements());
    }

    return SMESH_TEST_SUCCESS;
}

static int test_tet4_to_semistructured_hierarchical() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);

    auto ss_multi  = to_semistructured(2, multi, true, false);
    auto ss_single = to_semistructured(2, single, true, false);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    return SMESH_TEST_SUCCESS;
}

static int test_tet4_derefine() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    auto d_multi   = derefine(ss_multi, 1);
    auto d_single  = derefine(ss_single, 1);
    SMESH_TEST_ASSERT(d_multi != nullptr);
    SMESH_TEST_ASSERT(d_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(d_multi->n_blocks()), 2);
    SMESH_TEST_EQ(d_multi->n_nodes(), d_single->n_nodes());
    SMESH_TEST_EQ(d_multi->element_type(0), semistructured_type(TET4, 1));
    SMESH_TEST_EQ(d_multi->n_nodes(), single->n_nodes());
    return SMESH_TEST_SUCCESS;
}

static int test_quad4_to_semistructured() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_quad4_square(comm, 4, 4);
    auto multi  = split_first_half(single);
    SMESH_TEST_ASSERT(single != nullptr);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

    for (int L : {2, 4}) {
        auto ss_multi  = to_semistructured(L, multi);
        auto ss_single = to_semistructured(L, single);
        SMESH_TEST_ASSERT(ss_multi != nullptr);
        SMESH_TEST_ASSERT(ss_single != nullptr);
        SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
        SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
        const ptrdiff_t expect_nnodes = (static_cast<ptrdiff_t>(L) * 4 + 1) * (static_cast<ptrdiff_t>(L) * 4 + 1);
        SMESH_TEST_EQ(ss_single->n_nodes(), expect_nnodes);
        SMESH_TEST_EQ(ss_multi->element_type(0), semistructured_type(QUAD4, L));
        SMESH_TEST_EQ(ss_multi->element_type(1), semistructured_type(QUAD4, L));
        SMESH_TEST_EQ(static_cast<int>(ss_multi->block(0)->n_nodes_per_element()), ssquad4_nxe(L));
        SMESH_TEST_EQ(ss_multi->n_elements(0), multi->n_elements(0));
        SMESH_TEST_EQ(ss_multi->n_elements(1), multi->n_elements(1));
        SMESH_TEST_EQ(ss_multi->n_elements(), single->n_elements());

        std::set<idx_t> ids0;
        std::set<idx_t> ids1;
        auto            e0 = ss_multi->block(0)->elements()->data();
        auto            e1 = ss_multi->block(1)->elements()->data();
        const int       nxe = ssquad4_nxe(L);
        for (int v = 0; v < nxe; ++v) {
            for (ptrdiff_t e = 0; e < ss_multi->n_elements(0); ++e) {
                ids0.insert(e0[v][e]);
            }
            for (ptrdiff_t e = 0; e < ss_multi->n_elements(1); ++e) {
                ids1.insert(e1[v][e]);
            }
        }
        std::vector<idx_t> shared;
        std::set_intersection(ids0.begin(), ids0.end(), ids1.begin(), ids1.end(), std::back_inserter(shared));
        SMESH_TEST_ASSERT(!shared.empty());
    }

    return SMESH_TEST_SUCCESS;
}

static int test_quad4_to_semistructured_hierarchical() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_quad4_square(comm, 4, 4);
    auto multi  = split_first_half(single);

    auto ss_multi  = to_semistructured(2, multi, true, false);
    auto ss_single = to_semistructured(2, single, true, false);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    SMESH_TEST_EQ(ss_multi->element_type(0), semistructured_type(QUAD4, 2));
    return SMESH_TEST_SUCCESS;
}

static int test_quad4_derefine() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_quad4_square(comm, 4, 4);
    auto multi  = split_first_half(single);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    auto d_multi   = derefine(ss_multi, 1);
    auto d_single  = derefine(ss_single, 1);
    SMESH_TEST_ASSERT(d_multi != nullptr);
    SMESH_TEST_ASSERT(d_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(d_multi->n_blocks()), 2);
    SMESH_TEST_EQ(d_multi->n_nodes(), d_single->n_nodes());
    SMESH_TEST_EQ(d_multi->element_type(0), semistructured_type(QUAD4, 1));
    SMESH_TEST_EQ(d_multi->n_nodes(), single->n_nodes());
    return SMESH_TEST_SUCCESS;
}

static int test_quad_mixed_rejected() {
    auto comm = Communicator::self();
    auto quad = Mesh::create_quad4_square(comm, 2, 2);
    SMESH_TEST_ASSERT(quad != nullptr);
    auto qel = quad->elements(0)->data();
    auto tel = create_host_buffer<idx_t>(3, 1);
    tel->data()[0][0] = qel[0][0];
    tel->data()[1][0] = qel[1][0];
    tel->data()[2][0] = qel[2][0];
    std::vector<std::shared_ptr<Mesh::Block>> qt_blocks;
    qt_blocks.push_back(std::make_shared<Mesh::Block>(quad->block(0)->name(), QUAD4, quad->elements(0)));
    qt_blocks.push_back(std::make_shared<Mesh::Block>("tri", TRI3, tel));
    auto quad_tri = std::make_shared<Mesh>(comm, qt_blocks, quad->points());
    SMESH_TEST_ASSERT(to_semistructured(2, quad_tri) == nullptr);

    auto hex = Mesh::create_hex8_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(hex != nullptr);
    auto hexel = hex->elements(0)->data();
    auto qhex  = create_host_buffer<idx_t>(4, 1);
    for (int d = 0; d < 4; ++d) {
        qhex->data()[d][0] = hexel[d][0];
    }
    std::vector<std::shared_ptr<Mesh::Block>> qh_blocks;
    qh_blocks.push_back(std::make_shared<Mesh::Block>(hex->block(0)->name(), HEX8, hex->elements(0)));
    qh_blocks.push_back(std::make_shared<Mesh::Block>("quad", QUAD4, qhex));
    auto quad_hex = std::make_shared<Mesh>(comm, qh_blocks, hex->points());
    SMESH_TEST_ASSERT(to_semistructured(2, quad_hex) == nullptr);
    return SMESH_TEST_SUCCESS;
}

static int test_quadshell4_to_semistructured() {
    auto comm = Communicator::self();
    auto ring = Mesh::create_quad4_ring(comm, 0.5, 1.0, 2, 8);
    SMESH_TEST_ASSERT(ring != nullptr);
    ring->set_element_type(0, QUADSHELL4);
    SMESH_TEST_EQ(ring->element_type(0), QUADSHELL4);
    auto ss = to_semistructured(2, ring);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_EQ(ss->element_type(0), semistructured_type(QUADSHELL4, 2));
    SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), ssquad4_nxe(2));
    SMESH_TEST_EQ(ss->spatial_dimension(), 3);
    return SMESH_TEST_SUCCESS;
}

static int test_packed_checkerboard() {
    auto mesh = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto packed = PackedMesh<uint16_t>::create(mesh, {}, false);
    SMESH_TEST_ASSERT(packed != nullptr);
    SMESH_TEST_EQ(packed->n_blocks(), static_cast<ptrdiff_t>(2));
    SMESH_TEST_ASSERT(packed->block_name(0) == "white");
    SMESH_TEST_ASSERT(packed->block_name(1) == "black");
    for (int b = 0; b < 2; ++b) {
        auto elems = packed->elements(b);
        SMESH_TEST_ASSERT(elems != nullptr);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(elems->extent(1)), mesh->n_elements(static_cast<block_idx_t>(b)));
        SMESH_TEST_EQ(static_cast<int>(elems->extent(0)), 8);
    }
    SMESH_TEST_EQ(node_map_is_permutation(packed->node_map(), mesh->n_nodes()), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_packed_hex8_tet4() {
    auto mesh = create_hex8_tet4_serial(2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto packed = PackedMesh<uint16_t>::create(mesh, {}, false);
    SMESH_TEST_ASSERT(packed != nullptr);
    SMESH_TEST_EQ(packed->n_blocks(), static_cast<ptrdiff_t>(2));
    SMESH_TEST_EQ(static_cast<int>(packed->elements(0)->extent(0)), 8);
    SMESH_TEST_EQ(static_cast<int>(packed->elements(1)->extent(0)), 4);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(packed->elements(0)->extent(1)), mesh->n_elements(0));
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(packed->elements(1)->extent(1)), mesh->n_elements(1));
    SMESH_TEST_EQ(node_map_is_permutation(packed->node_map(), mesh->n_nodes()), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_checkerboard_to_semistructured);
    SMESH_RUN_TEST(test_checkerboard_to_semistructured_hierarchical);
    SMESH_RUN_TEST(test_checkerboard_derefine);
    SMESH_RUN_TEST(test_mixed_to_semistructured);
    SMESH_RUN_TEST(test_tet4_to_semistructured);
    SMESH_RUN_TEST(test_tet4_to_semistructured_hierarchical);
    SMESH_RUN_TEST(test_tet4_derefine);
    SMESH_RUN_TEST(test_quad4_to_semistructured);
    SMESH_RUN_TEST(test_quad4_to_semistructured_hierarchical);
    SMESH_RUN_TEST(test_quad4_derefine);
    SMESH_RUN_TEST(test_quad_mixed_rejected);
    SMESH_RUN_TEST(test_quadshell4_to_semistructured);
    SMESH_RUN_TEST(test_packed_checkerboard);
    SMESH_RUN_TEST(test_packed_hex8_tet4);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

