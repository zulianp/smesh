#include <algorithm>
#include <array>
#include <cmath>
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
#include "smesh_sideset.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sspyramid.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sswedge.hpp"
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

using NodeKey = std::array<long long, 3>;
using FaceKey = std::vector<NodeKey>;

static NodeKey quantize_node(const geom_t *const *pts, const idx_t node) {
    return {std::llround(static_cast<double>(pts[0][node]) * 1e9),
            std::llround(static_cast<double>(pts[1][node]) * 1e9),
            std::llround(static_cast<double>(pts[2][node]) * 1e9)};
}

static std::set<FaceKey> collect_surface_faces(const std::shared_ptr<Mesh>                 &mesh,
                                               const std::vector<std::shared_ptr<Sideset>> &sidesets,
                                               enum ElemType                               *type_out) {
    std::set<FaceKey> faces;
    if (type_out) {
        *type_out = INVALID;
    }
    auto pts = mesh->points()->data();
    for (const auto &ss : sidesets) {
        if (!ss || ss->size() == 0) {
            continue;
        }
        auto [st, surface] = create_surface_from_sideset(mesh, ss);
        if (type_out && *type_out == INVALID) {
            *type_out = st;
        }
        if (!surface) {
            continue;
        }
        const int nnxs = static_cast<int>(surface->extent(0));
        for (ptrdiff_t e = 0; e < static_cast<ptrdiff_t>(surface->extent(1)); ++e) {
            FaceKey key(static_cast<size_t>(nnxs));
            for (int ln = 0; ln < nnxs; ++ln) {
                key[static_cast<size_t>(ln)] = quantize_node(pts, surface->data()[ln][e]);
            }
            std::sort(key.begin(), key.end());
            faces.insert(std::move(key));
        }
    }
    return faces;
}

static int test_hex8_ss_sideset_surface_multiblock() {
    auto comm   = Communicator::self();
    auto multi  = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_ASSERT(single != nullptr);

    auto multi_ss  = Sideset::create_from_plane(multi, 1, 0, 0, 0.0);
    auto single_ss = Sideset::create_from_plane(single, 1, 0, 0, 0.0);
    SMESH_TEST_ASSERT(!multi_ss.empty());
    SMESH_TEST_ASSERT(!single_ss.empty());

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);

    enum ElemType multi_st  = INVALID;
    enum ElemType single_st = INVALID;
    auto          multi_faces  = collect_surface_faces(ss_multi, multi_ss, &multi_st);
    auto          single_faces = collect_surface_faces(ss_single, single_ss, &single_st);
    SMESH_TEST_EQ(multi_st, PROTEUS_QUADSHELL9);
    SMESH_TEST_EQ(single_st, PROTEUS_QUADSHELL9);
    SMESH_TEST_EQ(elem_num_nodes(multi_st), 9);
    SMESH_TEST_ASSERT(!multi_faces.empty());
    SMESH_TEST_EQ(multi_faces.size(), single_faces.size());
    SMESH_TEST_ASSERT(multi_faces == single_faces);
    SMESH_TEST_EQ(static_cast<int>(multi_faces.begin()->size()), 9);

    ptrdiff_t n_surf_multi = 0;
    for (const auto &ss : multi_ss) {
        n_surf_multi += ss->size();
        auto [st, surface] = create_surface_from_sideset(ss_multi, ss);
        SMESH_TEST_ASSERT(surface != nullptr);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(surface->extent(0)), static_cast<ptrdiff_t>(9));
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(surface->extent(1)), ss->size());
        auto nodeset = create_nodeset_from_sideset(ss_multi, ss);
        SMESH_TEST_ASSERT(nodeset != nullptr);
    }
    ptrdiff_t n_surf_single = 0;
    for (const auto &ss : single_ss) {
        n_surf_single += ss->size();
    }
    SMESH_TEST_EQ(n_surf_multi, n_surf_single);
    return SMESH_TEST_SUCCESS;
}

static int test_tet4_ss_sideset_surface_multiblock() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);
    SMESH_TEST_ASSERT(single != nullptr);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

    auto single_ss = Sideset::create_from_plane(single, 1, 0, 0, 0.0);
    auto multi_ss  = Sideset::create_from_plane(multi, 1, 0, 0, 0.0);
    SMESH_TEST_ASSERT(!single_ss.empty());
    SMESH_TEST_ASSERT(!multi_ss.empty());

    {
        auto ss_single = to_semistructured(1, single);
        SMESH_TEST_ASSERT(ss_single != nullptr);
        for (const auto &ss : single_ss) {
            if (ss->size() == 0) {
                continue;
            }
            auto [ust, usurf] = create_surface_from_sideset(single, ss);
            auto [sst, ssurf] = create_surface_from_sideset(ss_single, ss);
            SMESH_TEST_ASSERT(usurf != nullptr);
            SMESH_TEST_ASSERT(ssurf != nullptr);
            SMESH_TEST_EQ(sst, TRISHELL3);
            SMESH_TEST_EQ(static_cast<ptrdiff_t>(ssurf->extent(0)), static_cast<ptrdiff_t>(3));
            SMESH_TEST_EQ(static_cast<ptrdiff_t>(ssurf->extent(1)), ss->size());
            SMESH_TEST_EQ(static_cast<ptrdiff_t>(usurf->extent(1)), ss->size());
            for (ptrdiff_t e = 0; e < ss->size(); ++e) {
                for (int ln = 0; ln < 3; ++ln) {
                    SMESH_TEST_EQ(ssurf->data()[ln][e], usurf->data()[ln][e]);
                }
            }
        }
    }

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);

    enum ElemType multi_st  = INVALID;
    enum ElemType single_st = INVALID;
    auto          multi_faces  = collect_surface_faces(ss_multi, multi_ss, &multi_st);
    auto          single_faces = collect_surface_faces(ss_single, single_ss, &single_st);
    SMESH_TEST_EQ(multi_st, TRISHELL6);
    SMESH_TEST_EQ(single_st, TRISHELL6);
    SMESH_TEST_EQ(elem_num_nodes(multi_st), 6);
    SMESH_TEST_ASSERT(!multi_faces.empty());
    SMESH_TEST_EQ(multi_faces.size(), single_faces.size());
    SMESH_TEST_ASSERT(multi_faces == single_faces);
    SMESH_TEST_EQ(static_cast<int>(multi_faces.begin()->size()), 6);
    return SMESH_TEST_SUCCESS;
}

static int test_mixed_hex_tet_ss_sideset_surfaces_separate() {
    auto mixed = create_hex8_tet4_serial(2, 2, 2);
    SMESH_TEST_ASSERT(mixed != nullptr);

    auto hex_ss = Sideset::create_from_plane(mixed, 0, 0, 1, 0.0, 1e-6, {"hex"});
    auto tet_ss = Sideset::create_from_plane(mixed, 0, 0, 1, 1.0, 1e-6, {"tet"});
    if (hex_ss.empty() || hex_ss[0]->size() == 0) {
        hex_ss = Sideset::create_from_plane(mixed, 1, 0, 0, 0.0, 1e-6, {"hex"});
    }
    if (tet_ss.empty() || tet_ss[0]->size() == 0) {
        tet_ss = Sideset::create_from_plane(mixed, 1, 0, 0, 1.0, 1e-6, {"tet"});
    }
    SMESH_TEST_ASSERT(!hex_ss.empty());
    SMESH_TEST_ASSERT(!tet_ss.empty());
    SMESH_TEST_ASSERT(hex_ss[0]->size() > 0);
    SMESH_TEST_ASSERT(tet_ss[0]->size() > 0);
    SMESH_TEST_EQ(static_cast<int>(hex_ss[0]->block_id()), 0);
    SMESH_TEST_EQ(static_cast<int>(tet_ss[0]->block_id()), 1);

    auto ss = to_semistructured(2, mixed);
    SMESH_TEST_ASSERT(ss != nullptr);
    SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 2);

    auto [hex_st, hex_surf] = create_surface_from_sideset(ss, hex_ss[0]);
    auto [tet_st, tet_surf] = create_surface_from_sideset(ss, tet_ss[0]);
    SMESH_TEST_ASSERT(hex_surf != nullptr);
    SMESH_TEST_ASSERT(tet_surf != nullptr);
    SMESH_TEST_EQ(hex_st, PROTEUS_QUADSHELL9);
    SMESH_TEST_EQ(tet_st, TRISHELL6);
    SMESH_TEST_ASSERT(hex_st != tet_st);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(hex_surf->extent(0)), static_cast<ptrdiff_t>(9));
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(tet_surf->extent(0)), static_cast<ptrdiff_t>(6));
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(hex_surf->extent(1)), hex_ss[0]->size());
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(tet_surf->extent(1)), tet_ss[0]->size());
    return SMESH_TEST_SUCCESS;
}

static std::shared_ptr<Mesh> create_wedge6_serial() {
    auto tri = Mesh::create_tri3_square(Communicator::self(), 2, 2);
    if (tri == nullptr) {
        return nullptr;
    }
    auto wedge = extrude(tri, 1.0, 1);
    if (wedge == nullptr || wedge->element_type(0) != WEDGE6) {
        return nullptr;
    }
    return wedge;
}

static std::shared_ptr<Mesh> create_two_pyramids_serial() {
    auto points = create_host_buffer<geom_t>(3, 6);
    const geom_t coords[6][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}, {0.5, 0.5, -1}};
    for (int i = 0; i < 6; ++i) {
        points->data()[0][i] = coords[i][0];
        points->data()[1][i] = coords[i][1];
        points->data()[2][i] = coords[i][2];
    }
    auto elems = create_host_buffer<idx_t>(5, 2);
    const idx_t p0[5] = {0, 1, 2, 3, 4};
    const idx_t p1[5] = {0, 3, 2, 1, 5};
    for (int d = 0; d < 5; ++d) {
        elems->data()[d][0] = p0[d];
        elems->data()[d][1] = p1[d];
    }
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(std::make_shared<Mesh::Block>("pyramid", PYRAMID5, elems));
    return std::make_shared<Mesh>(Communicator::self(), blocks, points);
}

static std::set<idx_t> collect_ids(const idx_t *const *ss, const ptrdiff_t e, const std::vector<int> &lidx) {
    std::set<idx_t> ids;
    for (int id : lidx) {
        ids.insert(ss[id][e]);
    }
    return ids;
}

static std::vector<int> hex_quad_interior_lidx(const int L, const int face) {
    std::vector<int> ids;
    if (L < 2) {
        return ids;
    }
    for (int a = 1; a < L; ++a) {
        for (int b = 1; b < L; ++b) {
            if (face == 0) {
                ids.push_back(sshex8_lidx(L, a, 0, b));
            } else {
                ids.push_back(sshex8_lidx(L, L, a, b));
            }
        }
    }
    return ids;
}

static std::vector<int> wedge_quad0_interior_lidx(const int L) {
    std::vector<int> ids;
    for (int x = 1; x < L; ++x) {
        for (int z = 1; z < L; ++z) {
            ids.push_back(sswedge_lidx(L, x, 0, z));
        }
    }
    return ids;
}

static std::vector<int> pyr_base_interior_lidx(const int L) {
    std::vector<int> ids;
    for (int i = 1; i < L; ++i) {
        for (int j = 1; j < L; ++j) {
            ids.push_back(sspyramid_lidx(L, i, j, 0));
        }
    }
    return ids;
}

static std::vector<int> pyr_side0_interior_lidx(const int L) {
    std::vector<int> ids;
    for (int k = 1; k <= L - 2; ++k) {
        for (int i = 1; i <= L - 1 - k; ++i) {
            ids.push_back(sspyramid_lidx(L, i, 0, k));
        }
    }
    return ids;
}

static std::vector<int> tet_face012_interior_lidx(const int L) {
    std::vector<int> ids;
    for (int y = 1; y <= L - 2; ++y) {
        for (int x = 1; x <= L - 1 - y; ++x) {
            ids.push_back(sstet4_lidx(L, x, y, 0));
        }
    }
    return ids;
}

static int test_wedge_to_semistructured() {
    auto wedge = create_wedge6_serial();
    SMESH_TEST_ASSERT(to_semistructured(2, wedge, false, true) == nullptr);
    const ptrdiff_t n_coarse = wedge->n_nodes();
    for (int L : {1, 2, 4}) {
        auto ss = to_semistructured(L, wedge);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_EQ(ss->element_type(0), semistructured_type(WEDGE6, L));
        SMESH_TEST_EQ(ss->n_elements(0), wedge->n_elements());
        SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), sswedge_nxe(L));
        if (L == 1) {
            SMESH_TEST_EQ(ss->n_nodes(), n_coarse);
        }
        auto ss_h = to_semistructured(L, wedge, true, false);
        SMESH_TEST_ASSERT(ss_h != nullptr);
        SMESH_TEST_EQ(ss_h->n_nodes(), ss->n_nodes());
        auto d = derefine(ss, 1);
        SMESH_TEST_ASSERT(d != nullptr);
        SMESH_TEST_EQ(d->n_nodes(), n_coarse);
        SMESH_TEST_EQ(d->element_type(0), semistructured_type(WEDGE6, 1));
    }
    auto one = Mesh::create_tri3_square(Communicator::self(), 1, 1);
    auto one_w = extrude(one, 1.0, 1);
    SMESH_TEST_ASSERT(one_w != nullptr);
    SMESH_TEST_EQ(one_w->n_elements(), static_cast<ptrdiff_t>(2));
    auto ss1 = to_semistructured(2, one_w);
    SMESH_TEST_ASSERT(ss1 != nullptr);
    SMESH_TEST_ASSERT(ss1->n_nodes() < static_cast<ptrdiff_t>(2) * sswedge_nxe(2));
    SMESH_TEST_ASSERT(ss1->n_nodes() > sswedge_nxe(2));
    return SMESH_TEST_SUCCESS;
}

static int test_pyramid_to_semistructured() {
    auto pyr = create_two_pyramids_serial();
    SMESH_TEST_ASSERT(to_semistructured(2, pyr, false, true) == nullptr);
    const ptrdiff_t n_coarse = pyr->n_nodes();
    for (int L : {1, 2, 4}) {
        auto ss = to_semistructured(L, pyr);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_EQ(ss->element_type(0), semistructured_type(PYRAMID5, L));
        SMESH_TEST_EQ(ss->n_elements(0), static_cast<ptrdiff_t>(2));
        SMESH_TEST_EQ(static_cast<int>(ss->block(0)->n_nodes_per_element()), sspyramid_nxe(L));
        const ptrdiff_t expected = n_coarse + static_cast<ptrdiff_t>(12) * sspyramid_nxedge(L) +
                                   static_cast<ptrdiff_t>(8) * sspyramid_nx_tri_face(L) +
                                   static_cast<ptrdiff_t>(1) * sspyramid_nx_quad_face(L) +
                                   static_cast<ptrdiff_t>(2) * sspyramid_nxvol(L);
        SMESH_TEST_EQ(ss->n_nodes(), expected);
        auto ss_h = to_semistructured(L, pyr, true, false);
        SMESH_TEST_ASSERT(ss_h != nullptr);
        SMESH_TEST_EQ(ss_h->n_nodes(), ss->n_nodes());
        auto d = derefine(ss, 1);
        SMESH_TEST_ASSERT(d != nullptr);
        SMESH_TEST_EQ(d->n_nodes(), n_coarse);
        SMESH_TEST_EQ(d->element_type(0), semistructured_type(PYRAMID5, 1));
        if (L >= 2) {
            auto pyr_ss = ss->elements(0)->data();
            const auto a = collect_ids(pyr_ss, 0, pyr_base_interior_lidx(L));
            const auto b = collect_ids(pyr_ss, 1, pyr_base_interior_lidx(L));
            SMESH_TEST_EQ(a.size(), b.size());
            SMESH_TEST_ASSERT(a == b);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_hex_dominant_to_semistructured() {
    auto hexdom = Mesh::create_hex_dominant_serial(Communicator::self());
    SMESH_TEST_ASSERT(hexdom != nullptr);
    SMESH_TEST_EQ(static_cast<int>(hexdom->n_blocks()), 4);
    SMESH_TEST_ASSERT(to_semistructured(2, hexdom, false, true) == nullptr);

    auto mixed = create_hex8_tet4_serial(2, 2, 2);
    auto ss_ht = to_semistructured(2, mixed);
    SMESH_TEST_ASSERT(ss_ht != nullptr);
    const ptrdiff_t hex_tet_n_nodes = ss_ht->n_nodes();

    for (int L : {1, 2, 4}) {
        auto ss = to_semistructured(L, hexdom);
        SMESH_TEST_ASSERT(ss != nullptr);
        SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 4);
        SMESH_TEST_EQ(ss->element_type(0), semistructured_type(HEX8, L));
        SMESH_TEST_EQ(ss->element_type(1), semistructured_type(PYRAMID5, L));
        SMESH_TEST_EQ(ss->element_type(2), semistructured_type(TET4, L));
        SMESH_TEST_EQ(ss->element_type(3), semistructured_type(WEDGE6, L));
        if (L == 1) {
            SMESH_TEST_EQ(ss->n_nodes(), hexdom->n_nodes());
            continue;
        }
        auto hex_ss = ss->elements(0)->data();
        auto pyr_ss = ss->elements(1)->data();
        auto tet_ss = ss->elements(2)->data();
        auto wed_ss = ss->elements(3)->data();
        const auto hex_pyr = collect_ids(hex_ss, 0, hex_quad_interior_lidx(L, 1));
        const auto pyr_base = collect_ids(pyr_ss, 0, pyr_base_interior_lidx(L));
        SMESH_TEST_EQ(hex_pyr.size(), static_cast<size_t>((L - 1) * (L - 1)));
        SMESH_TEST_ASSERT(hex_pyr == pyr_base);
        const auto hex_wed = collect_ids(hex_ss, 0, hex_quad_interior_lidx(L, 0));
        const auto wed_q = collect_ids(wed_ss, 0, wedge_quad0_interior_lidx(L));
        SMESH_TEST_ASSERT(hex_wed == wed_q);
        if (L >= 3) {
            const auto pyr_tri = collect_ids(pyr_ss, 0, pyr_side0_interior_lidx(L));
            const auto tet_tri = collect_ids(tet_ss, 0, tet_face012_interior_lidx(L));
            SMESH_TEST_EQ(pyr_tri.size(), tet_tri.size());
            SMESH_TEST_ASSERT(pyr_tri == tet_tri);
        }
        auto ss_h = to_semistructured(L, hexdom, true, false);
        SMESH_TEST_ASSERT(ss_h != nullptr);
        SMESH_TEST_EQ(ss_h->n_nodes(), ss->n_nodes());
        auto d = derefine(ss, 1);
        SMESH_TEST_ASSERT(d != nullptr);
        SMESH_TEST_EQ(d->n_nodes(), hexdom->n_nodes());
    }
    SMESH_TEST_EQ(to_semistructured(2, mixed)->n_nodes(), hex_tet_n_nodes);
    return SMESH_TEST_SUCCESS;
}

static int test_wedge_pyramid_ss_sideset() {
    {
        auto wedge = create_wedge6_serial();
        const ptrdiff_t n = wedge->n_elements();
        auto parents = create_host_buffer<element_idx_t>(static_cast<size_t>(n));
        auto lfi     = create_host_buffer<i16>(static_cast<size_t>(n));
        for (ptrdiff_t e = 0; e < n; ++e) {
            parents->data()[e] = static_cast<element_idx_t>(e);
            lfi->data()[e]     = 0;
        }
        auto sides = Sideset::create(wedge->comm(), parents, lfi, 0);
        auto ss    = to_semistructured(2, wedge);
        SMESH_TEST_ASSERT(ss != nullptr);
        auto [st, surf] = create_surface_from_sideset(ss, sides);
        SMESH_TEST_ASSERT(surf != nullptr);
        SMESH_TEST_EQ(st, PROTEUS_QUADSHELL9);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(surf->extent(0)), static_cast<ptrdiff_t>(9));
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(surf->extent(1)), n);
        auto nodeset = create_nodeset_from_sideset(ss, sides);
        SMESH_TEST_ASSERT(nodeset != nullptr);
        for (ptrdiff_t e = 0; e < n; ++e) {
            lfi->data()[e] = 3;
        }
        auto tri_sides = Sideset::create(wedge->comm(), parents, lfi, 0);
        auto [tst, tsurf] = create_surface_from_sideset(ss, tri_sides);
        SMESH_TEST_ASSERT(tsurf != nullptr);
        SMESH_TEST_EQ(tst, TRISHELL6);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(tsurf->extent(0)), static_cast<ptrdiff_t>(6));
    }
    {
        auto pyr = create_two_pyramids_serial();
        auto parents = create_host_buffer<element_idx_t>(2);
        auto lfi     = create_host_buffer<i16>(2);
        parents->data()[0] = 0;
        parents->data()[1] = 1;
        lfi->data()[0]     = 4;
        lfi->data()[1]     = 4;
        auto sides = Sideset::create(pyr->comm(), parents, lfi, 0);
        auto ss    = to_semistructured(2, pyr);
        SMESH_TEST_ASSERT(ss != nullptr);
        auto [st, surf] = create_surface_from_sideset(ss, sides);
        SMESH_TEST_ASSERT(surf != nullptr);
        SMESH_TEST_EQ(st, PROTEUS_QUADSHELL9);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(surf->extent(0)), static_cast<ptrdiff_t>(9));
        lfi->data()[0] = 0;
        lfi->data()[1] = 0;
        auto tri_sides = Sideset::create(pyr->comm(), parents, lfi, 0);
        auto [tst, tsurf] = create_surface_from_sideset(ss, tri_sides);
        SMESH_TEST_ASSERT(tsurf != nullptr);
        SMESH_TEST_EQ(tst, TRISHELL6);
    }
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
    SMESH_RUN_TEST(test_hex8_ss_sideset_surface_multiblock);
    SMESH_RUN_TEST(test_tet4_ss_sideset_surface_multiblock);
    SMESH_RUN_TEST(test_mixed_hex_tet_ss_sideset_surfaces_separate);
    SMESH_RUN_TEST(test_wedge_to_semistructured);
    SMESH_RUN_TEST(test_pyramid_to_semistructured);
    SMESH_RUN_TEST(test_hex_dominant_to_semistructured);
    SMESH_RUN_TEST(test_wedge_pyramid_ss_sideset);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
