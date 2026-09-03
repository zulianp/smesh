#include <vector>

#include "smesh_buffer.hpp"
#include "smesh_sspyramid.hpp"
#include "smesh_sspyramid_graph.hpp"
#include "smesh_sswedge.hpp"
#include "smesh_sswedge_graph.hpp"
#include "smesh_sswedge_mesh.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static int test_sswedge_lidx_bijective() {
    for (int L = 1; L <= 5; ++L) {
        const int                      nxe = sswedge_nxe(L);
        std::vector<unsigned char>     seen(static_cast<size_t>(nxe), 0);
        int                            count = 0;
        for (int z = 0; z <= L; ++z) {
            for (int y = 0; y <= L; ++y) {
                for (int x = 0; x <= L - y; ++x) {
                    const int id = sswedge_lidx(L, x, y, z);
                    SMESH_TEST_ASSERT(id >= 0);
                    SMESH_TEST_ASSERT(id < nxe);
                    SMESH_TEST_ASSERT(seen[static_cast<size_t>(id)] == 0);
                    seen[static_cast<size_t>(id)] = 1;
                    count++;
                }
            }
        }
        SMESH_TEST_EQ(count, nxe);
    }
    SMESH_TEST_EQ(sswedge_nxe(1), 6);
    SMESH_TEST_EQ(sswedge_nxe(2), 18);
    SMESH_TEST_EQ(sswedge_nxe(3), 40);
    SMESH_TEST_EQ(sswedge_nxe(4), 75);
    SMESH_TEST_EQ(sswedge_nxe(8), 405);
    return SMESH_TEST_SUCCESS;
}

static int test_sspyramid_lidx_bijective() {
    for (int L = 1; L <= 5; ++L) {
        const int                      nxe = sspyramid_nxe(L);
        std::vector<unsigned char>     seen(static_cast<size_t>(nxe), 0);
        int                            count = 0;
        for (int k = 0; k <= L; ++k) {
            for (int j = 0; j <= L - k; ++j) {
                for (int i = 0; i <= L - k; ++i) {
                    const int id = sspyramid_lidx(L, i, j, k);
                    SMESH_TEST_ASSERT(id >= 0);
                    SMESH_TEST_ASSERT(id < nxe);
                    SMESH_TEST_ASSERT(seen[static_cast<size_t>(id)] == 0);
                    seen[static_cast<size_t>(id)] = 1;
                    count++;
                }
            }
        }
        SMESH_TEST_EQ(count, nxe);
    }
    SMESH_TEST_EQ(sspyramid_nxe(1), 5);
    SMESH_TEST_EQ(sspyramid_nxe(2), 14);
    SMESH_TEST_EQ(sspyramid_nxe(3), 30);
    SMESH_TEST_EQ(sspyramid_nxe(4), 55);
    SMESH_TEST_EQ(sspyramid_nxe(8), 285);
    return SMESH_TEST_SUCCESS;
}

static int test_sswedge_two_prisms_shared_quad() {
    const ptrdiff_t nelements = 2;
    const ptrdiff_t nnodes    = 8;
    auto            elements  = create_host_buffer<idx_t>(6, nelements);
    const idx_t     w0[6]     = {0, 1, 2, 4, 5, 6};
    const idx_t     w1[6]     = {1, 3, 2, 5, 7, 6};
    for (int d = 0; d < 6; ++d) {
        elements->data()[d][0] = w0[d];
        elements->data()[d][1] = w1[d];
    }
    for (int L : {1, 2, 3, 4}) {
        const int nxe      = sswedge_nxe(L);
        auto      ss_elems = create_host_buffer<idx_t>(nxe, nelements);
        ptrdiff_t n_unique = -1;
        ptrdiff_t interior = -1;
        SMESH_TEST_ASSERT(sswedge_generate_elements(L,
                                                    nelements,
                                                    nnodes,
                                                    elements->data(),
                                                    ss_elems->data(),
                                                    &n_unique,
                                                    &interior) == SMESH_SUCCESS);
        const int n_edges = 14;
        const int n_tri   = 4;
        const int n_quad  = 5;
        const ptrdiff_t expected = nnodes + static_cast<ptrdiff_t>(n_edges) * sswedge_nxedge(L) +
                                   static_cast<ptrdiff_t>(n_tri) * sswedge_nx_tri_face(L) +
                                   static_cast<ptrdiff_t>(n_quad) * sswedge_nx_quad_face(L) +
                                   nelements * sswedge_nxvol(L);
        SMESH_TEST_EQ(n_unique, expected);
        std::vector<unsigned char> seen(static_cast<size_t>(n_unique), 0);
        for (int i = 0; i < nxe; ++i) {
            for (ptrdiff_t e = 0; e < nelements; ++e) {
                const idx_t id = ss_elems->data()[i][e];
                SMESH_TEST_ASSERT(id >= 0);
                SMESH_TEST_ASSERT(id < n_unique);
                seen[static_cast<size_t>(id)] = 1;
            }
        }
        ptrdiff_t used = 0;
        for (ptrdiff_t i = 0; i < n_unique; ++i) {
            used += seen[static_cast<size_t>(i)];
        }
        SMESH_TEST_EQ(used, n_unique);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_sspyramid_two_shared_base() {
    const ptrdiff_t nelements = 2;
    const ptrdiff_t nnodes    = 6;
    auto            elements  = create_host_buffer<idx_t>(5, nelements);
    const idx_t     p0[5]     = {0, 1, 2, 3, 4};
    const idx_t     p1[5]     = {0, 3, 2, 1, 5};
    for (int d = 0; d < 5; ++d) {
        elements->data()[d][0] = p0[d];
        elements->data()[d][1] = p1[d];
    }
    for (int L : {1, 2, 3, 4}) {
        const int nxe      = sspyramid_nxe(L);
        auto      ss_elems = create_host_buffer<idx_t>(nxe, nelements);
        ptrdiff_t n_unique = -1;
        ptrdiff_t interior = -1;
        SMESH_TEST_ASSERT(sspyramid_generate_elements(L,
                                                      nelements,
                                                      nnodes,
                                                      elements->data(),
                                                      ss_elems->data(),
                                                      &n_unique,
                                                      &interior) == SMESH_SUCCESS);
        const int n_edges = 12;
        const int n_tri   = 8;
        const int n_quad  = 1;
        const ptrdiff_t expected = nnodes + static_cast<ptrdiff_t>(n_edges) * sspyramid_nxedge(L) +
                                   static_cast<ptrdiff_t>(n_tri) * sspyramid_nx_tri_face(L) +
                                   static_cast<ptrdiff_t>(n_quad) * sspyramid_nx_quad_face(L) +
                                   nelements * sspyramid_nxvol(L);
        SMESH_TEST_EQ(n_unique, expected);
        std::vector<unsigned char> seen(static_cast<size_t>(n_unique), 0);
        for (int i = 0; i < nxe; ++i) {
            for (ptrdiff_t e = 0; e < nelements; ++e) {
                const idx_t id = ss_elems->data()[i][e];
                SMESH_TEST_ASSERT(id >= 0);
                SMESH_TEST_ASSERT(id < n_unique);
                seen[static_cast<size_t>(id)] = 1;
            }
        }
        ptrdiff_t used = 0;
        for (ptrdiff_t i = 0; i < n_unique; ++i) {
            used += seen[static_cast<size_t>(i)];
        }
        SMESH_TEST_EQ(used, n_unique);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_sswedge_to_standard_wedge6() {
    const ptrdiff_t nelements = 1;
    auto            macro     = create_host_buffer<idx_t>(6, nelements);
    for (int d = 0; d < 6; ++d) {
        macro->data()[d][0] = static_cast<idx_t>(d);
    }
    for (int L : {1, 2, 3, 4}) {
        const int nxe      = sswedge_nxe(L);
        const int txe      = sswedge_txe(L);
        auto      ss_elems = create_host_buffer<idx_t>(nxe, nelements);
        ptrdiff_t n_unique = -1;
        ptrdiff_t interior = -1;
        SMESH_TEST_ASSERT(sswedge_generate_elements(L,
                                                    nelements,
                                                    6,
                                                    macro->data(),
                                                    ss_elems->data(),
                                                    &n_unique,
                                                    &interior) == SMESH_SUCCESS);
        SMESH_TEST_EQ(txe, L * L * L);
        auto wedge = create_host_buffer<idx_t>(6, txe);
        SMESH_TEST_ASSERT(sswedge_to_standard_wedge6_mesh(L, nelements, ss_elems->data(), wedge->data()) ==
                          SMESH_SUCCESS);
        if (L == 1) {
            for (int d = 0; d < 6; ++d) {
                SMESH_TEST_EQ(wedge->data()[d][0], static_cast<idx_t>(d));
            }
        }
        for (int d = 0; d < 6; ++d) {
            for (int e = 0; e < txe; ++e) {
                SMESH_TEST_ASSERT(wedge->data()[d][e] >= 0);
                SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(wedge->data()[d][e]) < n_unique);
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_sswedge_lidx_bijective);
    SMESH_RUN_TEST(test_sspyramid_lidx_bijective);
    SMESH_RUN_TEST(test_sswedge_two_prisms_shared_quad);
    SMESH_RUN_TEST(test_sspyramid_two_shared_base);
    SMESH_RUN_TEST(test_sswedge_to_standard_wedge6);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
