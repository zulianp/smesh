#include <vector>

#include "smesh_buffer.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sstet4_graph.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static int test_sstet4_lidx_bijective() {
    for (int L = 1; L <= 5; ++L) {
        const int nxe = sstet4_nxe(L);
        SMESH_TEST_EQ(nxe, sstet4_n_tet(L));
        std::vector<unsigned char> seen(static_cast<size_t>(nxe), 0);
        int                        count = 0;
        for (int z = 0; z <= L; ++z) {
            for (int y = 0; y <= L - z; ++y) {
                for (int x = 0; x <= L - z - y; ++x) {
                    const int id = sstet4_lidx(L, x, y, z);
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
    SMESH_TEST_EQ(sstet4_nxe(1), 4);
    SMESH_TEST_EQ(sstet4_nxe(2), 10);
    SMESH_TEST_EQ(sstet4_nxe(3), 20);
    SMESH_TEST_EQ(sstet4_nxe(4), 35);
    SMESH_TEST_EQ(sstet4_txe(2), 8);
    SMESH_TEST_EQ(sstet4_txe(3), 27);
    return SMESH_TEST_SUCCESS;
}

static int test_sstet4_two_tets_shared_face() {
    const ptrdiff_t nelements = 2;
    const ptrdiff_t nnodes    = 5;
    auto            elements  = create_host_buffer<idx_t>(4, nelements);

    const idx_t t0[4] = {0, 1, 2, 3};
    const idx_t t1[4] = {1, 2, 3, 4};
    for (int d = 0; d < 4; ++d) {
        elements->data()[d][0] = t0[d];
        elements->data()[d][1] = t1[d];
    }

    for (int L : {1, 2, 3, 4}) {
        const int nxe     = sstet4_nxe(L);
        auto      ss_elems = create_host_buffer<idx_t>(nxe, nelements);
        ptrdiff_t n_unique = -1;
        ptrdiff_t interior = -1;
        SMESH_TEST_ASSERT(sstet4_generate_elements(L,
                                                   nelements,
                                                   nnodes,
                                                   elements->data(),
                                                   ss_elems->data(),
                                                   &n_unique,
                                                   &interior) == SMESH_SUCCESS);

        const int n_edges = 9;
        const int n_faces = 7;
        const ptrdiff_t expected = nnodes + static_cast<ptrdiff_t>(n_edges) * sstet4_nxedge(L) +
                                   static_cast<ptrdiff_t>(n_faces) * sstet4_nxface(L) +
                                   nelements * sstet4_nxvol(L);
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

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_sstet4_lidx_bijective);
    SMESH_RUN_TEST(test_sstet4_two_tets_shared_face);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

