
#include <stdio.h>

#include "smesh_test.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_buffer.hpp"
#include "smesh_sshex8.hpp"

using namespace smesh;

int test_sshex8_hierarchical_renumbering() {
    const ptrdiff_t nelements = 2;
    const ptrdiff_t nnodes    = 27 + 18;

    auto elements = create_host_buffer<idx_t>(27, nelements);

    for (ptrdiff_t i = 0; i < 27; i++) {
        elements->data()[i][0] = i;
        elements->data()[i][1] = 18 + i;
    }

    int       L               = 24;
    const int nxe             = sshex8_nxe(L);
    auto      sshex8_elements = create_host_buffer<idx_t>(nxe, nelements);

    ptrdiff_t sshex_nnodes   = -1;
    ptrdiff_t interior_start = -1;

    SMESH_TEST_ASSERT(
            sshex8_generate_elements(
                    L, elements->extent(1), nnodes, elements->data(), sshex8_elements->data(), &sshex_nnodes, &interior_start) ==
            SMESH_SUCCESS);

    int  nlevels = sshex8_hierarchical_n_levels(L);
    auto levels  = create_host_buffer<int>(nlevels);
    sshex8_hierarchical_mesh_levels(L, nlevels, levels->data());

    SMESH_TEST_ASSERT(nlevels == 4);
    SMESH_TEST_ASSERT(levels->data()[0] == 1);
    SMESH_TEST_ASSERT(levels->data()[1] == 6);
    SMESH_TEST_ASSERT(levels->data()[2] == 12);
    SMESH_TEST_ASSERT(levels->data()[3] == 24);

    auto node_mapping = create_host_buffer<idx_t>(sshex_nnodes);
    SMESH_TEST_ASSERT(sshex8_hierarchical_renumbering(
                             L, nlevels, levels->data(), nelements, sshex_nnodes, sshex8_elements->data(), node_mapping->data()) == SMESH_SUCCESS);

    // Check that original nodes are in range
    for (int zi = 0; zi <= 1; zi++) {
        for (int yi = 0; yi <= 1; yi++) {
            for (int xi = 0; xi <= 1; xi++) {
                int v = sshex8_lidx(L, xi * L, yi * L, zi * L);

                for (ptrdiff_t e = 0; e < nelements; e++) {
                    SMESH_TEST_ASSERT(sshex8_elements->data()[v][e] < nnodes);
                }
            }
        }
    }

    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_sshex8_hierarchical_renumbering);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
