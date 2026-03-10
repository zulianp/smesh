#include <memory>


#include "smesh_test.hpp"
#include "smesh_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_ssquad4_prolongation.hpp"
#include "smesh_ssquad4_restriction.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

using namespace smesh;

int test_trace_space_operations(const std::shared_ptr<Mesh>                 &coarse_mesh,
                                const std::shared_ptr<Mesh>                 &fine_mesh,
                                const std::vector<std::shared_ptr<Sideset>> &sideset,
                                const std::string                           & /*name*/,
                                const int                                    block_size,
                                const ExecutionSpace                         es) {
    if (sideset.size() > 1) {
        SMESH_ERROR("Not implemented!\n");
    }

    auto coarse_x = create_host_buffer<real_t>(coarse_mesh->n_nodes() * block_size);
    auto fine_x   = create_buffer<real_t>(fine_mesh->n_nodes() * block_size, es);

    ptrdiff_t n_nodes{0};
    idx_t    *nodes{nullptr};
    SMESH_TEST_ASSERT(smesh::sshex8_extract_nodeset_from_sideset(smesh::semistructured_level(*coarse_mesh),
                                                                 coarse_mesh->elements(0)->data(),
                                                                 sideset[0]->parent()->size(),
                                                                 sideset[0]->parent()->data(),
                                                                 sideset[0]->lfi()->data(),
                                                                 &n_nodes,
                                                                 &nodes) == SMESH_SUCCESS);

    auto coarse_nodeset = manage_host_buffer(n_nodes, nodes);

    {
        const ptrdiff_t n      = coarse_nodeset->size();
        auto            idx    = coarse_nodeset->data();
        auto            data   = coarse_x->data();
        for (ptrdiff_t i = 0; i < n; i++) {
            data[idx[i] * block_size] = 1;
        }
    }

#ifdef SMESH_ENABLE_CUDA
    if (es == EXECUTION_SPACE_DEVICE) coarse_x = to_device(coarse_x);
#endif

    const int fine_level   = smesh::semistructured_level(*fine_mesh);
    const int coarse_level = smesh::semistructured_level(*coarse_mesh);
    const int nexs         = (fine_level + 1) * (fine_level + 1);
    auto      fine_sides   = create_host_buffer<idx_t>(nexs, sideset[0]->parent()->size());

    SMESH_TEST_ASSERT(smesh::sshex8_extract_surface_from_sideset(fine_level,
                                                                 fine_mesh->elements(0)->data(),
                                                                 sideset[0]->parent()->size(),
                                                                 sideset[0]->parent()->data(),
                                                                 sideset[0]->lfi()->data(),
                                                                 fine_sides->data()) == SMESH_SUCCESS);

    SMESH_TEST_ASSERT(ssquad4_prolongate(fine_sides->extent(1),      // nelements,
                                         coarse_level,               // rom_level
                                         fine_level / coarse_level,  // from_level_stride
                                         fine_sides->data(),         // from_elements
                                         fine_level,                 // to_level
                                         1,                          // to_level_stride
                                         fine_sides->data(),         // to_elements
                                         block_size,    // vec_size
                                         coarse_x->data(),
                                         fine_x->data()) == SMESH_SUCCESS);

    auto restricted_x = create_buffer<real_t>(coarse_mesh->n_nodes() * block_size, es);
    auto count        = create_host_buffer<uint16_t>(fine_mesh->n_nodes());

#ifdef SMESH_ENABLE_CUDA
    if (es == EXECUTION_SPACE_DEVICE) count = to_device(count);
#endif

    SMESH_TEST_ASSERT(ssquad4_element_node_incidence_count(
                              fine_level, 1, fine_sides->extent(1), fine_sides->data(), count->data()) == SMESH_SUCCESS);

    SMESH_TEST_ASSERT(ssquad4_restrict(fine_sides->extent(1),
                                       fine_level,
                                       1,
                                       fine_sides->data(),
                                       count->data(),
                                       coarse_level,
                                       fine_level / coarse_level,
                                       fine_sides->data(),
                                       block_size,
                                       fine_x->data(),
                                       restricted_x->data()) == SMESH_SUCCESS);

#ifdef SMESH_ENABLE_CUDA
    if (es == EXECUTION_SPACE_DEVICE) {
        restricted_x = to_host(restricted_x);
    }
#endif

    {
        const ptrdiff_t n      = coarse_nodeset->size();

        auto rx = restricted_x->data();
        auto cx = coarse_x->data();

        for (ptrdiff_t i = 0; i < n; i++) {
            SMESH_TEST_ASSERT((cx[i] != 0) == (rx[i] != 0));
        }
    }

    return SMESH_TEST_SUCCESS;
}

int test_trace_space_prolongation_restriction() {
    MPI_Comm comm = MPI_COMM_WORLD;
    auto     es   = EXECUTION_SPACE_HOST;

    const char *SMESH_EXECUTION_SPACE{nullptr};
    SMESH_READ_ENV(SMESH_EXECUTION_SPACE, );

    if (SMESH_EXECUTION_SPACE) {
        es = execution_space_from_string(SMESH_EXECUTION_SPACE);
    }

    int SMESH_ELEMENT_REFINE_LEVEL = 4;
    SMESH_READ_ENV(SMESH_ELEMENT_REFINE_LEVEL, atoi);

    int SMESH_BASE_RESOLUTION = 4;
    SMESH_READ_ENV(SMESH_BASE_RESOLUTION, atoi);

    geom_t Lx       = 1;
    auto   m        = Mesh::create_hex8_cube(Communicator::wrap(comm),
                                             SMESH_BASE_RESOLUTION * 1,
                                             SMESH_BASE_RESOLUTION * 1,
                                             SMESH_BASE_RESOLUTION * 1,
                                             0,
                                             0,
                                             0,
                                             Lx,
                                             1,
                                             1);
    m               = smesh::to_semistructured(SMESH_ELEMENT_REFINE_LEVEL, m, true, false);
    int  block_size = 1;
    auto sideset = Sideset::create_from_selector(
            m, [=](const geom_t /*x*/, const geom_t y, const geom_t /*z*/) -> bool { return y > -1e-5 && y < 1e-5; });

    return test_trace_space_operations(derefine(m, 2), m, {sideset}, "test_trace_space_prolongation_restriction", block_size, es);
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);

    SMESH_RUN_TEST(test_trace_space_prolongation_restriction);

    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
