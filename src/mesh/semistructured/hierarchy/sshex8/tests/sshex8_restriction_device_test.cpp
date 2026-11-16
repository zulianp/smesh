#include "smesh_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_restrict.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_test.hpp"

#ifdef SMESH_ENABLE_CUDA
#include "smesh_device_buffer.hpp"
#include "smesh_sshex8_restriction.cuh"
#endif

#include <cmath>

using namespace smesh;

#ifdef SMESH_ENABLE_CUDA

// Restriction from_level -> to_level on the device against the host, on random values. 10 macro
// elements: not a multiple of the 4 (or 3) macro elements a block of the node-centric 8 -> 4 kernel
// takes, so its last block is partial.
static int compare_restriction(const int from_level, const int to_level, const int block_size) {
    auto comm   = Communicator::self();
    auto linear = Mesh::create_hex8_cube(comm, 5, 2, 1, 0, 0, 0, 2.5, 1, 0.5);
    auto fine   = to_semistructured(from_level, linear, true, false);
    auto coarse = derefine(fine, to_level);

    auto host_restrict   = Restrict<real_t>::create(fine, coarse, EXECUTION_SPACE_HOST, block_size);
    auto device_restrict = Restrict<real_t>::create(fine, coarse, EXECUTION_SPACE_DEVICE, block_size);

    const ptrdiff_t n_from = fine->n_nodes() * block_size;
    const ptrdiff_t n_to   = coarse->n_nodes() * block_size;

    auto x      = create_host_buffer<real_t>(n_from);
    auto y_host = create_host_buffer<real_t>(n_to);
    auto y_zero = create_host_buffer<real_t>(n_to);

    // Deterministic values in [-1, 1)
    unsigned int state = 12345u;
    for (ptrdiff_t i = 0; i < n_from; i++) {
        state        = state * 1664525u + 1013904223u;
        x->data()[i] = (real_t)((state >> 8) * (1.0 / (1 << 24)) * 2 - 1);
    }

    for (ptrdiff_t i = 0; i < n_to; i++) {
        y_host->data()[i] = 0;
        y_zero->data()[i] = 0;
    }

    SMESH_TEST_ASSERT(host_restrict->apply(x->data(), y_host->data()) == SMESH_SUCCESS);

    auto d_x = to_device(x);
    auto d_y = to_device(y_zero);
    SMESH_TEST_ASSERT(device_restrict->apply(d_x->data(), d_y->data()) == SMESH_SUCCESS);
    SMESH_CUDA_CHECK(cudaDeviceSynchronize());
    SMESH_TEST_ASSERT(cudaPeekAtLastError() == cudaSuccess);
    auto y_device = to_host(d_y);

    const double tol  = sizeof(real_t) == 4 ? 1e-5 : 1e-12;
    double       norm = 0;
    for (ptrdiff_t i = 0; i < n_to; i++) {
        norm = fmax(norm, fabs((double)y_host->data()[i]));
    }

    SMESH_TEST_ASSERT(norm > 0);
    for (ptrdiff_t i = 0; i < n_to; i++) {
        const double diff = fabs((double)y_device->data()[i] - (double)y_host->data()[i]);
        if (diff > tol * norm) {
            fprintf(stderr,
                    "L%d -> L%d, block size %d: entry %ld device %g host %g\n",
                    from_level,
                    to_level,
                    block_size,
                    (long)i,
                    (double)y_device->data()[i],
                    (double)y_host->data()[i]);
            return SMESH_TEST_FAILURE;
        }
    }

    return SMESH_TEST_SUCCESS;
}

// Level 8 -> 4 with unit strides takes the node-centric kernel
static int test_restriction_8_to_4_node_centric() {
    SMESH_TEST_ASSERT(compare_restriction(8, 4, 1) == SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(compare_restriction(8, 4, 3) == SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

// Other levels take the tile-of-8 kernel
static int test_restriction_tile_of_8() {
    SMESH_TEST_ASSERT(compare_restriction(8, 2, 1) == SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(compare_restriction(4, 2, 3) == SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

// A rank without elements launches nothing and leaves no error behind
static int test_restriction_8_to_4_no_elements() {
    SMESH_TEST_ASSERT(cu_sshex8_restrict(0,
                                         8,
                                         1,
                                         nullptr,
                                         nullptr,
                                         4,
                                         1,
                                         nullptr,
                                         1,
                                         SMESH_DEFAULT,
                                         1,
                                         nullptr,
                                         SMESH_DEFAULT,
                                         1,
                                         nullptr,
                                         SMESH_DEFAULT_STREAM) == SMESH_SUCCESS);
    SMESH_CUDA_CHECK(cudaDeviceSynchronize());
    SMESH_TEST_ASSERT(cudaPeekAtLastError() == cudaSuccess);
    return SMESH_TEST_SUCCESS;
}

#else

static int test_restriction_device() { return SMESH_TEST_SKIPPED; }

#endif  // SMESH_ENABLE_CUDA

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);

#ifdef SMESH_ENABLE_CUDA
    SMESH_RUN_TEST(test_restriction_8_to_4_node_centric);
    SMESH_RUN_TEST(test_restriction_tile_of_8);
    SMESH_RUN_TEST(test_restriction_8_to_4_no_elements);
#else
    SMESH_RUN_TEST(test_restriction_device);
#endif

    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
