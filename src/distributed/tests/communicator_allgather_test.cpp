// Communicator::allgather / allgatherv.
//
// These exist because the two-step idiom they name -- allgather the per-rank counts, build
// displacements, allgatherv the payload -- was written out longhand at four separate sites
// (smesh_mesh.cpp twice, smesh_sideset.cpp, and the transforms test), each with its own
// spelling of the same MPI calls.
//
// SIZE 1 IS A REAL CASE HERE, which is why this test does not take the early return the other
// tests in this directory do. broadcast, sum and max are in-place, so a build without MPI can
// leave the buffer alone and be correct; a gather cannot. Its receive buffer is a different
// buffer, and a fallback that did nothing would hand back uninitialised memory. The one-rank
// path is therefore the path most likely to be wrong, and it is checked first.
#include <numeric>
#include <vector>

#include "smesh_communicator.hpp"
#include "smesh_context.hpp"
#include "smesh_test.hpp"

using namespace smesh;

// Fixed count: rank r contributes r*100 + i, so every entry names the rank that sent it.
static int test_allgather_fixed_count() {
    auto      comm = Communicator::world();
    const int rank = comm->rank();
    const int size = comm->size();

    const int        count = 3;
    std::vector<i64> send((size_t)count);
    for (int i = 0; i < count; ++i) send[(size_t)i] = (i64)rank * 100 + i;

    std::vector<i64> recv((size_t)count * (size_t)size, -1);
    comm->allgather(send.data(), recv.data(), count);

    // Rank-major: rank r's block starts at r * count.
    for (int r = 0; r < size; ++r)
        for (int i = 0; i < count; ++i)
            SMESH_TEST_EQ((i64)r * 100 + i, recv[(size_t)r * count + i]);

    return SMESH_TEST_SUCCESS;
}

// Variable count: rank r contributes r + 1 entries, the shape a coarse-matrix gather has,
// where no two ranks hold the same number of rows.
static int test_allgatherv_ragged() {
    auto      comm = Communicator::world();
    const int rank = comm->rank();
    const int size = comm->size();

    const int        sendcount = rank + 1;
    std::vector<f64> send((size_t)sendcount);
    for (int i = 0; i < sendcount; ++i) send[(size_t)i] = (f64)rank + (f64)i / 1000.0;

    // The counts come from the fixed-count form, which is how a caller that does not know
    // the other ranks' sizes builds the displacements.
    std::vector<int> counts((size_t)size, 0);
    comm->allgather(&sendcount, counts.data(), 1);

    std::vector<int> displs((size_t)size, 0);
    int              total = 0;
    for (int r = 0; r < size; ++r) {
        displs[(size_t)r] = total;
        total += counts[(size_t)r];
    }
    SMESH_TEST_EQ(size * (size + 1) / 2, total);

    std::vector<f64> recv((size_t)total, -1.0);
    comm->allgatherv(send.data(), sendcount, recv.data(), counts.data(), displs.data());

    for (int r = 0; r < size; ++r) {
        SMESH_TEST_EQ(r + 1, counts[(size_t)r]);
        for (int i = 0; i < counts[(size_t)r]; ++i) {
            const f64 expected = (f64)r + (f64)i / 1000.0;
            SMESH_TEST_APPROXEQ(expected, recv[(size_t)displs[(size_t)r] + i], 1e-12);
        }
    }

    return SMESH_TEST_SUCCESS;
}

// Every rank must come away with the same bytes; a gather that silently kept a rank's own
// contribution would pass the checks above on that rank alone.
static int test_allgather_agrees_across_ranks() {
    auto      comm = Communicator::world();
    const int size = comm->size();

    const int        count = 2;
    std::vector<i32> send((size_t)count);
    for (int i = 0; i < count; ++i) send[(size_t)i] = comm->rank() * 10 + i;

    std::vector<i32> recv((size_t)count * (size_t)size, 0);
    comm->allgather(send.data(), recv.data(), count);

    // Sum the gathered vector element-wise across ranks: if every rank holds the same
    // vector, each entry comes back multiplied by exactly size.
    std::vector<i32> check(recv);
    comm->sum(check.data(), (int)check.size(), TypeToEnum<i32>::value());
    for (size_t k = 0; k < recv.size(); ++k) SMESH_TEST_EQ(recv[k] * size, check[k]);

    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_allgather_fixed_count);
    SMESH_RUN_TEST(test_allgatherv_ragged);
    SMESH_RUN_TEST(test_allgather_agrees_across_ranks);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
