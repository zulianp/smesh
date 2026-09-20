// rank_split / rank_start / rank_owner must describe ONE partition of [0, n).
//
// These three are the only definition of how a global index space is cut across ranks, and
// every distributed path here depends on them agreeing: rank_start says where a rank's
// block begins, rank_split how long it is, and rank_owner which rank a given global index
// falls in. rank_owner is the inverse of the other two, and it was not: it divided the
// index by the uniform block size, which ignores that the first `n % comm_size` ranks hold
// one extra entry each, and so returned a rank past the end of the communicator.
//
// The consequence was not a wrong answer but heap corruption, because callers use the
// result as an array index with no bound of their own -- `send_displs[rank_owner(...) +
// 1]++` on a SMESH_CALLOC(comm_size + 1) buffer. A release build has every SMESH_ASSERT
// compiled out, so the first thing to notice was glibc, reporting "munmap_chunk(): invalid
// pointer" from an unrelated free much later.
//
// This test is arithmetic only. It runs at any communicator size, including one, and the
// rank counts it sweeps are independent of how many ranks it is launched with -- the point
// is to cover the (n, comm_size) pairs the distributed code actually forms, not to
// decompose anything. The pairs below are drawn from real runs: 512 and 729 are macro
// element and node counts of a cavity at N=8, 2601 and 4913 are semi-structured node counts
// at level 2, 208 is a coarse cube, and 125 and 64 are coarse multigrid levels. The rank
// counts go to 288, which is one per core on a GH200 node.

#include <vector>

#include "smesh_decompose.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static const ptrdiff_t kSizes[] = {64,  125, 208,  289,  300,  512,  577,
                                   729, 1000, 2601, 3179, 4913, 19652};
static const int kRanks[] = {1, 2, 3, 4, 8, 16, 32, 64, 128, 192, 256, 288};

// Every global index must be owned by a rank that exists.
//
// This is the property whose violation corrupted the heap: an out-of-range return is used
// directly as an array subscript by every caller.
static int test_rank_owner_is_in_range() {
    for (const ptrdiff_t n : kSizes) {
        for (const int p : kRanks) {
            for (ptrdiff_t g = 0; g < n; ++g) {
                const int owner = rank_owner(n, g, p);
                SMESH_TEST_ASSERT(owner >= 0);
                SMESH_TEST_ASSERT(owner < p);
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

// rank_owner must land inside the block rank_start and rank_split describe.
//
// Being in range is not enough: the owner has to be the rank that actually holds the index,
// or the alltoallv that follows sends data to a rank which is not expecting it.
static int test_rank_owner_inverts_rank_start() {
    for (const ptrdiff_t n : kSizes) {
        for (const int p : kRanks) {
            for (ptrdiff_t g = 0; g < n; ++g) {
                const int owner = rank_owner(n, g, p);
                const ptrdiff_t begin = rank_start(n, p, owner);
                const ptrdiff_t count = rank_split(n, p, owner);
                SMESH_TEST_ASSERT(g >= begin);
                SMESH_TEST_ASSERT(g < begin + count);
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

// The blocks must tile [0, n) exactly: no gaps, no overlaps, nothing left over.
static int test_ranges_tile_the_index_space() {
    for (const ptrdiff_t n : kSizes) {
        for (const int p : kRanks) {
            ptrdiff_t total = 0;
            for (int r = 0; r < p; ++r) {
                const ptrdiff_t begin = rank_start(n, p, r);
                const ptrdiff_t count = rank_split(n, p, r);
                SMESH_TEST_ASSERT(count >= 0);
                SMESH_TEST_EQ(begin, total);
                total += count;
            }
            SMESH_TEST_EQ(n, total);
        }
    }
    return SMESH_TEST_SUCCESS;
}

// Fewer indices than ranks is a real case, not a misuse.
//
// A coarse multigrid level has far fewer nodes and elements than a full node has ranks, and
// those counts reach these functions unchanged. The uniform block size is then zero, which
// is what the arithmetic used to divide by; the assertion that stood in for handling it is
// compiled out of every release build.
static int test_fewer_indices_than_ranks() {
    const ptrdiff_t small[] = {1, 2, 5, 17, 63, 64, 125};
    for (const ptrdiff_t n : small) {
        for (const int p : kRanks) {
            if (n >= p) continue;

            for (ptrdiff_t g = 0; g < n; ++g) {
                const int owner = rank_owner(n, g, p);
                SMESH_TEST_ASSERT(owner >= 0);
                SMESH_TEST_ASSERT(owner < p);
                SMESH_TEST_ASSERT(g >= rank_start(n, p, owner));
                SMESH_TEST_ASSERT(g < rank_start(n, p, owner) + rank_split(n, p, owner));
            }

            // The first n ranks hold one index each, the remaining ranks hold none.
            ptrdiff_t total = 0;
            for (int r = 0; r < p; ++r) {
                const ptrdiff_t count = rank_split(n, p, r);
                SMESH_TEST_EQ((ptrdiff_t)(r < n ? 1 : 0), count);
                total += count;
            }
            SMESH_TEST_EQ(n, total);
        }
    }
    return SMESH_TEST_SUCCESS;
}

// The pairs that used to fail, named so a regression is legible rather than a bare count.
//
// Each of these returned a rank equal to comm_size for at least one index: n=512 at 288
// ranks did so for 223 of its 512 indices, n=729 at 256 ranks for 215 of them. The failure
// is not monotonic in comm_size -- 4913 is wrong at 128, 192 and 256 but right at 288 --
// which is why it reads as a problem-size effect when it is arithmetic.
static int test_known_regressions() {
    struct Pair {
        ptrdiff_t n;
        int p;
    };
    const Pair pairs[] = {{512, 192}, {512, 288},  {729, 128},  {729, 192},
                          {729, 256}, {729, 288},  {2601, 64},  {2601, 128},
                          {2601, 192}, {3179, 128}, {4913, 128}, {4913, 192},
                          {4913, 256}, {1000, 256}, {208, 128},  {289, 192}};

    for (const Pair &pair : pairs) {
        std::vector<ptrdiff_t> seen((size_t)pair.p, 0);
        for (ptrdiff_t g = 0; g < pair.n; ++g) {
            const int owner = rank_owner(pair.n, g, pair.p);
            SMESH_TEST_ASSERT(owner >= 0);
            SMESH_TEST_ASSERT(owner < pair.p);
            seen[(size_t)owner]++;
        }
        for (int r = 0; r < pair.p; ++r) {
            SMESH_TEST_EQ(rank_split(pair.n, pair.p, r), seen[(size_t)r]);
        }
    }
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_rank_owner_is_in_range);
    SMESH_RUN_TEST(test_rank_owner_inverts_rank_start);
    SMESH_RUN_TEST(test_ranges_tile_the_index_space);
    SMESH_RUN_TEST(test_fewer_indices_than_ranks);
    SMESH_RUN_TEST(test_known_regressions);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
