// Collective .smz IO (smesh_smz_mpi.hpp), run at several rank counts (ctest *_np1.._np4):
//  * a file written by P ranks is byte-identical to the serial file (lossless and lossy);
//  * uniform reads give exactly matrix.io's split, any range can be read on any rank;
//  * mapped writes accept any order and duplicates; bound contributions combine by minimum;
//  * a damaged file fails on every rank (same status everywhere, no hang).

#include "smesh_byte_codec.hpp"
#include "smesh_codec_hash.hpp"
#include "smesh_smz_format.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_smz_mpi.hpp"
#include "smesh_test.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <vector>

using namespace smesh;
using namespace smesh::compression;
namespace fs = std::filesystem;

// One folder per rank count: ctest -j runs the _np1.._np4 variants concurrently.
static std::string kDir;

static double value_at(i64 g) { return std::sin(0.001 * double(g)) * 50 + 1e-3 * double(g % 13); }

static std::vector<u8> slurp(const std::string &p) {
    std::ifstream in(p, std::ios::binary);
    return std::vector<u8>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

static int rank_of(MPI_Comm c) {
    int r;
    MPI_Comm_rank(c, &r);
    return r;
}
static int size_of(MPI_Comm c) {
    int s;
    MPI_Comm_size(c, &s);
    return s;
}

// Contiguous blocks of uneven size, written in parallel == written serially.
static int test_parallel_equals_serial() {
    MPI_Comm   comm = MPI_COMM_WORLD;
    const int  rank = rank_of(comm), size = size_of(comm);
    const i64  n    = 123457;
    for (const char *spec : {"lossless,chunk:1000", "range:1e-5,chunk:4096", "pwrel:1e-3,zstd,chunk:777"}) {
        CompressionOptions o;
        SMESH_TEST_EQ(CompressionOptions::parse(spec, &o), SMESH_SUCCESS);
        if (!byte_stage_available(o.byte_stage)) continue;
        // Uneven blocks: rank r holds (r + 1) shares.
        const i64 shares = i64(size) * (size + 1) / 2;
        i64       begin  = n * (i64(rank) * (rank + 1) / 2) / shares;
        i64       end    = n * (i64(rank + 1) * (rank + 2) / 2) / shares;
        std::vector<f64> x(size_t(end - begin));
        for (i64 g = begin; g < end; g++) x[size_t(g - begin)] = value_at(g);
        const std::string par = kDir + "/par.float64.smz", ser = kDir + "/ser.float64.smz";
        SMESH_TEST_EQ(smz_write_contiguous<f64>(comm, Path(par), x.data(), ptrdiff_t(x.size()), ptrdiff_t(n), o),
                      SMESH_SUCCESS);
        if (rank == 0) {
            std::vector<f64> all(static_cast<size_t>(n));
            for (i64 g = 0; g < n; g++) all[size_t(g)] = value_at(g);
            SMESH_TEST_EQ(smz_write<f64>(Path(ser), all.data(), ptrdiff_t(n), o), SMESH_SUCCESS);
            const auto a = slurp(par), b = slurp(ser);
            if (a != b) {
                fprintf(stderr, "%s: parallel file (%zu B) differs from serial (%zu B)\n", spec, a.size(), b.size());
                return SMESH_TEST_FAILURE;
            }
        }
        MPI_Barrier(comm);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_uniform_read() {
    MPI_Comm  comm = MPI_COMM_WORLD;
    const int rank = rank_of(comm), size = size_of(comm);
    const i64 n    = 100003;
    const std::string p = kDir + "/u.float64.smz";
    if (rank == 0) {
        std::vector<f64> all(static_cast<size_t>(n));
        for (i64 g = 0; g < n; g++) all[size_t(g)] = value_at(g);
        CompressionOptions o;
        SMESH_TEST_EQ(CompressionOptions::parse("lossless,chunk:999", &o), SMESH_SUCCESS);
        SMESH_TEST_EQ(smz_write<f64>(Path(p), all.data(), ptrdiff_t(n), o), SMESH_SUCCESS);
    }
    MPI_Barrier(comm);
    f32      *y = nullptr;  // read with conversion
    ptrdiff_t nl = 0, ng = 0;
    SMESH_TEST_EQ(smz_read_uniform<f32>(comm, Path(p), &y, &nl, &ng), SMESH_SUCCESS);
    i64 b, c;
    uniform_split(n, size, rank, &b, &c);
    SMESH_TEST_EQ(i64(nl), c);
    SMESH_TEST_EQ(i64(ng), n);
    for (i64 i = 0; i < c; i++) SMESH_TEST_ASSERT(y[i] == f32(value_at(b + i)));
    free(y);
    // Arbitrary ranges, including empty ones.
    const i64        lo = (rank * 7919) % n, hi = std::min(n, lo + 5000 * (rank % 3));
    std::vector<f64> z(size_t(hi - lo));
    SMESH_TEST_EQ(smz_read_range<f64>(comm, Path(p), ptrdiff_t(lo), ptrdiff_t(hi), z.data()), SMESH_SUCCESS);
    for (i64 i = lo; i < hi; i++) SMESH_TEST_ASSERT(z[size_t(i - lo)] == value_at(i));
    return SMESH_TEST_SUCCESS;
}

// Each rank sends a scattered, overlapping set of (index, value) pairs; bounds from several ranks
// per value combine by minimum.
static int test_mapped_and_min_bounds() {
    MPI_Comm  comm = MPI_COMM_WORLD;
    const int rank = rank_of(comm), size = size_of(comm);
    const i64 n    = 50000;
    std::vector<i64> map;
    std::vector<f32> val;
    std::vector<double> bnd;
    for (i64 g = n - 1 - rank; g >= 0; g -= size) {  // owned, in reverse order
        map.push_back(g);
        val.push_back(f32(value_at(g)));
    }
    for (i64 g = (rank * 31) % size; g < n; g += 3) {  // duplicates owned by others
        map.push_back(g);
        val.push_back(f32(value_at(g)));
    }
    // Bound contributions: every rank contributes for every 2nd value; the true bound is the min
    // over ranks, 1e-4 * (1 + g % 5), contributed by the last rank.
    std::vector<i64> cmap;
    for (i64 g = 0; g < n; g += 2) {
        cmap.push_back(g);
        bnd.push_back(1e-4 * (1 + g % 5) * (rank == size - 1 ? 1.0 : 3.0 + rank));
    }
    for (i64 g = 1; g < n; g += 2) {  // odd values: only rank 0 contributes
        if (rank == 0) {
            cmap.push_back(g);
            bnd.push_back(1e-4 * (1 + g % 5));
        }
    }
    CompressionOptions o;
    SMESH_TEST_EQ(CompressionOptions::parse("h:1e-2,chunk:3000,verify", &o), SMESH_SUCCESS);
    const std::string p = kDir + "/m.float32.smz";
    const int wret = smz_write_mapped<f32, i64>(comm, Path(p), ptrdiff_t(map.size()), ptrdiff_t(n), map.data(),
                                                val.data(), o, ptrdiff_t(cmap.size()), cmap.data(), bnd.data());
    SMESH_TEST_EQ(wret, SMESH_SUCCESS);
    if (rank == 0) {
        f32      *y  = nullptr;
        ptrdiff_t ny = 0;
        SMESH_TEST_EQ(smz_read_convert<f32>(Path(p), &y, &ny), SMESH_SUCCESS);
        SMESH_TEST_EQ(i64(ny), n);
        for (i64 g = 0; g < n; g++) {
            const double err = std::fabs(double(y[g]) - double(f32(value_at(g))));
            SMESH_TEST_ASSERT(err <= 1e-4 * (1 + g % 5));
        }
        free(y);
    }
    MPI_Barrier(comm);
    return SMESH_TEST_SUCCESS;
}

static int test_damaged_file_fails_everywhere() {
    MPI_Comm          comm = MPI_COMM_WORLD;
    const int         rank = rank_of(comm);
    const std::string p    = kDir + "/bad.float64.smz";
    if (rank == 0) {
        std::vector<f64> all(20000);
        for (size_t g = 0; g < all.size(); g++) all[g] = value_at(i64(g));
        CompressionOptions o;
        SMESH_TEST_EQ(CompressionOptions::parse("lossless,chunk:1000", &o), SMESH_SUCCESS);
        SMESH_TEST_EQ(smz_write<f64>(Path(p), all.data(), 20000, o), SMESH_SUCCESS);
        std::vector<u8> b = slurp(p);
        b[b.size() - 10] ^= 0x40;  // the last chunk: only the last rank reads it
        std::ofstream(p, std::ios::binary).write(reinterpret_cast<const char *>(b.data()), std::streamsize(b.size()));
    }
    MPI_Barrier(comm);
    f64      *y = nullptr;
    ptrdiff_t nl = 0, ng = 0;
    SMESH_TEST_EQ(smz_read_uniform<f64>(comm, Path(p), &y, &nl, &ng), SMESH_FAILURE);
    SMESH_TEST_ASSERT(y == nullptr);
    SMESH_TEST_EQ(smz_read_uniform<f64>(comm, Path(kDir + "/missing.float64.smz"), &y, &nl, &ng), SMESH_FAILURE);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    const int rank = rank_of(MPI_COMM_WORLD);
    kDir           = "smz_mpi_test_tmp_np" + std::to_string(size_of(MPI_COMM_WORLD));
    if (rank == 0) {
        fs::remove_all(kDir);
        fs::create_directories(kDir);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    SMESH_RUN_TEST(test_parallel_equals_serial);
    SMESH_RUN_TEST(test_uniform_read);
    SMESH_RUN_TEST(test_mapped_and_min_bounds);
    SMESH_RUN_TEST(test_damaged_file_fails_everywhere);
    int any = err, all = 0;
    MPI_Allreduce(&any, &all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0 && !all) fs::remove_all(kDir);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
