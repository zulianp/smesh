// .smz files: round trips for every value type and file type conversion, partial reads across
// chunk boundaries, lossy bounds through files, identical bytes for any thread count, dispatch
// through array_{read,write}_convert_from_extension and Buffer, and rejection (never a crash or
// abort) of damaged files: every header byte flipped, payload bytes flipped, truncation, an
// unknown version.

#include "smesh_array_count.hpp"
#include "smesh_buffer.hpp"
#include "smesh_buffer.impl.hpp"
#include "smesh_byte_codec.hpp"
#include "smesh_codec_hash.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_read.hpp"
#include "smesh_read.impl.hpp"
#include "smesh_smz_format.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_test.hpp"
#include "smesh_write.hpp"
#include "smesh_write.impl.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace smesh;
using namespace smesh::compression;
namespace fs = std::filesystem;

static const fs::path kDir = "smz_file_test_tmp";

static Path tmp(const std::string &name) { return Path((kDir / name).string()); }

static std::vector<u8> slurp(const Path &p) {
    std::ifstream in(p.c_str(), std::ios::binary);
    return std::vector<u8>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
static void spit(const Path &p, const std::vector<u8> &b) {
    std::ofstream out(p.c_str(), std::ios::binary);
    out.write(reinterpret_cast<const char *>(b.data()), std::streamsize(b.size()));
}

template <class T>
static std::vector<T> sample(size_t n) {
    std::vector<T> v(n);
    for (size_t i = 0; i < n; i++) {
        if constexpr (is_float_v<T>) {
            v[i] = T(std::sin(0.01 * double(i)) * 100);
        } else {
            v[i] = T((i * 37) % 1000 + i / 64);
        }
    }
    return v;
}

template <class T>
static int roundtrip_type() {
    for (size_t n : {size_t(0), size_t(1), size_t(1000), size_t(70000)}) {
        for (i64 chunk : {i64(1), i64(7), i64(65536)}) {
            if (chunk == 1 && n > 1000) continue;
            const std::vector<T> x = sample<T>(n);
            CompressionOptions   o = CompressionOptions::lossless();
            o.chunk_values         = chunk;
            o.verify_write         = true;
            const Path p           = tmp(std::string("rt.") + to_string(DTypeOf<T>::value) + ".smz");
            SMESH_TEST_EQ(smz_write<T>(p, x.data(), ptrdiff_t(n), o), SMESH_SUCCESS);
            T        *y  = nullptr;
            ptrdiff_t ny = -1;
            SMESH_TEST_EQ(smz_read_convert<T>(p, &y, &ny), SMESH_SUCCESS);
            SMESH_TEST_EQ(ny, ptrdiff_t(n));
            SMESH_TEST_ASSERT(n == 0 || memcmp(y, x.data(), n * sizeof(T)) == 0);
            free(y);
            SMESH_TEST_EQ(array_count_from_file(p, sizeof(T)), ptrdiff_t(n));
            // A range across chunk boundaries.
            if (n >= 1000) {
                std::vector<T> part(500);
                SMESH_TEST_EQ(smz_read_range<T>(p, 123, 623, part.data()), SMESH_SUCCESS);
                SMESH_TEST_ASSERT(memcmp(part.data(), x.data() + 123, 500 * sizeof(T)) == 0);
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_roundtrip_types() {
    if (roundtrip_type<i8>() || roundtrip_type<u8>() || roundtrip_type<i16>() || roundtrip_type<u16>() ||
        roundtrip_type<i32>() || roundtrip_type<u32>() || roundtrip_type<i64>() || roundtrip_type<u64>() ||
        roundtrip_type<f32>() || roundtrip_type<f64>())
        return SMESH_TEST_FAILURE;
    return SMESH_TEST_SUCCESS;
}

// The file keeps its type; reading converts.
static int test_conversions() {
    const std::vector<f32> x = sample<f32>(5000);
    const Path             p = tmp("conv.float32.smz");
    SMESH_TEST_EQ(smz_write<f32>(p, x.data(), 5000, CompressionOptions::lossless()), SMESH_SUCCESS);
    f64      *y = nullptr;
    ptrdiff_t n = 0;
    SMESH_TEST_EQ(smz_read_convert<f64>(p, &y, &n), SMESH_SUCCESS);
    for (ptrdiff_t i = 0; i < n; i++) SMESH_TEST_ASSERT(y[i] == f64(x[size_t(i)]));
    free(y);
    // Written through the extension dispatch: idx-like i64 data stored as int32.
    std::vector<i64> idx(3000);
    for (size_t i = 0; i < idx.size(); i++) idx[i] = i64(i * 3);
    const Path q = tmp("i0.int32.smz");
    SMESH_TEST_EQ(array_write_convert_from_extension<i64>(q, idx.data(), 3000), SMESH_SUCCESS);
    SmzInfo info;
    SMESH_TEST_EQ(smz_read_info(q, &info), SMESH_SUCCESS);
    SMESH_TEST_ASSERT(info.plan.dtype == DType::I32);
    i64      *back = nullptr;
    ptrdiff_t nb   = 0;
    SMESH_TEST_EQ(array_read_convert_from_extension<i64>(q, &back, &nb), SMESH_SUCCESS);
    SMESH_TEST_EQ(nb, ptrdiff_t(3000));
    SMESH_TEST_ASSERT(memcmp(back, idx.data(), 3000 * sizeof(i64)) == 0);
    free(back);
    return SMESH_TEST_SUCCESS;
}

static int test_lossy_files() {
    const size_t     n = 100000;
    std::vector<f64> x(n);
    for (size_t i = 0; i < n; i++) x[i] = std::exp(std::sin(0.0003 * double(i))) - 1.0;
    for (const char *spec : {"abs:1e-4", "range:1e-5", "pwrel:1e-3", "abs:1e-4,zstd", "range:1e-3,chunk:1000,verify"}) {
        CompressionOptions o;
        SMESH_TEST_EQ(CompressionOptions::parse(spec, &o), SMESH_SUCCESS);
        if (!byte_stage_available(o.byte_stage)) continue;
        const Path p = tmp("lossy.float64.smz");
        SmzReport  r;
        SMESH_TEST_EQ(smz_write<f64>(p, x.data(), ptrdiff_t(n), o, nullptr, &r), SMESH_SUCCESS);
        f64      *y  = nullptr;
        ptrdiff_t ny = 0;
        SMESH_TEST_EQ(smz_read_convert<f64>(p, &y, &ny), SMESH_SUCCESS);
        double lo = x[0], hi = x[0];
        for (double v : x) lo = std::min(lo, v), hi = std::max(hi, v);
        double max_err = 0;
        for (size_t i = 0; i < n; i++) {
            const double err = std::fabs(y[i] - x[i]);
            const double b   = o.mode == BoundMode::Abs        ? o.eb
                               : o.mode == BoundMode::RelRange ? o.eb * (hi - lo)
                                                               : o.eb * std::fabs(x[i]);
            SMESH_TEST_ASSERT(err <= b);
            max_err = std::max(max_err, err);
        }
        free(y);
        SmzInfo info;
        SMESH_TEST_EQ(smz_read_info(p, &info), SMESH_SUCCESS);
        SMESH_TEST_ASSERT(info.max_abs_err == max_err);
        SMESH_TEST_ASSERT(r.ratio() > 2);
        printf("%-30s %.2fx  max err/bound %.3f  escaped %llu\n", spec, r.ratio(), r.max_ratio,
               (unsigned long long)r.n_escaped);
    }
    // Lossy bounds on integers are refused (not silently ignored).
    std::vector<i32>   iv(10, 1);
    CompressionOptions o;
    SMESH_TEST_EQ(CompressionOptions::parse("abs:1", &o), SMESH_SUCCESS);
    SMESH_TEST_EQ(smz_write<i32>(tmp("bad.int32.smz"), iv.data(), 10, o), SMESH_FAILURE);
    return SMESH_TEST_SUCCESS;
}

// Too many values off the grid: chunks fall back to lossless, values are still exact.
static int test_escape_fallback() {
    const size_t     n = 4096;
    std::vector<f32> x(n);
    std::mt19937     rng(1);
    for (size_t i = 0; i < n; i++) x[i] = f32(std::ldexp(1.0 + (rng() % 1000) / 1000.0, int(rng() % 200) - 100));
    CompressionOptions o;
    SMESH_TEST_EQ(CompressionOptions::parse("abs:1e-40", &o), SMESH_SUCCESS);
    const Path p = tmp("esc.float32.smz");
    SmzReport  r;
    SMESH_TEST_EQ(smz_write<f32>(p, x.data(), ptrdiff_t(n), o, nullptr, &r), SMESH_SUCCESS);
    f32      *y  = nullptr;
    ptrdiff_t ny = 0;
    SMESH_TEST_EQ(smz_read_convert<f32>(p, &y, &ny), SMESH_SUCCESS);
    SMESH_TEST_ASSERT(memcmp(y, x.data(), n * sizeof(f32)) == 0);
    free(y);
    return SMESH_TEST_SUCCESS;
}

static int test_thread_independence() {
#ifdef _OPENMP
    const size_t     n = 300000;
    std::vector<f64> x(n);
    for (size_t i = 0; i < n; i++) x[i] = std::cos(0.0001 * double(i)) + 1e-3 * double(i % 17);
    CompressionOptions o;
    SMESH_TEST_EQ(CompressionOptions::parse("range:1e-6,chunk:4096", &o), SMESH_SUCCESS);
    const int nt = omp_get_max_threads();
    omp_set_num_threads(1);
    SMESH_TEST_EQ(smz_write<f64>(tmp("t1.float64.smz"), x.data(), ptrdiff_t(n), o), SMESH_SUCCESS);
    omp_set_num_threads(nt > 1 ? nt : 4);
    SMESH_TEST_EQ(smz_write<f64>(tmp("tn.float64.smz"), x.data(), ptrdiff_t(n), o), SMESH_SUCCESS);
    omp_set_num_threads(nt);
    SMESH_TEST_ASSERT(slurp(tmp("t1.float64.smz")) == slurp(tmp("tn.float64.smz")));
    return SMESH_TEST_SUCCESS;
#else
    return SMESH_TEST_SKIPPED;
#endif
}

static int test_corruption() {
    const std::vector<f32> x = sample<f32>(20000);
    CompressionOptions     o;
    SMESH_TEST_EQ(CompressionOptions::parse("range:1e-4,chunk:2000", &o), SMESH_SUCCESS);
    const Path good = tmp("good.float32.smz");
    SMESH_TEST_EQ(smz_write<f32>(good, x.data(), 20000, o), SMESH_SUCCESS);
    const std::vector<u8> bytes = slurp(good);
    const Path            bad   = tmp("bad.float32.smz");
    f32                  *y     = nullptr;
    ptrdiff_t             n     = 0;

    for (size_t k = 0; k < kSmzHeaderBytes; k++) {  // every header byte
        std::vector<u8> b = bytes;
        b[k] ^= 0x01;
        spit(bad, b);
        SMESH_TEST_EQ(smz_read_convert<f32>(bad, &y, &n), SMESH_FAILURE);
    }
    for (size_t k = kSmzHeaderBytes; k < bytes.size(); k += 97) {  // index and payload bytes
        std::vector<u8> b = bytes;
        b[k] ^= 0x10;
        spit(bad, b);
        SMESH_TEST_EQ(smz_read_convert<f32>(bad, &y, &n), SMESH_FAILURE);
    }
    for (size_t cut : {size_t(0), size_t(64), size_t(130), bytes.size() / 2, bytes.size() - 1}) {
        std::vector<u8> b(bytes.begin(), bytes.begin() + ptrdiff_t(cut));
        spit(bad, b);
        SMESH_TEST_EQ(smz_read_convert<f32>(bad, &y, &n), SMESH_FAILURE);
    }
    {  // a future major version, with a valid header checksum
        std::vector<u8> b = bytes;
        b[8]              = 2;
        const u64 h       = xxh64(b.data(), 120);
        memcpy(b.data() + 120, &h, 8);
        spit(bad, b);
        SmzInfo info;
        SMESH_TEST_EQ(smz_read_info(bad, &info), SMESH_FAILURE);
    }
    SMESH_TEST_ASSERT(y == nullptr);
    return SMESH_TEST_SUCCESS;
}

// Mixed folders: raw and compressed files of one stem.
static int test_detect_files() {
    const fs::path d = kDir / "mixed";
    fs::create_directories(d);
    std::vector<i32> v(10, 3);
    const Path       folder(d.string());
    SMESH_TEST_EQ(array_write_convert_from_extension<i32>(folder / "i0.int32", v.data(), 10), SMESH_SUCCESS);
    SMESH_TEST_EQ(array_write_convert_from_extension<i32>(folder / "i1.int32.smz", v.data(), 10), SMESH_SUCCESS);
    std::ofstream(folder / "i2.int32.smz.tmp").put('x');  // leftover of an interrupted write
    auto files = detect_files(folder / "i*.*", {"raw", "int16", "int32", "int64"});
    SMESH_TEST_EQ(files.size(), size_t(2));
    for (const auto &f : files) {
        const std::string stem = typed_stem(f);
        SMESH_TEST_ASSERT(stem == "i0" || stem == "i1");
        SMESH_TEST_ASSERT(typed_extension(f) == "int32");
    }
    SMESH_TEST_EQ(typed_extension(std::string_view("a.b/x.float32.smz")), std::string_view("float32"));
    SMESH_TEST_EQ(typed_extension(std::string_view("x.smz")), std::string_view(""));
    // Raw wins over compressed for the same type, unless the compressed name is preferred.
    std::vector<Path> both = {Path("f/x.float32.smz"), Path("f/x.float32")};
    SMESH_TEST_EQ(select_one_typed_file(both, "").to_string(), std::string("f/x.float32"));
    SMESH_TEST_EQ(select_one_typed_file(both, "x.float32.smz").to_string(), std::string("f/x.float32.smz"));
    return SMESH_TEST_SUCCESS;
}

static int test_buffer_and_options() {
    auto b = create_host_buffer<f64>(1000);
    for (ptrdiff_t i = 0; i < 1000; i++) b->data()[i] = double(i) * 0.5;
    CompressionOptions o;
    SMESH_TEST_EQ(CompressionOptions::parse("lossless,chunk:100", &o), SMESH_SUCCESS);
    SMESH_TEST_EQ(b->to_file(tmp("buf.float64"), o), SMESH_SUCCESS);
    auto c = Buffer<f64>::from_file(tmp("buf.float64.smz"));
    SMESH_TEST_ASSERT(c && c->size() == 1000);
    SMESH_TEST_ASSERT(memcmp(c->data(), b->data(), 1000 * sizeof(f64)) == 0);

    CompressionOptions p;
    SMESH_TEST_EQ(CompressionOptions::parse("range:1e-3,zstd:5,shuffle,chunk:4096,verify", &p), SMESH_SUCCESS);
    SMESH_TEST_ASSERT(p.enabled && p.mode == BoundMode::RelRange && p.eb == 1e-3 && p.zstd_level == 5 && p.shuffle &&
                      p.chunk_values == 4096 && p.verify_write);
    SMESH_TEST_EQ(CompressionOptions::parse("off", &p), SMESH_SUCCESS);
    SMESH_TEST_ASSERT(!p.enabled);
    SMESH_TEST_EQ(CompressionOptions::parse("pwrel:2", &p), SMESH_FAILURE);
    SMESH_TEST_EQ(CompressionOptions::parse("h:0.5", &p), SMESH_FAILURE);
    SMESH_TEST_EQ(CompressionOptions::parse("abs", &p), SMESH_FAILURE);
    SMESH_TEST_EQ(CompressionOptions::parse("lossless,foo", &p), SMESH_FAILURE);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    fs::remove_all(kDir);
    fs::create_directories(kDir);
    SMESH_RUN_TEST(test_roundtrip_types);
    SMESH_RUN_TEST(test_conversions);
    SMESH_RUN_TEST(test_lossy_files);
    SMESH_RUN_TEST(test_escape_fallback);
    SMESH_RUN_TEST(test_thread_independence);
    SMESH_RUN_TEST(test_corruption);
    SMESH_RUN_TEST(test_detect_files);
    SMESH_RUN_TEST(test_buffer_and_options);
    if (!err) fs::remove_all(kDir);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
