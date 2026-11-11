// Lossless kernels of src/compression: bit packing at every width, the ordered-word map, the
// FOR/DFOR block codec for every value type, the checked decoders on damaged input, byte
// shuffle, XXH64 against reference values, and the zstd / lz4 stage when built.
//
// Every round trip must be bit-exact; the checked decoders must return nullptr, never read past
// the end, on truncated or corrupted streams.

#include "smesh_byte_codec.hpp"
#include "smesh_codec_float.hpp"
#include "smesh_codec_hash.hpp"
#include "smesh_codec_int.hpp"
#include "smesh_test.hpp"

#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

using namespace smesh;
using namespace smesh::compression;

static const size_t kSizes[] = {0, 1, 2, 7, 8, 9, 127, 128, 129, 255, 1000, 65537};

template <class U>
static int check_bitpack() {
    std::mt19937_64 rng(42);
    for (int w = 0; w <= kWordBits<U>; w++) {
        for (size_t n : kSizes) {
            std::vector<U> v(n);
            const U        m = mask_of<U>(w);
            for (auto &x : v) x = U(rng()) & m;
            if (w > 0 && n > 0) v[0] = m;  // the widest value
            ByteSink s;
            s.pack<U>(v.data(), n, w);
            SMESH_TEST_EQ(s.size(), packed_bytes(n, w));
            s.bytes.resize(s.size() + kSlack);
            std::vector<U> out(n + 1, U(0x5a));
            const size_t   used = bp_unpack<U>(s.bytes.data(), n, w, out.data());
            SMESH_TEST_EQ(used, packed_bytes(n, w));
            for (size_t i = 0; i < n; i++) SMESH_TEST_EQ(out[i], v[i]);
            std::vector<U> gen(n + 1);
            bp_unpack_generic<U>(s.bytes.data(), n, w, gen.data());
            for (size_t i = 0; i < n; i++) SMESH_TEST_EQ(gen[i], v[i]);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_bitpack_u32() { return check_bitpack<u32>(); }
static int test_bitpack_u64() { return check_bitpack<u64>(); }

template <class T>
static bool bits_equal(const T &a, const T &b) {
    return memcmp(&a, &b, sizeof(T)) == 0;
}

// Ordered words: bijective on the bit patterns, and monotone for ordered (non-NaN) values.
template <class T>
static int check_word_map() {
    std::mt19937_64 rng(7);
    for (int k = 0; k < 20000; k++) {
        T a, b;
        u64 ra = rng(), rb = rng();
        memcpy(&a, &ra, sizeof(T));
        memcpy(&b, &rb, sizeof(T));
        SMESH_TEST_ASSERT(bits_equal(from_word<T>(to_word<T>(a)), a));
        if constexpr (is_float_v<T>) {
            const double da = double(a), db = double(b);
            if (std::isnan(da) || std::isnan(db) || da == db) continue;
            SMESH_TEST_ASSERT((da < db) == (to_word<T>(a) < to_word<T>(b)));
        } else {
            if (a == b) continue;
            SMESH_TEST_ASSERT((a < b) == (to_word<T>(a) < to_word<T>(b)));
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_word_map() {
    if (check_word_map<i8>() || check_word_map<u8>() || check_word_map<i16>() || check_word_map<u16>() ||
        check_word_map<i32>() || check_word_map<u32>() || check_word_map<i64>() || check_word_map<u64>() ||
        check_word_map<f16>() || check_word_map<f32>() || check_word_map<f64>())
        return SMESH_TEST_FAILURE;
    return SMESH_TEST_SUCCESS;
}

template <class T>
static std::vector<T> make_values(int kind, size_t n, std::mt19937_64 &rng) {
    std::vector<T> v(n);
    for (size_t i = 0; i < n; i++) {
        u64 r = rng();
        switch (kind) {
            case 0: memcpy(&v[i], &r, sizeof(T)); break;  // random bit patterns
            case 1: v[i] = T(i / 3); break;               // sorted with repeats
            case 2: v[i] = T(1000 + (r % 64)); break;     // clustered
            case 3: v[i] = T(7); break;                   // constant
            case 4: {                                     // extremes, alternating
                const T lo = std::numeric_limits<T>::lowest(), hi = std::numeric_limits<T>::max();
                v[i]       = (i & 1) ? hi : lo;
                break;
            }
            default: v[i] = T(0);
        }
    }
    return v;
}

template <class T>
static int check_packed_roundtrip() {
    std::mt19937_64 rng(11);
    for (int kind = 0; kind < 5; kind++) {
        for (size_t n : kSizes) {
            const std::vector<T> v = make_values<T>(kind, n, rng);
            ByteSink             s;
            enc_packed<T>(s, v.data(), n);
            const size_t enc = s.size();
            s.bytes.resize(enc + kSlack);
            std::vector<T> out(n + 1);
            const u8      *end = s.bytes.data() + enc;
            const u8      *p   = dec_packed<T, true>(s.bytes.data(), end, n, out.data());
            SMESH_TEST_ASSERT(p == end);
            for (size_t i = 0; i < n; i++) SMESH_TEST_ASSERT(bits_equal(out[i], v[i]));
            const u8 *q = dec_packed<T, false>(s.bytes.data(), end, n, out.data());
            SMESH_TEST_ASSERT(q == end);
            for (size_t i = 0; i < n; i++) SMESH_TEST_ASSERT(bits_equal(out[i], v[i]));
            if (kind == 3 && n == 65537) {  // constant: one header per block, no bits
                SMESH_TEST_ASSERT(enc == ((n + kBlock - 1) / kBlock) * (2 + sizeof(word_t<T>)));
            }
            // Truncated streams are rejected without reading past the end.
            if (n > 0) {
                for (size_t cut : {size_t(0), enc / 2, enc - 1}) {
                    std::vector<u8> t(s.bytes.begin(), s.bytes.begin() + ptrdiff_t(cut));
                    t.resize(cut + kSlack, 0xff);
                    SMESH_TEST_ASSERT((dec_packed<T, true>(t.data(), t.data() + cut, n, out.data()) == nullptr));
                }
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_packed_roundtrip() {
    if (check_packed_roundtrip<i8>() || check_packed_roundtrip<u8>() || check_packed_roundtrip<i16>() ||
        check_packed_roundtrip<u16>() || check_packed_roundtrip<i32>() || check_packed_roundtrip<u32>() ||
        check_packed_roundtrip<i64>() || check_packed_roundtrip<u64>() || check_packed_roundtrip<f32>() ||
        check_packed_roundtrip<f64>())
        return SMESH_TEST_FAILURE;
    return SMESH_TEST_SUCCESS;
}

// Floats: NaN payloads, signed zeros, infinities, denormals survive bit-exactly.
template <class T>
static int check_float_specials() {
    std::vector<T> v;
    using B = typename uint_of_size<sizeof(T)>::type;
    for (B bits : {B(0), B(B(1) << (8 * sizeof(T) - 1)), B(1), B(~B(0)), B(~B(0) >> 1), B(B(~B(0)) - 1)}) {
        T x;
        memcpy(&x, &bits, sizeof(T));
        v.push_back(x);
    }
    v.push_back(std::numeric_limits<T>::infinity());
    v.push_back(-std::numeric_limits<T>::infinity());
    v.push_back(std::numeric_limits<T>::quiet_NaN());
    v.push_back(std::numeric_limits<T>::denorm_min());
    v.push_back(-std::numeric_limits<T>::denorm_min());
    v.push_back(std::numeric_limits<T>::max());
    for (int i = 0; i < 300; i++) v.push_back(T(std::sin(0.01 * i)));
    ByteSink s;
    enc_packed<T>(s, v.data(), v.size());
    const size_t enc = s.size();
    s.bytes.resize(enc + kSlack);
    std::vector<T> out(v.size());
    SMESH_TEST_ASSERT((dec_packed<T, true>(s.bytes.data(), s.bytes.data() + enc, v.size(), out.data()) ==
                       s.bytes.data() + enc));
    for (size_t i = 0; i < v.size(); i++) SMESH_TEST_ASSERT(bits_equal(out[i], v[i]));
    return SMESH_TEST_SUCCESS;
}

static int test_float_specials() {
    if (check_float_specials<f32>() || check_float_specials<f64>()) return SMESH_TEST_FAILURE;
    return SMESH_TEST_SUCCESS;
}

// A corrupted width byte or tag is rejected.
static int test_corrupt_block() {
    std::vector<u32> v(300);
    for (size_t i = 0; i < v.size(); i++) v[i] = u32(i * 7);
    ByteSink s;
    enc_packed<u32>(s, v.data(), v.size());
    const size_t enc = s.size();
    s.bytes.resize(enc + kSlack);
    std::vector<u32> out(v.size());
    std::vector<u8>  bad = s.bytes;
    bad[0]               = 9;  // unknown tag
    SMESH_TEST_ASSERT((dec_packed<u32, true>(bad.data(), bad.data() + enc, v.size(), out.data()) == nullptr));
    bad    = s.bytes;
    bad[5] = 200;  // FOR width > 32 (tag 0, lo 4 bytes, then the width)
    if (s.bytes[0] == 0) {
        SMESH_TEST_ASSERT((dec_packed<u32, true>(bad.data(), bad.data() + enc, v.size(), out.data()) == nullptr));
    }
    return SMESH_TEST_SUCCESS;
}

static int test_byte_shuffle() {
    std::vector<f64> v(1001);
    for (size_t i = 0; i < v.size(); i++) v[i] = std::cos(0.001 * double(i));
    std::vector<u8> sh(v.size() * 8), back(v.size() * 8);
    byte_shuffle((const u8 *)v.data(), v.size(), 8, sh.data());
    byte_unshuffle(sh.data(), v.size(), 8, back.data());
    SMESH_TEST_ASSERT(memcmp(back.data(), v.data(), back.size()) == 0);
    return SMESH_TEST_SUCCESS;
}

// Reference values of XXH64 with seed 0.
static int test_xxh64() {
    SMESH_TEST_EQ(xxh64("", 0), u64(0xEF46DB3751D8E999ull));
    SMESH_TEST_EQ(xxh64("a", 1), u64(0xD24EC4F1A98C6E5Bull));
    SMESH_TEST_EQ(xxh64("abc", 3), u64(0x44BC2CF5AD770999ull));
    const char *s = "Nobody inspects the spammish repetition";
    SMESH_TEST_EQ(xxh64(s, strlen(s)), u64(0xFBCEA83C8A378BF1ull));
    return SMESH_TEST_SUCCESS;
}

static int test_byte_stages() {
    std::vector<u8> src(100000);
    for (size_t i = 0; i < src.size(); i++) src[i] = u8((i * i) >> 7);
    for (ByteStage st : {ByteStage::None, ByteStage::Zstd, ByteStage::Lz4}) {
        if (!byte_stage_available(st)) {
            fprintf(stderr, "byte stage %s not built: skipped\n", to_string(st));
            continue;
        }
        ByteEncoder     enc;
        ByteDecoder     dec;
        std::vector<u8> c;
        SMESH_TEST_EQ(enc.compress(st, 3, src.data(), src.size(), c), SMESH_SUCCESS);
        std::vector<u8> back(src.size());
        SMESH_TEST_EQ(dec.decompress(st, c.data(), c.size(), back.data(), back.size()), SMESH_SUCCESS);
        SMESH_TEST_ASSERT(back == src);
        if (st != ByteStage::None && !c.empty()) {
            c[c.size() / 2] ^= 0x55;  // damaged input must fail or at least not crash
            (void)dec.decompress(st, c.data(), c.size(), back.data(), back.size());
        }
    }
    return SMESH_TEST_SUCCESS;
}

// Decode speed of the packed codec on index-like int32 data (reported, not asserted).
static int test_decode_speed() {
    const size_t     n = size_t(1) << 22;
    std::vector<i32> v(n);
    std::mt19937     rng(3);
    for (size_t i = 0; i < n; i++) v[i] = i32(i / 4 + (rng() % 512));
    ByteSink s;
    enc_packed<i32>(s, v.data(), n);
    const size_t enc = s.size();
    s.bytes.resize(enc + kSlack);
    std::vector<i32> out(n);
    double           best = 1e30;
    for (int r = 0; r < 5; r++) {
        const auto t0 = std::chrono::steady_clock::now();
        dec_packed<i32, true>(s.bytes.data(), s.bytes.data() + enc, n, out.data());
        best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    SMESH_TEST_ASSERT(out == v);
    printf("packed int32 decode: %.2f GB/s (decoded bytes, one core), %.2f bits/value\n", 4.0 * n / best / 1e9,
           8.0 * double(enc) / double(n));
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_bitpack_u32);
    SMESH_RUN_TEST(test_bitpack_u64);
    SMESH_RUN_TEST(test_word_map);
    SMESH_RUN_TEST(test_packed_roundtrip);
    SMESH_RUN_TEST(test_float_specials);
    SMESH_RUN_TEST(test_corrupt_block);
    SMESH_RUN_TEST(test_byte_shuffle);
    SMESH_RUN_TEST(test_xxh64);
    SMESH_RUN_TEST(test_byte_stages);
    SMESH_RUN_TEST(test_decode_speed);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
