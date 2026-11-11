// Error-bounded quantisers (smesh_codec_quant.hpp): the bound must hold for every value, for
// every mode, value type and bound, on inputs chosen to break it: values on grid midpoints,
// huge dynamic range, values near the largest float, bounds below one ulp, NaN / Inf mixed in,
// constants, signed zeros. Non-finite values must come back bit-exact; the statistics the
// encoder reports must match a recomputation from the decoded values.

#include "smesh_codec_quant.hpp"
#include "smesh_test.hpp"

#include <cmath>
#include <limits>
#include <random>
#include <vector>

using namespace smesh;
using namespace smesh::compression;

static const double kBounds[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12};

template <class T>
static std::vector<T> make_input(int kind, size_t n, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<T>                         v(n);
    const double                           big = double(std::numeric_limits<T>::max());
    for (size_t i = 0; i < n; i++) {
        const double t = double(i) / double(n);
        switch (kind) {
            case 0: v[i] = T(std::sin(20 * t) + 0.1 * std::cos(300 * t)); break;  // smooth
            case 1: v[i] = T(u(rng) - 0.5); break;                                // noise
            case 2: v[i] = T((double(i % 97) + 0.5) * 0.0625); break;             // grid midpoints
            case 3: v[i] = T(std::pow(10.0, -30 + 60 * u(rng)) * (u(rng) < 0.5 ? -1 : 1)); break;  // 1e-30..1e30
            case 4: v[i] = T(big * (0.9 + 0.1 * u(rng))); break;                  // near the largest value
            case 5: v[i] = T(3.25); break;                                        // constant
            case 6: {                                                             // specials mixed in
                const double r = u(rng);
                v[i]           = r < 0.05   ? std::numeric_limits<T>::quiet_NaN()
                                 : r < 0.1  ? std::numeric_limits<T>::infinity()
                                 : r < 0.15 ? -std::numeric_limits<T>::infinity()
                                 : r < 0.2  ? T(-0.0)
                                 : r < 0.25 ? T(0.0)
                                            : T(std::sin(7 * t));
                break;
            }
            default: v[i] = T(0);
        }
    }
    return v;
}

template <class T>
static void finite_range(const std::vector<T> &v, double *lo, double *hi, double *amax, double *amin_nz) {
    *lo = INFINITY, *hi = -INFINITY, *amax = 0, *amin_nz = INFINITY;
    for (T x : v) {
        const double d = double(x);
        if (!std::isfinite(d)) continue;
        *lo   = std::min(*lo, d);
        *hi   = std::max(*hi, d);
        *amax = std::max(*amax, std::fabs(d));
        if (d != 0) *amin_nz = std::min(*amin_nz, std::fabs(d));
    }
}

// Check decoded values against bounds b(i); non-finite inputs bit-exact. Returns the recomputed
// max abs error over finite values.
template <class T, class B>
static int check_decoded(const std::vector<T> &x, const std::vector<T> &y, B bound, double *max_err, double *max_ratio) {
    *max_err = 0, *max_ratio = 0;
    for (size_t i = 0; i < x.size(); i++) {
        const double xi = double(x[i]), yi = double(y[i]);
        if (!std::isfinite(xi)) {
            SMESH_TEST_ASSERT(memcmp(&x[i], &y[i], sizeof(T)) == 0);
            continue;
        }
        const double err = std::fabs(yi - xi), b = bound(i);
        if (!(err <= b)) {
            fprintf(stderr, "bound violated at %zu: x=%.17g y=%.17g err=%g bound=%g\n", i, xi, yi, err, b);
            return SMESH_TEST_FAILURE;
        }
        *max_err   = std::max(*max_err, err);
        *max_ratio = std::max(*max_ratio, b > 0 ? err / b : 0.0);
    }
    return SMESH_TEST_SUCCESS;
}

template <class T>
static int check_fixed(BoundMode mode) {
    std::mt19937_64 rng(5);
    size_t          n_cases = 0, n_lossless = 0;
    for (int kind = 0; kind < 7; kind++) {
        for (size_t n : {size_t(1), size_t(129), size_t(5000)}) {
            const std::vector<T> x = make_input<T>(kind, n, rng);
            double lo, hi, amax, amin_nz;
            finite_range(x, &lo, &hi, &amax, &amin_nz);
            for (double eb : kBounds) {
                n_cases++;
                const double eb_abs = mode == BoundMode::RelRange ? eb * (hi - lo) : eb;
                QuantParams  p;
                p.mode = mode;
                p.eb   = eb;
                if (!std::isfinite(lo) || !plan_quant_fixed(eb_abs, lo, hi, p)) {
                    n_lossless++;  // the chunked layer writes such arrays losslessly
                    continue;
                }
                QuantStats st;
                ByteSink   s;
                enc_quant_fixed<T>(s, x.data(), n, p, st);
                const size_t enc = s.size();
                s.bytes.resize(enc + kSlack);
                std::vector<T>   y(n);
                std::vector<u32> q32(n), pos(n);
                std::vector<u64> q64(n);
                const u8        *end = s.bytes.data() + enc;
                SMESH_TEST_ASSERT((dec_quant_fixed<T, true>(s.bytes.data(), end, n, p, y.data(), q32.data(), q64.data(),
                                                            pos.data()) == end));
                double max_err, max_ratio;
                if (check_decoded(x, y, [&](size_t) { return p.eb_abs; }, &max_err, &max_ratio)) {
                    fprintf(stderr, "fixed %s kind=%d n=%zu eb=%g\n", to_string(mode), kind, n, eb);
                    return SMESH_TEST_FAILURE;
                }
                SMESH_TEST_ASSERT(max_err == st.max_abs_err);
                SMESH_TEST_ASSERT(st.max_ratio <= 1.0);
                SMESH_TEST_ASSERT(st.n_escaped <= n);
                // Truncation is detected.
                if (enc > 4) {
                    SMESH_TEST_ASSERT((dec_quant_fixed<T, true>(s.bytes.data(), s.bytes.data() + enc - 1, n, p, y.data(),
                                                                q32.data(), q64.data(), pos.data()) == nullptr));
                }
            }
        }
    }
    printf("fixed %-6s %s: %zu cases, %zu not representable (written losslessly)\n", to_string(mode),
           to_string(DTypeOf<T>::value), n_cases, n_lossless);
    return SMESH_TEST_SUCCESS;
}

template <class T>
static int check_var(BoundMode mode) {
    std::mt19937_64 rng(9);
    for (int kind = 0; kind < 7; kind++) {
        for (size_t n : {size_t(1), size_t(129), size_t(5000)}) {
            const std::vector<T> x = make_input<T>(kind, n, rng);
            double lo, hi, amax, amin_nz;
            finite_range(x, &lo, &hi, &amax, &amin_nz);
            // NodeHeight-like bounds: a per-value "element size" spanning 4 orders of magnitude.
            std::vector<double> h(n);
            for (size_t i = 0; i < n; i++) h[i] = std::pow(10.0, -4.0 * double((i * 7919) % 1000) / 1000.0);
            for (double eb : kBounds) {
                std::vector<double> b(n);
                double              bmin = INFINITY;
                for (size_t i = 0; i < n; i++) {
                    b[i] = mode == BoundMode::PointwiseRel ? eb * std::fabs(double(x[i])) : eb * h[i];
                    if (b[i] > 0 && std::isfinite(b[i])) bmin = std::min(bmin, b[i]);
                }
                QuantParams p;
                p.mode = mode;
                p.eb   = eb;
                const double origin = mode == BoundMode::PointwiseRel ? 0.0 : (std::isfinite(lo) ? lo : 0.0);
                const double reach  = mode == BoundMode::PointwiseRel ? amax : (std::isfinite(lo) ? hi - lo : 0.0);
                plan_quant_var(origin, reach, std::isfinite(bmin) ? bmin : 0.0, p);
                QuantStats st;
                ByteSink   s;
                enc_quant_var<T>(s, x.data(), n, p, b.data(), st);
                const size_t enc = s.size();
                s.bytes.resize(enc + kSlack);
                std::vector<T>   y(n);
                std::vector<u32> lv(n), r(n);
                const u8        *end = s.bytes.data() + enc;
                SMESH_TEST_ASSERT((dec_quant_var<T, true>(s.bytes.data(), end, n, p, y.data(), lv.data(), r.data()) == end));
                double max_err, max_ratio;
                auto   bound = [&](size_t i) {
                    // Exact zeros and escapes have zero error; NaN bounds only occur for non-finite x.
                    return std::isfinite(b[i]) ? b[i] : 0.0;
                };
                if (check_decoded(x, y, bound, &max_err, &max_ratio)) {
                    fprintf(stderr, "var %s kind=%d n=%zu eb=%g\n", to_string(mode), kind, n, eb);
                    return SMESH_TEST_FAILURE;
                }
                SMESH_TEST_ASSERT(max_err == st.max_abs_err);
                if (mode == BoundMode::PointwiseRel) {
                    for (size_t i = 0; i < n; i++) {  // signs (and signed zeros) preserved
                        if (std::isfinite(double(x[i])))
                            SMESH_TEST_ASSERT(std::signbit(double(x[i])) == std::signbit(double(y[i])));
                    }
                }
                if (enc > 0) {
                    SMESH_TEST_ASSERT((dec_quant_var<T, true>(s.bytes.data(), s.bytes.data() + enc - 1, n, p, y.data(),
                                                              lv.data(), r.data()) == nullptr));
                }
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_level_for_bound() {
    SMESH_TEST_EQ(level_for_bound(0.5), i32(0));        // step 1 <= 2 * 0.5
    SMESH_TEST_EQ(level_for_bound(0.49), i32(-1));      // step 0.5
    SMESH_TEST_EQ(level_for_bound(1.0), i32(1));        // step 2
    SMESH_TEST_EQ(level_for_bound(0.0), i32(INT32_MIN));
    SMESH_TEST_EQ(level_for_bound(-1.0), i32(INT32_MIN));
    for (double b = 1e-300; b < 1e300; b *= 3.7) {
        const i32 L = level_for_bound(b);
        SMESH_TEST_ASSERT(pow2i(L) <= 2 * b && 2 * b < pow2i(L + 1));
    }
    return SMESH_TEST_SUCCESS;
}

static int test_predict() {
    // Exact for shifts up, round half up for shifts down, clamped to 2^30.
    SMESH_TEST_EQ(quant_predict(5, 3, 1), i64(20));
    SMESH_TEST_EQ(quant_predict(5, 1, 2), i64(3));    // 2.5 -> 3
    SMESH_TEST_EQ(quant_predict(-5, 1, 2), i64(-2));  // -2.5 -> -2 (half up)
    SMESH_TEST_EQ(quant_predict(-7, 0, 2), i64(-2));  // -1.75 -> -2
    SMESH_TEST_EQ(quant_predict(1, 0, -40), i64(1) << 30);
    SMESH_TEST_EQ(quant_predict(-1, 0, -40), -(i64(1) << 30));
    return SMESH_TEST_SUCCESS;
}

static int test_fixed() {
    if (check_fixed<f32>(BoundMode::Abs) || check_fixed<f64>(BoundMode::Abs) || check_fixed<f32>(BoundMode::RelRange) ||
        check_fixed<f64>(BoundMode::RelRange))
        return SMESH_TEST_FAILURE;
    return SMESH_TEST_SUCCESS;
}

static int test_var() {
    if (check_var<f32>(BoundMode::PointwiseRel) || check_var<f64>(BoundMode::PointwiseRel) ||
        check_var<f32>(BoundMode::NodeHeight) || check_var<f64>(BoundMode::NodeHeight))
        return SMESH_TEST_FAILURE;
    return SMESH_TEST_SUCCESS;
}

// Compression is real on smooth data: f32 sine at eb 1e-3 (relative to range) takes well under
// the 32 bits per value.
static int test_ratio_smooth() {
    const size_t   n = 1 << 16;
    std::vector<f32> x(n);
    for (size_t i = 0; i < n; i++) x[i] = f32(std::sin(0.001 * double(i)));
    QuantParams p;
    p.mode = BoundMode::RelRange;
    SMESH_TEST_ASSERT(plan_quant_fixed(1e-3 * 2.0, -1.0, 1.0, p));
    QuantStats st;
    ByteSink   s;
    enc_quant_fixed<f32>(s, x.data(), n, p, st);
    const double bits = 8.0 * double(s.size()) / double(n);
    printf("smooth f32, range-relative 1e-3: %.2f bits/value, max err/bound %.3f\n", bits, st.max_ratio);
    SMESH_TEST_ASSERT(bits < 8);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_level_for_bound);
    SMESH_RUN_TEST(test_predict);
    SMESH_RUN_TEST(test_fixed);
    SMESH_RUN_TEST(test_var);
    SMESH_RUN_TEST(test_ratio_smooth);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
