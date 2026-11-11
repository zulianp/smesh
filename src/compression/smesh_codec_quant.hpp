#ifndef SMESH_CODEC_QUANT_HPP
#define SMESH_CODEC_QUANT_HPP

// Error-bounded lossy codecs for floating-point values (the SZ recipe: quantise on a grid,
// predict, store the integer residuals losslessly), with two properties:
//
//  * Exact, portable decoding. Every grid step is a power of two, x' = T(origin + q * 2^L) with
//    an integer q: q * 2^L is exact, so the only rounding is one addition (fused or not), and
//    every CPU, GPU and compiler decodes the same bits.
//  * A guaranteed bound. The encoder computes each decoded value with the decoder's own
//    function, including the final cast to T, and checks |x' - x| <= bound in double. A value
//    that fails, and every NaN or Inf, is stored verbatim ("escaped").
//
// QuantFixed (BoundMode Abs, RelRange): one step 2^L <= 2 eb_abs for the whole array, origin
// the global minimum. Chunk layout:
//   [n_esc u32][escape positions: words][q: words (u32 or u64)][escaped values: T each]
// The words go through the FOR / DFOR block codec, so DFOR acts as a previous-value predictor.
//
// QuantVar (BoundMode PointwiseRel, NodeHeight): a level per value, L_i = floor(log2(2 b_i)),
// never below L_min (L_min keeps |q| <= 2^30 for the whole array). Chunk layout:
//   [lv: words, 0 = escaped, else L_i - L_min + 1][r: words][escaped values: T each]
// with r_i = zz(q_i - predict(q_{i-1} rescaled to L_i)) in 32-bit modular arithmetic.
// The previous value in the stream is the predictor; escaped values leave it unchanged.

#include "smesh_codec_int.hpp"

#include <cmath>
#include <vector>

namespace smesh {
namespace compression {

/// Parameters of a lossy encoding, fixed for the whole array (all chunks, all ranks).
struct QuantParams {
    BoundMode mode{BoundMode::Lossless};
    double    eb{0};       // bound as requested
    double    eb_abs{0};   // QuantFixed: absolute bound used (eb, or eb * range)
    double    origin{0};   // grid origin
    i32       L{0};        // QuantFixed: level of the step
    i32       L_min{0};    // QuantVar: finest level
    u8        wide{0};     // QuantFixed: q stored as u64 words
};

/// What an encoding achieved (exact maxima and counts, so they can be reduced in any order).
struct QuantStats {
    double max_abs_err{0};  // max |x' - x| over the finite values that were quantised
    double max_ratio{0};    // max |x' - x| / bound
    u64    n_escaped{0};

    void merge(const QuantStats &o) {
        max_abs_err = o.max_abs_err > max_abs_err ? o.max_abs_err : max_abs_err;
        max_ratio   = o.max_ratio > max_ratio ? o.max_ratio : max_ratio;
        n_escaped += o.n_escaped;
    }
};

/// Grid levels are kept in [kLevelLo, kLevelHi], so 2^L and 2^-L are finite doubles.
static constexpr int kLevelLo = -1000;
static constexpr int kLevelHi = 1000;

/// 2^L, exactly.
SMESH_HOST_DEVICE SMESH_INLINE double pow2i(int L) { return ldexp(1.0, L); }

/// floor(log2(2 b)) computed exactly (frexp), for b > 0 finite; INT32_MIN otherwise.
SMESH_HOST_DEVICE SMESH_INLINE i32 level_for_bound(double b) {
    if (!(b > 0) || !std::isfinite(b)) return INT32_MIN;
    int e = 0;
    (void)frexp(2 * b, &e);  // 2b = m 2^e, m in [0.5, 1)
    return i32(e - 1);
}

/// The decoded value of grid index q at level L (shared by encoder and decoder).
template <class T>
SMESH_HOST_DEVICE SMESH_INLINE T quant_value(double origin, i64 q, int L) {
    return T(origin + double(q) * pow2i(L));
}

/// Previous index rescaled from level Lp to level L, in integers only (encoder and decoder agree
/// bit for bit); clamped to [-2^30, 2^30].
SMESH_HOST_DEVICE SMESH_INLINE i64 quant_predict(i64 q_prev, int Lp, int L) {
    constexpr i64 lim = i64(1) << 30;
    const int     s   = Lp - L;
    if (s >= 0) {
        if (q_prev == 0) return 0;
        if (s >= 31) return q_prev > 0 ? lim : -lim;
        const i64 p = q_prev * (i64(1) << s);
        return p > lim ? lim : (p < -lim ? -lim : p);
    }
    const int r = -s;
    if (r >= 62) return 0;
    const i64 num = q_prev + (i64(1) << (r - 1));  // round half up: floor((q + 2^(r-1)) / 2^r)
    return num >= 0 ? (num >> r) : -((-num + (i64(1) << r) - 1) >> r);
}

template <class T>
SMESH_HOST_DEVICE SMESH_INLINE bool quant_finite(T x) {
    return std::isfinite(double(x));
}

/// QuantFixed parameters for |x' - x| <= eb_abs on values in [vmin, vmax] (finite). Returns
/// false if the bound is not representable (eb_abs <= 0, or more than 2^52 steps in the range):
/// the caller then encodes losslessly.
static inline bool plan_quant_fixed(double eb_abs, double vmin, double vmax, QuantParams &p) {
    const i32 L = level_for_bound(eb_abs);
    if (L == INT32_MIN || L < kLevelLo || L > kLevelHi || !(vmax >= vmin)) return false;
    const double steps = (vmax - vmin) * pow2i(-L) + 1;
    if (!(steps < 4503599627370496.0)) return false;  // 2^52
    p.eb_abs = eb_abs;
    p.origin = vmin;
    p.L      = L;
    p.wide   = steps >= 4294967295.0 ? 1 : 0;
    return true;
}

/// QuantVar parameters: values within amax of `origin`, smallest positive bound bmin (0 if
/// none). L_min is the coarsest level that still keeps |q| <= 2^30 and serves the smallest
/// bound; a value whose bound needs a finer level is escaped.
static inline void plan_quant_var(double origin, double amax, double bmin, QuantParams &p) {
    int L_range = kLevelLo;
    if (amax > 0 && std::isfinite(amax)) {
        int e = 0;
        (void)frexp(amax, &e);  // amax < 2^e, so |x - origin| / 2^(e-30) < 2^30
        L_range = e - 30;
    }
    const i32 Lb = level_for_bound(bmin);
    int       L  = Lb == INT32_MIN ? L_range : (int(Lb) > L_range ? int(Lb) : L_range);
    L            = L < kLevelLo ? kLevelLo : (L > kLevelHi ? kLevelHi : L);
    p.origin     = origin;
    p.L_min      = L;
}

namespace detail {
template <class T>
static inline double quant_err(T xd, T x) {
    return std::fabs(double(xd) - double(x));
}
template <class T>
static inline void track(QuantStats &st, double err, double bound) {
    st.max_abs_err = err > st.max_abs_err ? err : st.max_abs_err;
    const double ratio = bound > 0 ? err / bound : 0.0;
    st.max_ratio       = ratio > st.max_ratio ? ratio : st.max_ratio;
}
}  // namespace detail

// ------------------------------------------------------------------------------ QuantFixed
/// Encode n values with one step 2^p.L (p.origin, p.eb_abs, p.wide set by the caller).
template <class T>
static inline void enc_quant_fixed(ByteSink &out, const T *x, size_t n, const QuantParams &p, QuantStats &st) {
    std::vector<u32> esc;
    std::vector<u64> q(n);
    const double     inv = pow2i(-p.L);
    const double     qmax = p.wide ? 9007199254740992.0 : 4294967295.0;  // 2^53 or 2^32 - 1
    u64              prev = 0;
    for (size_t i = 0; i < n; i++) {
        bool ok = quant_finite(x[i]);
        if (ok) {
            const double qd = std::nearbyint((double(x[i]) - p.origin) * inv);
            ok              = qd >= 0 && qd <= qmax;
            if (ok) {
                const T      xd  = quant_value<T>(p.origin, i64(qd), p.L);
                const double err = detail::quant_err(xd, x[i]);
                ok               = err <= p.eb_abs && quant_finite(xd);
                if (ok) {
                    q[i] = u64(qd);
                    detail::track<T>(st, err, p.eb_abs);
                }
            }
        }
        if (!ok) {
            esc.push_back(u32(i));
            q[i] = prev;  // costs nothing in DFOR
        }
        prev = q[i];
    }
    st.n_escaped += esc.size();
    out.put<u32>(u32(esc.size()));
    enc_words<u32>(out, esc.data(), esc.size());
    if (p.wide) {
        enc_words<u64>(out, q.data(), n);
    } else {
        std::vector<u32> q32(q.begin(), q.end());
        enc_words<u32>(out, q32.data(), n);
    }
    for (u32 i : esc) out.put<T>(x[i]);
}

/// Decode n values. Scratch: qw of n words (u32, or u64 if p.wide: pass both, one is used),
/// pos of n u32.
template <class T, bool Checked>
SMESH_HOST_DEVICE SMESH_INLINE const u8 *dec_quant_fixed(const u8 *p, const u8 *end, size_t n, const QuantParams &qp,
                                                         T *SMESH_RESTRICT out, u32 *SMESH_RESTRICT q32,
                                                         u64 *SMESH_RESTRICT q64, u32 *SMESH_RESTRICT pos) {
    if constexpr (Checked) {
        if (size_t(end - p) < 4) return nullptr;
    }
    const u32 n_esc = take<u32>(p);
    if constexpr (Checked) {
        if (n_esc > n) return nullptr;
    }
    p = dec_words<u32, Checked>(p, end, n_esc, pos);
    if constexpr (Checked) {
        if (!p) return nullptr;
    }
    if (qp.wide) {
        p = dec_words<u64, Checked>(p, end, n, q64);
        if constexpr (Checked) {
            if (!p) return nullptr;
        }
        for (size_t i = 0; i < n; i++) out[i] = quant_value<T>(qp.origin, i64(q64[i]), qp.L);
    } else {
        p = dec_words<u32, Checked>(p, end, n, q32);
        if constexpr (Checked) {
            if (!p) return nullptr;
        }
        for (size_t i = 0; i < n; i++) out[i] = quant_value<T>(qp.origin, i64(q32[i]), qp.L);
    }
    if constexpr (Checked) {
        if (size_t(end - p) < size_t(n_esc) * sizeof(T)) return nullptr;
    }
    for (u32 k = 0; k < n_esc; k++) {
        if constexpr (Checked) {
            if (pos[k] >= n) return nullptr;
        }
        memcpy(&out[pos[k]], p, sizeof(T));
        p += sizeof(T);
    }
    return p;
}

// ------------------------------------------------------------------------------ QuantVar
/// Encode n values x[src(i)] with per-value bounds bound(i) (both functors of the position i in
/// the stream, so sfem can encode a pack's node list from global arrays).
template <class T, class Src, class Bound>
static inline void enc_quant_var(ByteSink &out, const T *x, size_t n, const QuantParams &p, Src src, Bound bound,
                                 QuantStats &st) {
    std::vector<u32> lv(n), r(n);
    std::vector<T>   esc;
    i64              q_prev = 0;
    int              L_prev = p.L_min;
    for (size_t i = 0; i < n; i++) {
        const T      xi = x[src(i)];
        const double b  = bound(i);
        bool         ok = quant_finite(xi) && b >= 0 && std::isfinite(b);
        int          L  = p.L_min;
        i64          q  = 0;
        if (ok) {
            const i32 Lb = level_for_bound(b);
            const int L0 = Lb < p.L_min ? p.L_min : (Lb > kLevelHi ? kLevelHi : int(Lb));
            ok           = false;
            // The rounding error is at most half a step; the cast to T can add up to half an ulp,
            // so near T's resolution a finer level may still meet the bound.
            for (int k = 0; k < 3 && !ok && L0 - k >= p.L_min; k++) {
                L               = L0 - k;
                const double qd = std::nearbyint((double(xi) - p.origin) * pow2i(-L));
                if (!(std::fabs(qd) <= 1073741824.0)) break;  // 2^30
                q                = i64(qd);
                const T      xd  = quant_value<T>(p.origin, q, L);
                const double err = detail::quant_err(xd, xi);
                ok = err <= b && quant_finite(xd) &&
                     (p.mode != BoundMode::PointwiseRel || std::signbit(double(xd)) == std::signbit(double(xi)));
                if (ok) detail::track<T>(st, err, b);
            }
        }
        if (!ok) {
            lv[i] = 0;
            r[i]  = 0;
            esc.push_back(xi);
            continue;
        }
        lv[i]        = u32(L - p.L_min + 1);
        const i64 pr = quant_predict(q_prev, L_prev, L);
        r[i]         = zz_enc<u32>(u32(u64(q) - u64(pr)));
        q_prev       = q;
        L_prev       = L;
    }
    st.n_escaped += esc.size();
    enc_words<u32>(out, lv.data(), n);
    enc_words<u32>(out, r.data(), n);
    out.append(esc.data(), esc.size() * sizeof(T));
}

/// Encode n contiguous values with per-value bounds b[i].
template <class T>
static inline void enc_quant_var(ByteSink &out, const T *x, size_t n, const QuantParams &p, const double *b,
                                 QuantStats &st) {
    enc_quant_var<T>(out, x, n, p, [](size_t i) { return i; }, [b](size_t i) { return b[i]; }, st);
}

/// Decode n values. Scratch: lv and r of n u32 each.
template <class T, bool Checked>
SMESH_HOST_DEVICE SMESH_INLINE const u8 *dec_quant_var(const u8 *p, const u8 *end, size_t n, const QuantParams &qp,
                                                       T *SMESH_RESTRICT out, u32 *SMESH_RESTRICT lv,
                                                       u32 *SMESH_RESTRICT r) {
    p = dec_words<u32, Checked>(p, end, n, lv);
    if constexpr (Checked) {
        if (!p) return nullptr;
    }
    p = dec_words<u32, Checked>(p, end, n, r);
    if constexpr (Checked) {
        if (!p) return nullptr;
    }
    size_t n_esc = 0;
    for (size_t i = 0; i < n; i++) n_esc += lv[i] == 0;
    if constexpr (Checked) {
        if (size_t(end - p) < n_esc * sizeof(T)) return nullptr;
        const u32 lv_max = u32(kLevelHi - qp.L_min + 1);
        for (size_t i = 0; i < n; i++)
            if (lv[i] > lv_max) return nullptr;
    }
    i64 q_prev = 0;
    int L_prev = qp.L_min;
    for (size_t i = 0; i < n; i++) {
        if (lv[i] == 0) {
            memcpy(&out[i], p, sizeof(T));
            p += sizeof(T);
            continue;
        }
        const int L  = qp.L_min + int(lv[i]) - 1;
        const i64 pr = quant_predict(q_prev, L_prev, L);
        const i64 q  = i64(i32(u32(u64(pr) + u64(zz_dec<u32>(r[i])))));
        out[i]       = quant_value<T>(qp.origin, q, L);
        q_prev       = q;
        L_prev       = L;
    }
    return p;
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CODEC_QUANT_HPP
