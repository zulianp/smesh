#ifndef SMESH_CODEC_INT_HPP
#define SMESH_CODEC_INT_HPP

// Lossless integer codecs on ordered words (smesh_compression_types.hpp: to_word / from_word),
// so they serve every value type: integers of any width and sign, and float bit patterns.
//
// A stream of n values is cut into blocks of kBlock = 128, each encoded as the smallest of
//   FOR    [0 u8][lo U][w u8][n x w bits: r_i = v_i - lo]
//   DFOR   [1 u8][first U][dlo U][w u8][(n-1) x w bits: r_i = zz(v_i - v_{i-1}) - dlo]
//   PFOR   [2 u8][lo U]              [patched r, n values]
//   PDFOR  [3 u8][first U][dlo U]    [patched r, n-1 values]
// with "patched r" = [w u8][n_exc u8][wh u8][m x w bits: low w bits of r][n_exc positions u8]
//                    [n_exc x wh bits: r >> w of the exceptions]
// (U = word_t<T>, u32 or u64). FOR suits clustered values, DFOR sorted or smooth ones; the
// patched forms keep a few outliers (row breaks of a mesh numbering, jumps) from widening the
// whole block. A constant block costs 1 + sizeof(U) + 1 bytes.
//
// Decoders are SMESH_HOST_DEVICE and allocation-free. With Checked = true they validate every
// block against `end` and return nullptr on malformed input instead of reading past it (the
// buffer must still have kSlack readable bytes after `end`); Checked = false is for trusted
// data in hot loops.

#include "smesh_codec_bits.hpp"

namespace smesh {
namespace compression {

namespace detail {
static constexpr u8 kTagFor   = 0;
static constexpr u8 kTagDFor  = 1;
static constexpr u8 kTagPFor  = 2;
static constexpr u8 kTagPDFor = 3;

/// Best patch width for m residuals: returns the encoded size of the patched part and sets
/// w, wh and the number of exceptions.
template <class U>
static inline size_t plan_patch(const U *r, size_t m, int *w_out, int *wh_out, size_t *exc_out) {
    size_t hist[65] = {0};
    int    maxbw    = 0;
    for (size_t i = 0; i < m; i++) {
        const int b = bit_width<U>(r[i]);
        hist[b]++;
        maxbw = b > maxbw ? b : maxbw;
    }
    size_t best = SIZE_MAX, exc = 0;
    for (int w = maxbw; w >= 0; w--) {
        if (w < maxbw) exc += hist[w + 1];
        if (exc > 255) break;
        const int    wh   = exc ? maxbw - w : 0;
        const size_t cost = 3 + packed_bytes(m, w) + exc + packed_bytes(exc, wh);
        if (cost < best) {
            best     = cost;
            *w_out   = w;
            *wh_out  = wh;
            *exc_out = exc;
        }
    }
    return best;
}

template <class U>
static inline void put_patched(ByteSink &out, const U *r, size_t m, int w, int wh, size_t n_exc) {
    U   low[kBlock], high[kBlock];
    u8  pos[kBlock];
    size_t k = 0;
    const U mask = mask_of<U>(w);
    for (size_t i = 0; i < m; i++) {
        low[i] = U(r[i] & mask);
        if (w < kWordBits<U> && (r[i] >> w) != 0) {
            pos[k]  = u8(i);
            high[k] = U(r[i] >> w);
            k++;
        }
    }
    out.put<u8>(u8(w));
    out.put<u8>(u8(n_exc));
    out.put<u8>(u8(wh));
    out.pack<U>(low, m, w);
    out.append(pos, k);
    out.pack<U>(high, k, wh);
}

/// Decode a patched residual block of m values into out (low bits | exceptions).
template <class U, bool Checked>
SMESH_HOST_DEVICE SMESH_INLINE const u8 *get_patched(const u8 *p, const u8 *end, size_t m, U *SMESH_RESTRICT out) {
    if constexpr (Checked) {
        if (size_t(end - p) < 3) return nullptr;
    }
    const int    w     = take<u8>(p);
    const size_t n_exc = take<u8>(p);
    const int    wh    = take<u8>(p);
    if constexpr (Checked) {
        if (w > kWordBits<U> || wh > kWordBits<U> || w + wh > kWordBits<U> || n_exc > m ||
            size_t(end - p) < packed_bytes(m, w) + n_exc + packed_bytes(n_exc, wh))
            return nullptr;
    }
    p += bp_unpack<U>(p, m, w, out);
    const u8 *pos = p;
    p += n_exc;
    U high[kBlock];
    p += bp_unpack<U>(p, n_exc, wh, high);
    for (size_t k = 0; k < n_exc; k++) {
        if constexpr (Checked) {
            if (pos[k] >= m) return nullptr;
        }
        out[pos[k]] = U(out[pos[k]] | U(high[k] << w));
    }
    return p;
}

/// Encoded size of a FOR / DFOR block (header + packed bits).
template <class U>
static inline size_t for_block_bytes(size_t n, int w) {
    return 1 + sizeof(U) + 1 + packed_bytes(n, w);
}
template <class U>
static inline size_t dfor_block_bytes(size_t n, int w) {
    return 1 + 2 * sizeof(U) + 1 + packed_bytes(n - 1, w);
}
}  // namespace detail

/// Encode one block of 1 <= n <= kBlock words (the smaller of FOR and DFOR).
template <class U>
static inline void enc_word_block(ByteSink &out, const U *v, size_t n) {
    U lo = v[0], hi = v[0];
    for (size_t i = 1; i < n; i++) {
        lo = v[i] < lo ? v[i] : lo;
        hi = v[i] > hi ? v[i] : hi;
    }
    const int w_for = bit_width<U>(U(hi - lo));
    U         rf[kBlock];
    for (size_t i = 0; i < n; i++) rf[i] = U(v[i] - lo);

    U   d[kBlock];
    U   dlo = 0, dhi = 0;
    int w_dfor = 0;
    if (n > 1) {
        for (size_t i = 1; i < n; i++) d[i - 1] = zz_enc<U>(U(v[i] - v[i - 1]));
        dlo = dhi = d[0];
        for (size_t i = 1; i + 1 < n; i++) {
            dlo = d[i] < dlo ? d[i] : dlo;
            dhi = d[i] > dhi ? d[i] : dhi;
        }
        w_dfor = bit_width<U>(U(dhi - dlo));
        for (size_t i = 0; i + 1 < n; i++) d[i] = U(d[i] - dlo);
    }

    // The smallest of the four forms.
    int          pw = 0, pwh = 0, dw = 0, dwh = 0;
    size_t       pexc = 0, dexc = 0;
    const size_t c_for   = detail::for_block_bytes<U>(n, w_for);
    const size_t c_pfor  = 1 + sizeof(U) + detail::plan_patch<U>(rf, n, &pw, &pwh, &pexc);
    const size_t c_dfor  = n > 1 ? detail::dfor_block_bytes<U>(n, w_dfor) : SIZE_MAX;
    const size_t c_pdfor = n > 1 ? 1 + 2 * sizeof(U) + detail::plan_patch<U>(d, n - 1, &dw, &dwh, &dexc) : SIZE_MAX;
    size_t       best    = c_for;
    u8           tag     = detail::kTagFor;
    if (c_dfor < best) best = c_dfor, tag = detail::kTagDFor;
    if (c_pfor < best) best = c_pfor, tag = detail::kTagPFor;
    if (c_pdfor < best) best = c_pdfor, tag = detail::kTagPDFor;

    out.put<u8>(tag);
    switch (tag) {
        case detail::kTagFor:
            out.put<U>(lo);
            out.put<u8>(u8(w_for));
            out.pack<U>(rf, n, w_for);
            break;
        case detail::kTagDFor:
            out.put<U>(v[0]);
            out.put<U>(dlo);
            out.put<u8>(u8(w_dfor));
            out.pack<U>(d, n - 1, w_dfor);
            break;
        case detail::kTagPFor:
            out.put<U>(lo);
            detail::put_patched<U>(out, rf, n, pw, pwh, pexc);
            break;
        default:
            out.put<U>(v[0]);
            out.put<U>(dlo);
            detail::put_patched<U>(out, d, n - 1, dw, dwh, dexc);
            break;
    }
}

/// Decode one block of n words into out (out holds >= n words). Returns the position after the
/// block, or nullptr (Checked) if the block is malformed or runs past `end`.
template <class U, bool Checked>
SMESH_HOST_DEVICE SMESH_INLINE const u8 *dec_word_block(const u8 *p, const u8 *end, size_t n, U *SMESH_RESTRICT out) {
    if constexpr (Checked) {
        if (p == nullptr || size_t(end - p) < 1 + sizeof(U) + 1) return nullptr;
    }
    const u8 tag = take<u8>(p);
    if (tag == detail::kTagPFor) {
        const U lo = take<U>(p);
        p          = detail::get_patched<U, Checked>(p, end, n, out);
        if constexpr (Checked) {
            if (!p) return nullptr;
        }
        for (size_t i = 0; i < n; i++) out[i] = U(out[i] + lo);
        return p;
    }
    if (tag == detail::kTagPDFor) {
        if constexpr (Checked) {
            if (size_t(end - p) < 2 * sizeof(U)) return nullptr;
        }
        const U first = take<U>(p);
        const U dlo   = take<U>(p);
        out[0]        = first;
        if (n > 1) {
            p = detail::get_patched<U, Checked>(p, end, n - 1, out + 1);
            if constexpr (Checked) {
                if (!p) return nullptr;
            }
            U acc = first;
            for (size_t i = 1; i < n; i++) {
                acc    = U(acc + zz_dec<U>(U(out[i] + dlo)));
                out[i] = acc;
            }
        }
        return p;
    }
    if (tag == detail::kTagFor) {
        const U   lo = take<U>(p);
        const int w  = take<u8>(p);
        if constexpr (Checked) {
            if (w > kWordBits<U> || size_t(end - p) < packed_bytes(n, w)) return nullptr;
        }
        p += bp_unpack<U>(p, n, w, out);
        for (size_t i = 0; i < n; i++) out[i] = U(out[i] + lo);
        return p;
    }
    if constexpr (Checked) {
        if (tag != detail::kTagDFor || size_t(end - p) < 2 * sizeof(U) + 1) return nullptr;
    }
    const U   first = take<U>(p);
    const U   dlo   = take<U>(p);
    const int w     = take<u8>(p);
    if constexpr (Checked) {
        if (w > kWordBits<U> || size_t(end - p) < packed_bytes(n - 1, w)) return nullptr;
    }
    out[0] = first;
    if (n == 1) return p;
    p += bp_unpack<U>(p, n - 1, w, out + 1);
    U acc = first;
    for (size_t i = 1; i < n; i++) {
        acc    = U(acc + zz_dec<U>(U(out[i] + dlo)));
        out[i] = acc;
    }
    return p;
}

/// Encode n values of any value type T (lossless).
template <class T>
static inline void enc_packed(ByteSink &out, const T *v, size_t n) {
    using U = word_t<T>;
    U w[kBlock];
    for (size_t b = 0; b < n; b += kBlock) {
        const size_t m = n - b < kBlock ? n - b : kBlock;
        for (size_t i = 0; i < m; i++) w[i] = to_word<T>(v[b + i]);
        enc_word_block<U>(out, w, m);
    }
}

/// Decode n values of type T written by enc_packed.
template <class T, bool Checked>
SMESH_HOST_DEVICE SMESH_INLINE const u8 *dec_packed(const u8 *p, const u8 *end, size_t n, T *SMESH_RESTRICT out) {
    using U = word_t<T>;
    U w[kBlock];
    for (size_t b = 0; b < n; b += kBlock) {
        const size_t m = n - b < kBlock ? n - b : kBlock;
        p              = dec_word_block<U, Checked>(p, end, m, w);
        if constexpr (Checked) {
            if (!p) return nullptr;
        }
        for (size_t i = 0; i < m; i++) out[b + i] = from_word<T>(w[i]);
    }
    return p;
}

/// Encode n words that the caller already mapped (e.g. quantizer residuals): FOR/DFOR blocks.
template <class U>
static inline void enc_words(ByteSink &out, const U *v, size_t n) {
    for (size_t b = 0; b < n; b += kBlock) enc_word_block<U>(out, v + b, n - b < kBlock ? n - b : kBlock);
}

template <class U, bool Checked>
SMESH_HOST_DEVICE SMESH_INLINE const u8 *dec_words(const u8 *p, const u8 *end, size_t n, U *SMESH_RESTRICT out) {
    for (size_t b = 0; b < n; b += kBlock) {
        p = dec_word_block<U, Checked>(p, end, n - b < kBlock ? n - b : kBlock, out + b);
        if constexpr (Checked) {
            if (!p) return nullptr;
        }
    }
    return p;
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CODEC_INT_HPP
