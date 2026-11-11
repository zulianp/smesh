#ifndef SMESH_CODEC_BITS_HPP
#define SMESH_CODEC_BITS_HPP

// Bit-level primitives: widths, zig-zag, bit packing of n values of w bits (little endian,
// value i at bits [i*w, (i+1)*w)), and a byte sink for the (host-only) encoders.
// Decoders read 64-bit words and may touch up to kSlack bytes past the packed data.

#include "smesh_compression_types.hpp"

#include <array>
#include <utility>
#include <vector>

namespace smesh {
namespace compression {

template <class U>
static constexpr int kWordBits = int(8 * sizeof(U));

template <class U>
SMESH_HOST_DEVICE SMESH_INLINE int bit_width(U x) {
#ifdef __CUDA_ARCH__
    if constexpr (sizeof(U) == 8) {
        return x ? 64 - __clzll((long long)x) : 0;
    } else {
        return x ? 32 - __clz((int)x) : 0;
    }
#else
    if constexpr (sizeof(U) == 8) {
        return x ? 64 - __builtin_clzll(x) : 0;
    } else {
        return x ? 32 - __builtin_clz(x) : 0;
    }
#endif
}

/// Zig-zag of a difference d, taken as a two's complement number of U's width.
template <class U>
SMESH_HOST_DEVICE SMESH_INLINE U zz_enc(U d) {
    using S = std::make_signed_t<U>;
    return U(d << 1) ^ U(S(d) >> (kWordBits<U> - 1));
}
template <class U>
SMESH_HOST_DEVICE SMESH_INLINE U zz_dec(U z) {
    return U(z >> 1) ^ U(U(0) - U(z & 1));
}

SMESH_HOST_DEVICE SMESH_INLINE size_t packed_bytes(size_t n, int w) { return (n * size_t(w) + 7) / 8; }

SMESH_HOST_DEVICE SMESH_INLINE u64 load_u64(const u8 *p) {
    u64 x;
    memcpy(&x, p, 8);
    return x;
}

template <class U>
SMESH_HOST_DEVICE SMESH_INLINE U mask_of(int w) {
    return w >= kWordBits<U> ? U(~U(0)) : U((U(1) << w) - 1);
}

/// Pack n values of width w (0 <= w <= bits of U) into dst, which must hold
/// packed_bytes(n, w) + kSlack zeroed bytes. Values must fit in w bits.
template <class U>
static inline void bp_pack(u8 *dst, const U *v, size_t n, int w) {
    if (w == 0) {
        return;
    }
    for (size_t i = 0; i < n; i++) {
        const size_t bit = i * size_t(w);
        const size_t b   = bit >> 3;
        const int    sh  = int(bit & 7);
        u64          x   = load_u64(dst + b);
        x |= u64(v[i]) << sh;
        memcpy(dst + b, &x, 8);
        if (sh + w > 64) {  // only u64 values wider than 56 bits spill into a ninth byte
            u64 y = load_u64(dst + b + 8);
            y |= u64(v[i]) >> (64 - sh);
            memcpy(dst + b + 8, &y, 8);
        }
    }
}

/// Generic unpacking for any width (device and u64 path).
template <class U>
SMESH_HOST_DEVICE SMESH_INLINE void bp_unpack_generic(const u8 *SMESH_RESTRICT src, size_t n, int w,
                                                      U *SMESH_RESTRICT out) {
    if (w == 0) {
        for (size_t i = 0; i < n; i++) out[i] = 0;
        return;
    }
    const U m = mask_of<U>(w);
    for (size_t i = 0; i < n; i++) {
        const size_t bit = i * size_t(w);
        const size_t b   = bit >> 3;
        const int    sh  = int(bit & 7);
        u64          x   = load_u64(src + b) >> sh;
        if (sh + w > 64) {
            x |= load_u64(src + b + 8) << (64 - sh);
        }
        out[i] = U(x) & m;
    }
}

#ifndef __CUDA_ARCH__
namespace detail {
// u32 unpacking specialised on the width: 8 values occupy exactly W bytes, so every shift is a
// compile-time constant and the loop vectorises.
template <int W>
static void unpack32_w(const u8 *SMESH_RESTRICT src, size_t n, u32 *SMESH_RESTRICT out) {
    if constexpr (W == 0) {
        for (size_t i = 0; i < n; i++) out[i] = 0;
    } else {
        constexpr u64 mask = W == 32 ? 0xffffffffull : ((1ull << W) - 1);
        size_t        i    = 0;
        for (; i + 8 <= n; i += 8) {
            const u8 *b = src + (i / 8) * W;
            for (int j = 0; j < 8; j++) {
                const int bit = j * W;
                out[i + j]    = u32((load_u64(b + (bit >> 3)) >> (bit & 7)) & mask);
            }
        }
        for (; i < n; i++) {
            const size_t bit = i * size_t(W);
            out[i]           = u32((load_u64(src + (bit >> 3)) >> (bit & 7)) & mask);
        }
    }
}
using Unpack32Fn = void (*)(const u8 *, size_t, u32 *);
template <size_t... Ws>
static constexpr std::array<Unpack32Fn, sizeof...(Ws)> make_unpack32_table(std::index_sequence<Ws...>) {
    return {&unpack32_w<int(Ws)>...};
}
static constexpr auto unpack32_table = make_unpack32_table(std::make_index_sequence<33>{});
}  // namespace detail
#endif

/// Unpack n values of width w; returns the number of bytes consumed (packed_bytes(n, w)).
template <class U>
SMESH_HOST_DEVICE SMESH_INLINE size_t bp_unpack(const u8 *SMESH_RESTRICT src, size_t n, int w,
                                                U *SMESH_RESTRICT out) {
#ifndef __CUDA_ARCH__
    if constexpr (sizeof(U) == 4) {
        detail::unpack32_table[size_t(w)](src, n, out);
        return packed_bytes(n, w);
    }
#endif
    bp_unpack_generic<U>(src, n, w, out);
    return packed_bytes(n, w);
}

// ------------------------------------------------------------------------------ byte sink
/// Append-only byte buffer for the encoders (host only).
class ByteSink {
public:
    std::vector<u8> bytes;

    size_t size() const { return bytes.size(); }

    template <class T>
    void put(T v) {
        const size_t p = bytes.size();
        bytes.resize(p + sizeof(T));
        memcpy(bytes.data() + p, &v, sizeof(T));
    }

    void append(const void *src, size_t n) {
        if (n == 0) return;
        const size_t p = bytes.size();
        bytes.resize(p + n);
        memcpy(bytes.data() + p, src, n);
    }

    /// Append n values of width w, bit-packed.
    template <class U>
    void pack(const U *v, size_t n, int w) {
        if (w == 0 || n == 0) return;
        const size_t p  = bytes.size();
        const size_t nb = packed_bytes(n, w);
        bytes.resize(p + nb + kSlack, 0);
        bp_pack<U>(bytes.data() + p, v, n, w);
        bytes.resize(p + nb);
    }
};

/// Read a value and advance (no bounds check: callers check before).
template <class T>
SMESH_HOST_DEVICE SMESH_INLINE T take(const u8 *&p) {
    T v;
    memcpy(&v, p, sizeof(T));
    p += sizeof(T);
    return v;
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CODEC_BITS_HPP
