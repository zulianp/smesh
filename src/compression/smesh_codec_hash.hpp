#ifndef SMESH_CODEC_HASH_HPP
#define SMESH_CODEC_HASH_HPP

// XXH64 (the public xxHash algorithm, seed 0 by default): checksums of .smz headers and chunks.
// Own implementation, bit-compatible with the reference, so files can be checked by any xxhash
// tool (e.g. Python's xxhash.xxh64).

#include "smesh_compression_types.hpp"

namespace smesh {
namespace compression {

namespace detail {
static constexpr u64 kP1 = 0x9E3779B185EBCA87ull;
static constexpr u64 kP2 = 0xC2B2AE3D27D4EB4Full;
static constexpr u64 kP3 = 0x165667B19E3779F9ull;
static constexpr u64 kP4 = 0x85EBCA77C2B2AE63ull;
static constexpr u64 kP5 = 0x27D4EB2F165667C5ull;

SMESH_HOST_DEVICE SMESH_INLINE u64 rotl64(u64 x, int r) { return (x << r) | (x >> (64 - r)); }
SMESH_HOST_DEVICE SMESH_INLINE u64 xxh_round(u64 acc, u64 in) {
    acc += in * kP2;
    acc = rotl64(acc, 31);
    return acc * kP1;
}
SMESH_HOST_DEVICE SMESH_INLINE u64 xxh_merge(u64 acc, u64 v) {
    acc ^= xxh_round(0, v);
    return acc * kP1 + kP4;
}
SMESH_HOST_DEVICE SMESH_INLINE u64 read64(const u8 *p) {
    u64 v;
    memcpy(&v, p, 8);
    return v;  // little-endian hosts only (checked when files are opened)
}
SMESH_HOST_DEVICE SMESH_INLINE u32 read32(const u8 *p) {
    u32 v;
    memcpy(&v, p, 4);
    return v;
}
}  // namespace detail

SMESH_HOST_DEVICE SMESH_INLINE u64 xxh64(const void *data, size_t len, u64 seed = 0) {
    using namespace detail;
    const u8 *p   = static_cast<const u8 *>(data);
    const u8 *end = p + len;
    u64       h;
    if (len >= 32) {
        u64 v1 = seed + kP1 + kP2, v2 = seed + kP2, v3 = seed, v4 = seed - kP1;
        const u8 *limit = end - 32;
        do {
            v1 = xxh_round(v1, read64(p));
            v2 = xxh_round(v2, read64(p + 8));
            v3 = xxh_round(v3, read64(p + 16));
            v4 = xxh_round(v4, read64(p + 24));
            p += 32;
        } while (p <= limit);
        h = rotl64(v1, 1) + rotl64(v2, 7) + rotl64(v3, 12) + rotl64(v4, 18);
        h = xxh_merge(h, v1);
        h = xxh_merge(h, v2);
        h = xxh_merge(h, v3);
        h = xxh_merge(h, v4);
    } else {
        h = seed + kP5;
    }
    h += u64(len);
    while (p + 8 <= end) {
        h ^= xxh_round(0, read64(p));
        h = rotl64(h, 27) * kP1 + kP4;
        p += 8;
    }
    if (p + 4 <= end) {
        h ^= u64(read32(p)) * kP1;
        h = rotl64(h, 23) * kP2 + kP3;
        p += 4;
    }
    while (p < end) {
        h ^= u64(*p) * kP5;
        h = rotl64(h, 11) * kP1;
        p++;
    }
    h ^= h >> 33;
    h *= kP2;
    h ^= h >> 29;
    h *= kP3;
    h ^= h >> 32;
    return h;
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CODEC_HASH_HPP
