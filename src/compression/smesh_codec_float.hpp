#ifndef SMESH_CODEC_FLOAT_HPP
#define SMESH_CODEC_FLOAT_HPP

// Lossless float encodings. The default, Codec::Packed, is enc_packed on ordered words
// (smesh_codec_int.hpp), bit-exact for every pattern. Codec::Shuffle regroups the bytes of n
// values by significance (all byte 0, then all byte 1, ...), which a following zstd / lz4 stage
// compresses well for smooth fields; on its own it does not shrink anything.

#include "smesh_compression_types.hpp"

namespace smesh {
namespace compression {

/// dst[k * n + i] = byte k of value i (value size `vb`).
static inline void byte_shuffle(const u8 *SMESH_RESTRICT src, size_t n, size_t vb, u8 *SMESH_RESTRICT dst) {
    for (size_t i = 0; i < n; i++)
        for (size_t k = 0; k < vb; k++) dst[k * n + i] = src[i * vb + k];
}

SMESH_HOST_DEVICE SMESH_INLINE void byte_unshuffle(const u8 *SMESH_RESTRICT src, size_t n, size_t vb,
                                                   u8 *SMESH_RESTRICT dst) {
    for (size_t k = 0; k < vb; k++)
        for (size_t i = 0; i < n; i++) dst[i * vb + k] = src[k * n + i];
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CODEC_FLOAT_HPP
