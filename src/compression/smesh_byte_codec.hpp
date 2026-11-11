#ifndef SMESH_BYTE_CODEC_HPP
#define SMESH_BYTE_CODEC_HPP

// Optional general-purpose stage on encoded chunks: zstd (SMESH_ENABLE_ZSTD) or lz4
// (SMESH_ENABLE_LZ4). Contexts are per thread: create one ByteEncoder / ByteDecoder per thread
// and reuse it. The library headers stay private to smesh_byte_codec.cpp.

#include "smesh_compression_types.hpp"

#include <memory>
#include <vector>

namespace smesh {
namespace compression {

/// Whether this build can encode and decode the stage.
bool byte_stage_available(ByteStage stage);

class ByteEncoder {
public:
    ByteEncoder();
    ~ByteEncoder();
    ByteEncoder(const ByteEncoder &)            = delete;
    ByteEncoder &operator=(const ByteEncoder &) = delete;

    /// Compress src[0, n) and append the result to `out`. Returns SMESH_SUCCESS, or
    /// SMESH_FAILURE if the stage is unavailable or the library reports an error.
    int compress(ByteStage stage, int level, const u8 *src, size_t n, std::vector<u8> &out);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class ByteDecoder {
public:
    ByteDecoder();
    ~ByteDecoder();
    ByteDecoder(const ByteDecoder &)            = delete;
    ByteDecoder &operator=(const ByteDecoder &) = delete;

    /// Decompress src[0, n) into dst, which must receive exactly dst_n bytes.
    int decompress(ByteStage stage, const u8 *src, size_t n, u8 *dst, size_t dst_n);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_BYTE_CODEC_HPP
