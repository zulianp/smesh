#include "smesh_byte_codec.hpp"

#include "smesh_config.hpp"

#include <climits>
#include <cstdio>

#ifdef SMESH_ENABLE_ZSTD
#include <zstd.h>
#endif
#ifdef SMESH_ENABLE_LZ4
#include <lz4.h>
#endif

namespace smesh {
namespace compression {

bool byte_stage_available(ByteStage stage) {
    switch (stage) {
        case ByteStage::None: return true;
        case ByteStage::Zstd:
#ifdef SMESH_ENABLE_ZSTD
            return true;
#else
            return false;
#endif
        case ByteStage::Lz4:
#ifdef SMESH_ENABLE_LZ4
            return true;
#else
            return false;
#endif
    }
    return false;
}

struct ByteEncoder::Impl {
#ifdef SMESH_ENABLE_ZSTD
    ZSTD_CCtx *cctx{nullptr};
    ~Impl() { ZSTD_freeCCtx(cctx); }
#endif
};

ByteEncoder::ByteEncoder() : impl_(std::make_unique<Impl>()) {}
ByteEncoder::~ByteEncoder() = default;

int ByteEncoder::compress(ByteStage stage, int level, const u8 *src, size_t n, std::vector<u8> &out) {
    (void)level;
    switch (stage) {
        case ByteStage::None:
            out.insert(out.end(), src, src + n);
            return SMESH_SUCCESS;
        case ByteStage::Zstd: {
#ifdef SMESH_ENABLE_ZSTD
            if (!impl_->cctx) impl_->cctx = ZSTD_createCCtx();
            const size_t p = out.size();
            out.resize(p + ZSTD_compressBound(n));
            const size_t r = ZSTD_compressCCtx(impl_->cctx, out.data() + p, out.size() - p, src, n, level);
            if (ZSTD_isError(r)) {
                fprintf(stderr, "zstd: %s\n", ZSTD_getErrorName(r));
                out.resize(p);
                return SMESH_FAILURE;
            }
            out.resize(p + r);
            return SMESH_SUCCESS;
#else
            break;
#endif
        }
        case ByteStage::Lz4: {
#ifdef SMESH_ENABLE_LZ4
            if (n > size_t(LZ4_MAX_INPUT_SIZE)) return SMESH_FAILURE;
            const size_t p     = out.size();
            const int    bound = LZ4_compressBound(int(n));
            out.resize(p + size_t(bound));
            const int r = LZ4_compress_default((const char *)src, (char *)out.data() + p, int(n), bound);
            if (r <= 0) {
                out.resize(p);
                return SMESH_FAILURE;
            }
            out.resize(p + size_t(r));
            return SMESH_SUCCESS;
#else
            break;
#endif
        }
    }
    fprintf(stderr, "compression: byte stage %s is not available in this build\n", to_string(stage));
    return SMESH_FAILURE;
}

struct ByteDecoder::Impl {
#ifdef SMESH_ENABLE_ZSTD
    ZSTD_DCtx *dctx{nullptr};
    ~Impl() { ZSTD_freeDCtx(dctx); }
#endif
};

ByteDecoder::ByteDecoder() : impl_(std::make_unique<Impl>()) {}
ByteDecoder::~ByteDecoder() = default;

int ByteDecoder::decompress(ByteStage stage, const u8 *src, size_t n, u8 *dst, size_t dst_n) {
    switch (stage) {
        case ByteStage::None:
            if (n != dst_n) return SMESH_FAILURE;
            memcpy(dst, src, n);
            return SMESH_SUCCESS;
        case ByteStage::Zstd: {
#ifdef SMESH_ENABLE_ZSTD
            if (!impl_->dctx) impl_->dctx = ZSTD_createDCtx();
            const size_t r = ZSTD_decompressDCtx(impl_->dctx, dst, dst_n, src, n);
            return (ZSTD_isError(r) || r != dst_n) ? SMESH_FAILURE : SMESH_SUCCESS;
#else
            break;
#endif
        }
        case ByteStage::Lz4: {
#ifdef SMESH_ENABLE_LZ4
            if (n > size_t(INT_MAX) || dst_n > size_t(INT_MAX)) return SMESH_FAILURE;
            const int r = LZ4_decompress_safe((const char *)src, (char *)dst, int(n), int(dst_n));
            return (r < 0 || size_t(r) != dst_n) ? SMESH_FAILURE : SMESH_SUCCESS;
#else
            break;
#endif
        }
    }
    fprintf(stderr, "compression: byte stage %s is not available in this build\n", to_string(stage));
    return SMESH_FAILURE;
}

}  // namespace compression
}  // namespace smesh
