#include "smesh_chunked_codec.impl.hpp"

#include "smesh_byte_codec.hpp"

namespace smesh {
namespace compression {

void ArrayStats::merge(const ArrayStats &o) {
    vmin    = o.vmin < vmin ? o.vmin : vmin;
    vmax    = o.vmax > vmax ? o.vmax : vmax;
    amax    = o.amax > amax ? o.amax : amax;
    amin_nz = o.amin_nz < amin_nz ? o.amin_nz : amin_nz;
    bmin    = o.bmin < bmin ? o.bmin : bmin;
    n_nonfinite += o.n_nonfinite;
}

void chunks_covering(const ArrayPlan &plan, u64 begin, u64 end, u64 *c0, u64 *c1) {
    if (end <= begin) {
        *c0 = *c1 = begin / plan.chunk_values;
        return;
    }
    *c0 = begin / plan.chunk_values;
    *c1 = (end + plan.chunk_values - 1) / plan.chunk_values;
}

int plan_array(const CompressionOptions &opts, DType dtype, u64 n_values, const ArrayStats &g, ArrayPlan *plan) {
    ArrayPlan p;
    p.dtype               = dtype;
    p.n_values            = n_values;
    p.chunk_values        = u64(opts.chunk_values > 0 ? opts.chunk_values : (i64(1) << 16));
    p.n_chunks            = (n_values + p.chunk_values - 1) / p.chunk_values;
    p.requested           = opts.enabled ? opts.mode : BoundMode::Lossless;
    p.byte_stage          = opts.byte_stage;
    p.zstd_level          = opts.zstd_level;
    p.checksums           = opts.checksums;
    p.verify_write        = opts.verify_write;
    p.max_escape_fraction = opts.max_escape_fraction;
    p.quant.mode          = p.requested;
    p.quant.eb            = opts.eb;

    if (!byte_stage_available(p.byte_stage)) {
        fprintf(stderr, "compression: byte stage %s is not available in this build\n", to_string(p.byte_stage));
        return SMESH_FAILURE;
    }
    if (p.requested != BoundMode::Lossless && !dtype_is_float(dtype)) {
        fprintf(stderr, "compression: a lossy bound (%s) needs floating-point values, not %s\n",
                to_string(p.requested), to_string(dtype));
        return SMESH_FAILURE;
    }

    const bool finite = std::isfinite(g.vmin) && std::isfinite(g.vmax);
    p.codec           = opts.shuffle && dtype_is_float(dtype) ? Codec::Shuffle : Codec::Packed;
    switch (p.requested) {
        case BoundMode::Lossless: break;
        case BoundMode::Abs:
            if (finite && plan_quant_fixed(opts.eb, g.vmin, g.vmax, p.quant)) p.codec = Codec::QuantFixed;
            break;
        case BoundMode::RelRange:
            if (finite && g.vmax > g.vmin && plan_quant_fixed(opts.eb * (g.vmax - g.vmin), g.vmin, g.vmax, p.quant))
                p.codec = Codec::QuantFixed;
            break;
        case BoundMode::PointwiseRel:
            if (finite && g.amax > 0) {
                plan_quant_var(0.0, g.amax, std::isfinite(g.amin_nz) ? opts.eb * g.amin_nz : 0.0, p.quant);
                p.codec = Codec::QuantVar;
            }
            break;
        case BoundMode::NodeHeight:
            if (finite) {
                plan_quant_var(g.vmin, g.vmax - g.vmin, std::isfinite(g.bmin) ? g.bmin : 0.0, p.quant);
                p.codec = Codec::QuantVar;
            }
            break;
    }
    *plan = p;
    return SMESH_SUCCESS;
}

#define SMESH_CHUNKED_INSTANTIATE_T(T)                                                                                 \
    template ArrayStats local_stats<T>(const T *, ptrdiff_t, const double *);                                         \
    template int        encode_chunks<T>(const ArrayPlan &, const T *, u64, u64, const double *, EncodedChunks *);

#define SMESH_CHUNKED_INSTANTIATE_FT(FT, T)                                                                            \
    template int decode_chunks<FT, T>(const ArrayPlan &, const u8 *, const ChunkEntry *, u64, u64, u64, u64, T *, bool);

#define SMESH_CHUNKED_INSTANTIATE_ALL_FT(T)  \
    SMESH_CHUNKED_INSTANTIATE_FT(i8, T)      \
    SMESH_CHUNKED_INSTANTIATE_FT(u8, T)      \
    SMESH_CHUNKED_INSTANTIATE_FT(i16, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(u16, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(i32, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(u32, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(i64, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(u64, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(f16, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(f32, T)     \
    SMESH_CHUNKED_INSTANTIATE_FT(f64, T)

#define SMESH_CHUNKED_INSTANTIATE(T) SMESH_CHUNKED_INSTANTIATE_T(T) SMESH_CHUNKED_INSTANTIATE_ALL_FT(T)

SMESH_CHUNKED_INSTANTIATE(i8)
SMESH_CHUNKED_INSTANTIATE(u8)
SMESH_CHUNKED_INSTANTIATE(i16)
SMESH_CHUNKED_INSTANTIATE(u16)
SMESH_CHUNKED_INSTANTIATE(i32)
SMESH_CHUNKED_INSTANTIATE(u32)
SMESH_CHUNKED_INSTANTIATE(i64)
SMESH_CHUNKED_INSTANTIATE(u64)
SMESH_CHUNKED_INSTANTIATE(f16)
SMESH_CHUNKED_INSTANTIATE(f32)
SMESH_CHUNKED_INSTANTIATE(f64)

}  // namespace compression
}  // namespace smesh
