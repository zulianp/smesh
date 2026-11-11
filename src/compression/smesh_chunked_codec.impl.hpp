#ifndef SMESH_CHUNKED_CODEC_IMPL_HPP
#define SMESH_CHUNKED_CODEC_IMPL_HPP

#include "smesh_byte_codec.hpp"
#include "smesh_chunked_codec.hpp"
#include "smesh_codec_float.hpp"
#include "smesh_codec_hash.hpp"
#include "smesh_codec_int.hpp"
#include "smesh_codec_quant.hpp"

#include <cmath>
#include <cstdio>

namespace smesh {
namespace compression {

template <class T>
ArrayStats local_stats(const T *x, ptrdiff_t n, const double *bound) {
    double vmin = INFINITY, vmax = -INFINITY, amax = 0, amin_nz = INFINITY, bmin = INFINITY;
    u64    nnf  = 0;
#pragma omp parallel for reduction(min : vmin, amin_nz, bmin) reduction(max : vmax, amax) reduction(+ : nnf)
    for (ptrdiff_t i = 0; i < n; i++) {
        const double v = double(x[i]);
        if (!std::isfinite(v)) {
            nnf++;
            continue;
        }
        vmin           = v < vmin ? v : vmin;
        vmax           = v > vmax ? v : vmax;
        const double a = std::fabs(v);
        amax           = a > amax ? a : amax;
        if (a > 0) amin_nz = a < amin_nz ? a : amin_nz;
        if (bound) {
            const double b = bound[i];
            if (b > 0 && std::isfinite(b)) bmin = b < bmin ? b : bmin;
        }
    }
    ArrayStats s;
    s.vmin        = vmin;
    s.vmax        = vmax;
    s.amax        = amax;
    s.amin_nz     = amin_nz;
    s.bmin        = bmin;
    s.n_nonfinite = nnf;
    return s;
}

namespace detail {

/// Per-thread buffers for decoding one chunk.
template <class FT>
struct ChunkScratch {
    std::vector<u8>  stage;  // byte stage output
    std::vector<FT>  values;
    std::vector<u32> a, b;
    std::vector<u64> w;
    ByteDecoder      bdec;

    void reserve(size_t m) {
        if (values.size() < m) {
            values.resize(m);
            a.resize(m);
            b.resize(m);
            w.resize(m);
        }
    }
};

/// Decode one chunk of m values (stored bytes `src`) into out. Returns SMESH_SUCCESS or
/// SMESH_FAILURE on malformed data.
template <class FT>
int decode_one(const ArrayPlan &plan, const ChunkEntry &e, const u8 *src, size_t m, FT *out, ChunkScratch<FT> &sc) {
    const u8 *p   = src;
    size_t    len = e.stored_bytes;
    if (e.byte_stage != ByteStage::None) {
        if (sc.stage.size() < size_t(e.stage1_bytes) + kSlack) sc.stage.resize(size_t(e.stage1_bytes) + kSlack);
        if (sc.bdec.decompress(e.byte_stage, src, len, sc.stage.data(), e.stage1_bytes) != SMESH_SUCCESS)
            return SMESH_FAILURE;
        p   = sc.stage.data();
        len = e.stage1_bytes;
    } else if (e.stage1_bytes != e.stored_bytes) {
        return SMESH_FAILURE;
    }
    const u8 *end = p + len;
    sc.reserve(m);
    switch (e.codec) {
        case Codec::Stored:
            if (len != m * sizeof(FT)) return SMESH_FAILURE;
            memcpy(out, p, len);
            return SMESH_SUCCESS;
        case Codec::Packed: return dec_packed<FT, true>(p, end, m, out) == end ? SMESH_SUCCESS : SMESH_FAILURE;
        case Codec::Shuffle:
            if (len != m * sizeof(FT)) return SMESH_FAILURE;
            byte_unshuffle(p, m, sizeof(FT), reinterpret_cast<u8 *>(out));
            return SMESH_SUCCESS;
        case Codec::QuantFixed:
            if constexpr (is_float_v<FT>) {
                return dec_quant_fixed<FT, true>(p, end, m, plan.quant, out, sc.a.data(), sc.w.data(), sc.b.data()) == end
                               ? SMESH_SUCCESS
                               : SMESH_FAILURE;
            }
            return SMESH_FAILURE;
        case Codec::QuantVar:
            if constexpr (is_float_v<FT>) {
                return dec_quant_var<FT, true>(p, end, m, plan.quant, out, sc.a.data(), sc.b.data()) == end
                               ? SMESH_SUCCESS
                               : SMESH_FAILURE;
            }
            return SMESH_FAILURE;
    }
    return SMESH_FAILURE;
}

template <class T>
static inline bool same_bits(const T &a, const T &b) {
    return memcmp(&a, &b, sizeof(T)) == 0;
}

/// The bound a decoded value of a chunk must meet (verification after encoding).
template <class T>
static inline double chunk_bound(const ArrayPlan &plan, Codec used, const T &x, const double *bound, size_t i) {
    if (used == Codec::QuantFixed) return plan.quant.eb_abs;
    if (used == Codec::QuantVar) {
        return plan.requested == BoundMode::PointwiseRel ? plan.quant.eb * std::fabs(double(x)) : bound[i];
    }
    return 0;
}

/// Encode one chunk of m values (bound aligned with x, may be null). Fills e (all but offset).
template <class T>
int encode_one(const ArrayPlan &plan, const T *x, size_t m, const double *bound, ByteEncoder &benc,
               std::vector<u8> &stored, ChunkEntry &e, QuantStats &cs, u64 &fallback) {
    ByteSink   s;
    Codec      used = plan.codec;
    QuantStats qs;
    if constexpr (is_float_v<T>) {
        if (used == Codec::QuantFixed) {
            enc_quant_fixed<T>(s, x, m, plan.quant, qs);
        } else if (used == Codec::QuantVar) {
            if (plan.requested == BoundMode::PointwiseRel) {
                const double eb = plan.quant.eb;
                enc_quant_var<T>(
                        s, x, m, plan.quant, [](size_t i) { return i; },
                        [x, eb](size_t i) { return eb * std::fabs(double(x[i])); }, qs);
            } else {
                enc_quant_var<T>(s, x, m, plan.quant, bound, qs);
            }
        }
        if ((used == Codec::QuantFixed || used == Codec::QuantVar) &&
            double(qs.n_escaped) > plan.max_escape_fraction * double(m)) {
            s.bytes.clear();
            qs   = QuantStats();
            used = Codec::Packed;
            fallback++;
        }
    } else if (used == Codec::QuantFixed || used == Codec::QuantVar) {
        return SMESH_FAILURE;  // plan_array never does this
    }
    if (used == Codec::Shuffle) {
        s.bytes.resize(m * sizeof(T));
        byte_shuffle(reinterpret_cast<const u8 *>(x), m, sizeof(T), s.bytes.data());
    } else if (used == Codec::Packed) {
        enc_packed<T>(s, x, m);
    }
    if (s.size() >= m * sizeof(T) && used != Codec::Shuffle) {  // incompressible: the raw bytes
        s.bytes.assign(reinterpret_cast<const u8 *>(x), reinterpret_cast<const u8 *>(x) + m * sizeof(T));
        qs   = QuantStats();
        used = Codec::Stored;
    }

    stored.clear();
    ByteStage stage = ByteStage::None;
    if (plan.byte_stage != ByteStage::None) {
        if (benc.compress(plan.byte_stage, plan.zstd_level, s.bytes.data(), s.size(), stored) != SMESH_SUCCESS)
            return SMESH_FAILURE;
        if (double(stored.size()) <= 0.95 * double(s.size())) {
            stage = plan.byte_stage;
        } else {
            stored.clear();
        }
    }
    if (stage == ByteStage::None) stored.swap(s.bytes);

    e.stored_bytes = u32(stored.size());
    e.stage1_bytes = u32(stage == ByteStage::None ? stored.size() : s.size());
    e.n_escaped    = u32(qs.n_escaped);
    e.codec        = used;
    e.byte_stage   = stage;
    e.hash         = plan.checksums ? xxh64(stored.data(), stored.size()) : 0;
    cs             = qs;

    if (plan.verify_write) {
        ChunkScratch<T> sc;
        std::vector<T>  y(m);
        std::vector<u8> buf(stored);
        buf.resize(stored.size() + kSlack);
        if (decode_one<T>(plan, e, buf.data(), m, y.data(), sc) != SMESH_SUCCESS) {
            fprintf(stderr, "compression: verify: chunk does not decode\n");
            return SMESH_FAILURE;
        }
        for (size_t i = 0; i < m; i++) {
            if (used == Codec::QuantFixed || used == Codec::QuantVar) {
                if (std::isfinite(double(x[i]))) {
                    if (!(std::fabs(double(y[i]) - double(x[i])) <= chunk_bound(plan, used, x[i], bound, i))) {
                        fprintf(stderr, "compression: verify: bound violated at value %zu of a chunk\n", i);
                        return SMESH_FAILURE;
                    }
                    continue;
                }
            }
            if (!same_bits(y[i], x[i])) {
                fprintf(stderr, "compression: verify: value %zu of a chunk differs\n", i);
                return SMESH_FAILURE;
            }
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace detail

template <class T>
int encode_chunks(const ArrayPlan &plan, const T *x, u64 first_value, u64 n, const double *bound, EncodedChunks *out) {
    const u64 cv = plan.chunk_values;
    if (first_value % cv != 0 || first_value + n > plan.n_values ||
        ((first_value + n) % cv != 0 && first_value + n != plan.n_values)) {
        fprintf(stderr, "compression: encode_chunks: range [%llu, %llu) is not chunk-aligned\n",
                (unsigned long long)first_value, (unsigned long long)(first_value + n));
        return SMESH_FAILURE;
    }
    if (plan.codec == Codec::QuantVar && plan.requested == BoundMode::NodeHeight && !bound && n > 0) {
        fprintf(stderr, "compression: bound mode h needs per-node bounds (element heights)\n");
        return SMESH_FAILURE;
    }
    const u64       c0 = first_value / cv;
    const ptrdiff_t nc = ptrdiff_t((n + cv - 1) / cv);

    const size_t                 n_chunks = static_cast<size_t>(nc);
    std::vector<std::vector<u8>> bufs(n_chunks);
    std::vector<ChunkEntry>      idx(n_chunks);
    std::vector<QuantStats>      stats(n_chunks);
    int                          err      = 0;
    u64                          fallback = 0;
#pragma omp parallel
    {
        ByteEncoder benc;
#pragma omp for schedule(dynamic, 1) reduction(max : err) reduction(+ : fallback)
        for (ptrdiff_t k = 0; k < nc; k++) {
            const u64 c   = c0 + u64(k);
            const u64 off = plan.chunk_begin(c) - first_value;
            const u64 m   = plan.chunk_size(c);
            if (detail::encode_one<T>(plan, x + off, size_t(m), bound ? bound + off : nullptr, benc, bufs[size_t(k)],
                                      idx[size_t(k)], stats[size_t(k)], fallback) != SMESH_SUCCESS) {
                err = 1;
            }
        }
    }
    if (err) return SMESH_FAILURE;

    // Join in chunk order.
    u64 total = 0;
    for (ptrdiff_t k = 0; k < nc; k++) {
        idx[size_t(k)].offset = total;
        total += idx[size_t(k)].stored_bytes;
    }
    out->payload.resize(size_t(total));
#pragma omp parallel for schedule(static)
    for (ptrdiff_t k = 0; k < nc; k++) {
        const auto &b = bufs[size_t(k)];
        if (!b.empty()) memcpy(out->payload.data() + idx[size_t(k)].offset, b.data(), b.size());
    }
    out->index = std::move(idx);
    out->stats = QuantStats();
    for (const auto &s : stats) out->stats.merge(s);
    out->n_lossy_fallback = fallback;
    return SMESH_SUCCESS;
}

template <class FT, class T>
int decode_chunks(const ArrayPlan &plan, const u8 *payload, const ChunkEntry *index, u64 c0, u64 c1, u64 begin,
                  u64 end, T *out, bool verify_checksums) {
    if (c1 <= c0) return SMESH_SUCCESS;
    const u64 base = index[c0].offset;
    int       err  = 0;
#pragma omp parallel
    {
        detail::ChunkScratch<FT> sc;
#pragma omp for schedule(dynamic, 1) reduction(max : err)
        for (ptrdiff_t k = 0; k < ptrdiff_t(c1 - c0); k++) {
            const u64         c   = c0 + u64(k);
            const ChunkEntry &e   = index[c];
            const u8         *src = payload + (e.offset - base);
            const u64         vb  = plan.chunk_begin(c);
            const u64         m   = plan.chunk_size(c);
            if (verify_checksums && plan.checksums && xxh64(src, e.stored_bytes) != e.hash) {
                fprintf(stderr, "compression: checksum mismatch in chunk %llu\n", (unsigned long long)c);
                err = 1;
                continue;
            }
            const u64 lo = begin > vb ? begin : vb;
            const u64 hi = end < vb + m ? end : vb + m;
            if (lo >= hi) continue;
            if constexpr (std::is_same_v<FT, T>) {
                if (lo == vb && hi == vb + m) {
                    if (detail::decode_one<FT>(plan, e, src, size_t(m), out + (vb - begin), sc) != SMESH_SUCCESS) {
                        fprintf(stderr, "compression: corrupt chunk %llu\n", (unsigned long long)c);
                        err = 1;
                    }
                    continue;
                }
            }
            sc.reserve(size_t(m));
            std::vector<FT> &tmp = sc.values;
            // decode_one uses sc.a / sc.b / sc.w as scratch, not sc.values
            if (detail::decode_one<FT>(plan, e, src, size_t(m), tmp.data(), sc) != SMESH_SUCCESS) {
                fprintf(stderr, "compression: corrupt chunk %llu\n", (unsigned long long)c);
                err = 1;
                continue;
            }
            for (u64 i = lo; i < hi; i++) out[i - begin] = T(tmp[size_t(i - vb)]);
        }
    }
    return err ? SMESH_FAILURE : SMESH_SUCCESS;
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CHUNKED_CODEC_IMPL_HPP
