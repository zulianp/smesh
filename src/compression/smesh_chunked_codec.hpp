#ifndef SMESH_CHUNKED_CODEC_HPP
#define SMESH_CHUNKED_CODEC_HPP

// Whole arrays as a sequence of independent chunks of `chunk_values` values (the last may be
// shorter). Encoding has three steps, so that serial and MPI writers share them and produce the
// same bytes:
//
//   local_stats -> (merge across ranks) -> plan_array -> encode_chunks
//
// The plan fixes everything that is global (codec, grid origin and step, ...); a chunk's bytes
// then depend only on its own values and the plan. Chunks are encoded in parallel (OpenMP),
// each into its own buffer, and joined in chunk order: the output does not depend on the number
// of threads or ranks.
//
// Per chunk, the encoder keeps the smallest of: the plan's codec (lossy chunks with too many
// escapes fall back to lossless), the lossless codec, and the raw bytes; then the byte stage if
// it saves at least 5%. Every chunk records what it used, its sizes and an XXH64 checksum.

#include "smesh_codec_quant.hpp"
#include "smesh_compression_options.hpp"

#include <vector>

namespace smesh {
namespace compression {

/// Statistics of (part of) an array, reducible across ranks with merge() (MPI: MIN / MAX / SUM).
struct ArrayStats {
    double vmin{INFINITY};     // finite values only
    double vmax{-INFINITY};
    double amax{0};            // max |x|
    double amin_nz{INFINITY};  // min |x| over nonzero finite values
    double bmin{INFINITY};     // NodeHeight: min positive finite bound
    u64    n_nonfinite{0};

    void merge(const ArrayStats &o);
};

/// Parameters fixed for a whole array.
struct ArrayPlan {
    DType       dtype{DType::F64};
    u64         n_values{0};
    u64         chunk_values{u64(1) << 16};
    u64         n_chunks{0};
    BoundMode   requested{BoundMode::Lossless};  // what the options asked for
    Codec       codec{Codec::Packed};            // what chunks try first
    ByteStage   byte_stage{ByteStage::None};
    int         zstd_level{3};
    bool        checksums{true};
    bool        verify_write{false};
    double      max_escape_fraction{0.5};
    QuantParams quant;  // lossy codecs only

    u64 chunk_begin(u64 c) const { return c * chunk_values; }
    u64 chunk_size(u64 c) const {
        const u64 b = c * chunk_values;
        return n_values - b < chunk_values ? n_values - b : chunk_values;
    }
};

/// One chunk of an encoded array (one index entry of an .smz file).
struct ChunkEntry {
    u64   offset{0};        // from the start of the payload
    u64   hash{0};          // XXH64 of the stored bytes (0 without checksums)
    u32   stored_bytes{0};  // bytes in the file
    u32   stage1_bytes{0};  // bytes before the byte stage (= stored_bytes without one)
    u32   n_escaped{0};
    Codec codec{Codec::Stored};
    ByteStage byte_stage{ByteStage::None};
    u16   flags{0};
};

struct EncodedChunks {
    std::vector<u8>         payload;
    std::vector<ChunkEntry> index;  // offsets relative to payload.data()
    QuantStats              stats;
    u64                     n_lossy_fallback{0};  // lossy chunks written losslessly
};

/// Statistics of n values; `bound` (NodeHeight: per-value absolute bound) may be null.
template <class T>
ArrayStats local_stats(const T *x, ptrdiff_t n, const double *bound);

/// Fix the parameters for an array of n_values with global statistics `g`. Returns
/// SMESH_FAILURE (with a message) for impossible requests, e.g. a lossy bound on integers.
int plan_array(const CompressionOptions &opts, DType dtype, u64 n_values, const ArrayStats &g, ArrayPlan *plan);

/// Encode values [first_value, first_value + n) of the array; first_value must be a chunk
/// boundary and n must end on one (or at the end of the array). `bound`: per-value absolute
/// bounds for NodeHeight (aligned with x), else null.
template <class T>
int encode_chunks(const ArrayPlan &plan, const T *x, u64 first_value, u64 n, const double *bound, EncodedChunks *out);

/// Decode values [begin, end) of an array stored with element type FT into out (converted to T).
/// `payload` holds the stored bytes of chunks [c0, c1) contiguously, starting at the file offset
/// index[c0].offset, followed by kSlack readable bytes; `index` is the array's full index.
/// Returns SMESH_FAILURE (with a message naming the chunk) on corrupt data.
template <class FT, class T>
int decode_chunks(const ArrayPlan &plan, const u8 *payload, const ChunkEntry *index, u64 c0, u64 c1, u64 begin,
                  u64 end, T *out, bool verify_checksums);

/// The chunks [c0, c1) that cover values [begin, end).
void chunks_covering(const ArrayPlan &plan, u64 begin, u64 end, u64 *c0, u64 *c1);

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_CHUNKED_CODEC_HPP
