#ifndef SMESH_COMPRESSION_TYPES_HPP
#define SMESH_COMPRESSION_TYPES_HPP

// Shared vocabulary of the compression module (src/compression): value types with stable
// on-disk codes, bound modes, codecs, byte stages, and the mapping of every value type to an
// order-preserving unsigned word, on which all integer codecs work.
//
// Layering (see docs/compression.md):
//   smesh_codec_*.hpp       header-only kernels: no MPI, no IO, no allocation in decoders,
//                           SMESH_HOST_DEVICE decoders (usable inside sfem's operator kernels)
//   smesh_byte_codec.*      optional zstd / lz4 stage
//   smesh_chunked_codec.*   whole arrays in self-contained chunks, OpenMP, deterministic
//   smesh_smz_*.*           the .smz file format and its serial IO; MPI IO in src/distributed

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <cstring>
#include <type_traits>

namespace smesh {
namespace compression {

/// Decoders may read up to this many bytes past the end of the encoded data (64-bit loads):
/// every buffer handed to a decoder must have them allocated (their content does not matter).
static constexpr size_t kSlack = 8;

/// Values per block of the integer codecs: each block has its own reference and bit width.
static constexpr size_t kBlock = 128;

/// Value types, with codes that are stored in files and must never change.
enum class DType : u8 { I8 = 1, U8 = 2, I16 = 3, U16 = 4, I32 = 5, U32 = 6, I64 = 7, U64 = 8, F16 = 9, F32 = 10, F64 = 11 };

/// Error bound of a lossy encoding (codes stored in files).
///   Lossless      bit-exact
///   Abs           |x' - x| <= eb
///   RelRange      |x' - x| <= eb * (max - min), over the finite values of the whole array
///   PointwiseRel  |x' - x| <= eb * |x|  (signs and zeros exact)
///   NodeHeight    |x' - x| <= eb * h_n, h_n the smallest height of the elements at node n
enum class BoundMode : u8 { Lossless = 0, Abs = 1, RelRange = 2, PointwiseRel = 3, NodeHeight = 4 };

/// How a chunk is encoded (codes stored in files).
///   Stored      the values' bytes, as they are
///   Packed      lossless: values mapped to ordered words, per block FOR or delta-FOR, bit-packed
///   Shuffle     lossless: byte planes (meant to be followed by zstd / lz4)
///   QuantFixed  lossy, one power-of-two step for the whole array (Abs, RelRange)
///   QuantVar    lossy, a power-of-two step per value (PointwiseRel, NodeHeight), with escapes
enum class Codec : u8 { Stored = 0, Packed = 1, Shuffle = 2, QuantFixed = 3, QuantVar = 4 };

/// General-purpose stage applied to an encoded chunk (codes stored in files).
enum class ByteStage : u8 { None = 0, Zstd = 1, Lz4 = 2 };

static inline const char *to_string(DType t) {
    switch (t) {
        case DType::I8: return "int8";
        case DType::U8: return "uint8";
        case DType::I16: return "int16";
        case DType::U16: return "uint16";
        case DType::I32: return "int32";
        case DType::U32: return "uint32";
        case DType::I64: return "int64";
        case DType::U64: return "uint64";
        case DType::F16: return "float16";
        case DType::F32: return "float32";
        case DType::F64: return "float64";
    }
    return "?";
}

static inline const char *to_string(BoundMode m) {
    switch (m) {
        case BoundMode::Lossless: return "lossless";
        case BoundMode::Abs: return "abs";
        case BoundMode::RelRange: return "range";
        case BoundMode::PointwiseRel: return "pwrel";
        case BoundMode::NodeHeight: return "h";
    }
    return "?";
}

static inline const char *to_string(Codec c) {
    switch (c) {
        case Codec::Stored: return "stored";
        case Codec::Packed: return "packed";
        case Codec::Shuffle: return "shuffle";
        case Codec::QuantFixed: return "quant";
        case Codec::QuantVar: return "quantvar";
    }
    return "?";
}

static inline const char *to_string(ByteStage s) {
    switch (s) {
        case ByteStage::None: return "none";
        case ByteStage::Zstd: return "zstd";
        case ByteStage::Lz4: return "lz4";
    }
    return "?";
}

static inline size_t dtype_bytes(DType t) {
    switch (t) {
        case DType::I8:
        case DType::U8: return 1;
        case DType::I16:
        case DType::U16:
        case DType::F16: return 2;
        case DType::I32:
        case DType::U32:
        case DType::F32: return 4;
        case DType::I64:
        case DType::U64:
        case DType::F64: return 8;
    }
    return 0;
}

static inline bool dtype_valid(u8 code) { return code >= 1 && code <= 11; }
static inline bool dtype_is_float(DType t) { return t == DType::F16 || t == DType::F32 || t == DType::F64; }

template <class T>
struct DTypeOf;
template <> struct DTypeOf<i8> { static constexpr DType value = DType::I8; };
template <> struct DTypeOf<u8> { static constexpr DType value = DType::U8; };
template <> struct DTypeOf<i16> { static constexpr DType value = DType::I16; };
template <> struct DTypeOf<u16> { static constexpr DType value = DType::U16; };
template <> struct DTypeOf<i32> { static constexpr DType value = DType::I32; };
template <> struct DTypeOf<u32> { static constexpr DType value = DType::U32; };
template <> struct DTypeOf<i64> { static constexpr DType value = DType::I64; };
template <> struct DTypeOf<u64> { static constexpr DType value = DType::U64; };
template <> struct DTypeOf<f16> { static constexpr DType value = DType::F16; };
template <> struct DTypeOf<f32> { static constexpr DType value = DType::F32; };
template <> struct DTypeOf<f64> { static constexpr DType value = DType::F64; };

template <class T>
static constexpr bool is_float_v = std::is_same_v<T, f16> || std::is_same_v<T, f32> || std::is_same_v<T, f64>;

// ------------------------------------------------------------------------------ ordered words
// Every value type maps to an unsigned word (u32 for types up to 4 bytes, u64 otherwise) such
// that the map is a bijection on the type's bit patterns and preserves order:
//   unsigned  v
//   signed    v with its sign bit flipped (two's complement -> offset binary)
//   float     bits with the sign bit set if non-negative, all bits flipped if negative
// so the integer codecs below see small ranges for clustered values, and floats round-trip
// bit-exactly (NaN payloads, -0, denormals).

template <class T>
using word_t = std::conditional_t<(sizeof(T) <= 4), u32, u64>;

template <size_t N>
struct uint_of_size;
template <> struct uint_of_size<1> { using type = u8; };
template <> struct uint_of_size<2> { using type = u16; };
template <> struct uint_of_size<4> { using type = u32; };
template <> struct uint_of_size<8> { using type = u64; };

template <class T>
SMESH_HOST_DEVICE SMESH_INLINE word_t<T> to_word(T v) {
    using B                = typename uint_of_size<sizeof(T)>::type;
    constexpr B sign       = B(B(1) << (8 * sizeof(T) - 1));
    B           b;
    memcpy(&b, &v, sizeof(T));
    if constexpr (is_float_v<T>) {
        b = (b & sign) ? B(~b) : B(b | sign);
    } else if constexpr (std::is_signed_v<T>) {
        b = B(b ^ sign);
    }
    return word_t<T>(b);
}

template <class T>
SMESH_HOST_DEVICE SMESH_INLINE T from_word(word_t<T> w) {
    using B          = typename uint_of_size<sizeof(T)>::type;
    constexpr B sign = B(B(1) << (8 * sizeof(T) - 1));
    B           b    = B(w);
    if constexpr (is_float_v<T>) {
        b = (b & sign) ? B(b & B(~sign)) : B(~b);
    } else if constexpr (std::is_signed_v<T>) {
        b = B(b ^ sign);
    }
    T v;
    memcpy(&v, &b, sizeof(T));
    return v;
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_COMPRESSION_TYPES_HPP
