#ifndef SMESH_SMZ_FORMAT_HPP
#define SMESH_SMZ_FORMAT_HPP

// The .smz file: one compressed array, self-describing. Little-endian throughout.
//
//   [header: 128 B][index: n_chunks x 32 B][payload: the chunks' stored bytes, in order]
//
// Header (offset, type, field):
//     0  u8[8] magic "\x89SMZ\r\n\x1a\n"
//     8  u16   version_major (1)          10  u16  version_minor (0)
//    12  u16   header_bytes (128)         14  u16  index_entry_bytes (32)
//    16  u8    dtype (DType)              17  u8   bound mode requested (BoundMode)
//    18  u8    codec planned (Codec)      19  u8   byte stage planned (ByteStage)
//    20  u8    checksums (0 none, 1 XXH64) 21 u8   log2 of the block size (7)
//    22  u8    q stored as u64 words      23  u8   reserved (0)
//    24  u64   n_values                   32  u64  chunk_values
//    40  u64   n_chunks                   48  u64  payload_bytes
//    56  f64   eb (as requested)          64  f64  eb_abs (QuantFixed)
//    72  f64   grid origin                80  i32  L (QuantFixed)   84 i32 L_min (QuantVar)
//    88  f64   min of the finite values   96  f64  max of the finite values
//   104  f64   achieved max |x' - x|     112  f32  achieved max |x' - x| / bound
//   116  u32   flags (bit 0: verified on write)
//   120  u64   XXH64 of bytes [0, 120)
//
// Index entry (32 B): u64 offset in the payload, u64 XXH64 of the stored bytes, u32 stored
// bytes, u32 bytes before the byte stage, u32 escaped values, u8 codec, u8 byte stage, u16 flags.
//
// Readers check the header hash and version, that the chunks tile the payload exactly, and
// that the file is header + index + payload long; every chunk's hash is checked on decode.

#include "smesh_chunked_codec.hpp"

#include <string>

namespace smesh {
namespace compression {

static constexpr size_t kSmzHeaderBytes     = 128;
static constexpr size_t kSmzIndexEntryBytes = 32;
static constexpr u16    kSmzVersionMajor    = 1;
static constexpr u16    kSmzVersionMinor    = 0;

/// Everything in an .smz header.
struct SmzInfo {
    u16       version_major{kSmzVersionMajor};
    u16       version_minor{kSmzVersionMinor};
    ArrayPlan plan;
    u64       payload_bytes{0};
    double    vmin{0}, vmax{0};
    double    max_abs_err{0};
    double    max_ratio{0};
    u32       flags{0};

    u64 index_offset() const { return kSmzHeaderBytes; }
    u64 payload_offset() const { return kSmzHeaderBytes + kSmzIndexEntryBytes * plan.n_chunks; }
    u64 file_bytes() const { return payload_offset() + payload_bytes; }
};

bool smz_host_is_little_endian();

void smz_encode_header(const SmzInfo &info, u8 out[kSmzHeaderBytes]);
/// SMESH_FAILURE with a reason in *why for anything that is not a valid v1 header.
int smz_decode_header(const u8 in[kSmzHeaderBytes], SmzInfo *info, std::string *why);

void       smz_encode_index_entry(const ChunkEntry &e, u8 out[kSmzIndexEntryBytes]);
ChunkEntry smz_decode_index_entry(const u8 in[kSmzIndexEntryBytes]);

/// Index entries [c0, c0 + n) must tile the payload in order: offset[c0] = `expected_offset`,
/// offset[c+1] = offset[c] + stored[c], codes valid. With c0 + n = n_chunks the last one must
/// end at payload_bytes.
int smz_validate_index(const SmzInfo &info, const ChunkEntry *entries, u64 c0, u64 n, u64 expected_offset,
                       std::string *why);

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_SMZ_FORMAT_HPP
