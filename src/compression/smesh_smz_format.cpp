#include "smesh_smz_format.hpp"

#include "smesh_codec_hash.hpp"

#include <cstring>

namespace smesh {
namespace compression {

namespace {
const u8 kMagic[8] = {0x89, 'S', 'M', 'Z', '\r', '\n', 0x1a, '\n'};

template <class T>
void put_at(u8 *buf, size_t off, T v) {
    memcpy(buf + off, &v, sizeof(T));
}
template <class T>
T get_at(const u8 *buf, size_t off) {
    T v;
    memcpy(&v, buf + off, sizeof(T));
    return v;
}
bool codec_valid(u8 c) { return c <= u8(Codec::QuantVar); }
bool stage_valid(u8 s) { return s <= u8(ByteStage::Lz4); }
bool mode_valid(u8 m) { return m <= u8(BoundMode::NodeHeight); }
}  // namespace

bool smz_host_is_little_endian() {
    const u16 one = 1;
    u8        b;
    memcpy(&b, &one, 1);
    return b == 1;
}

void smz_encode_header(const SmzInfo &info, u8 out[kSmzHeaderBytes]) {
    memset(out, 0, kSmzHeaderBytes);
    const ArrayPlan &p = info.plan;
    memcpy(out, kMagic, 8);
    put_at<u16>(out, 8, info.version_major);
    put_at<u16>(out, 10, info.version_minor);
    put_at<u16>(out, 12, u16(kSmzHeaderBytes));
    put_at<u16>(out, 14, u16(kSmzIndexEntryBytes));
    put_at<u8>(out, 16, u8(p.dtype));
    put_at<u8>(out, 17, u8(p.requested));
    put_at<u8>(out, 18, u8(p.codec));
    put_at<u8>(out, 19, u8(p.byte_stage));
    put_at<u8>(out, 20, p.checksums ? 1 : 0);
    put_at<u8>(out, 21, 7);
    put_at<u8>(out, 22, p.quant.wide);
    put_at<u64>(out, 24, p.n_values);
    put_at<u64>(out, 32, p.chunk_values);
    put_at<u64>(out, 40, p.n_chunks);
    put_at<u64>(out, 48, info.payload_bytes);
    put_at<f64>(out, 56, p.quant.eb);
    put_at<f64>(out, 64, p.quant.eb_abs);
    put_at<f64>(out, 72, p.quant.origin);
    put_at<i32>(out, 80, p.quant.L);
    put_at<i32>(out, 84, p.quant.L_min);
    put_at<f64>(out, 88, info.vmin);
    put_at<f64>(out, 96, info.vmax);
    put_at<f64>(out, 104, info.max_abs_err);
    put_at<f32>(out, 112, f32(info.max_ratio));
    put_at<u32>(out, 116, info.flags);
    put_at<u64>(out, 120, xxh64(out, 120));
}

int smz_decode_header(const u8 in[kSmzHeaderBytes], SmzInfo *info, std::string *why) {
    auto fail = [why](const char *msg) {
        if (why) *why = msg;
        return SMESH_FAILURE;
    };
    if (!smz_host_is_little_endian()) return fail("big-endian hosts are not supported");
    if (memcmp(in, kMagic, 8) != 0) return fail("not an .smz file (bad magic)");
    if (get_at<u64>(in, 120) != xxh64(in, 120)) return fail("header checksum mismatch");
    SmzInfo i;
    i.version_major = get_at<u16>(in, 8);
    i.version_minor = get_at<u16>(in, 10);
    if (i.version_major != kSmzVersionMajor) return fail("unsupported major version");
    if (get_at<u16>(in, 12) != kSmzHeaderBytes || get_at<u16>(in, 14) != kSmzIndexEntryBytes)
        return fail("unexpected header or index entry size");
    const u8 dt = get_at<u8>(in, 16), md = get_at<u8>(in, 17), cd = get_at<u8>(in, 18), st = get_at<u8>(in, 19);
    if (!dtype_valid(dt) || !mode_valid(md) || !codec_valid(cd) || !stage_valid(st)) return fail("unknown type, mode or codec");
    if (get_at<u8>(in, 21) != 7) return fail("unsupported block size");
    ArrayPlan &p      = i.plan;
    p.dtype           = DType(dt);
    p.requested       = BoundMode(md);
    p.codec           = Codec(cd);
    p.byte_stage      = ByteStage(st);
    p.checksums       = get_at<u8>(in, 20) == 1;
    p.quant.wide      = get_at<u8>(in, 22);
    p.n_values        = get_at<u64>(in, 24);
    p.chunk_values    = get_at<u64>(in, 32);
    p.n_chunks        = get_at<u64>(in, 40);
    i.payload_bytes   = get_at<u64>(in, 48);
    p.quant.mode      = p.requested;
    p.quant.eb        = get_at<f64>(in, 56);
    p.quant.eb_abs    = get_at<f64>(in, 64);
    p.quant.origin    = get_at<f64>(in, 72);
    p.quant.L         = get_at<i32>(in, 80);
    p.quant.L_min     = get_at<i32>(in, 84);
    i.vmin            = get_at<f64>(in, 88);
    i.vmax            = get_at<f64>(in, 96);
    i.max_abs_err     = get_at<f64>(in, 104);
    i.max_ratio       = double(get_at<f32>(in, 112));
    i.flags           = get_at<u32>(in, 116);
    if (p.chunk_values == 0 || p.chunk_values > (u64(1) << 28)) return fail("bad chunk size");
    if (p.n_chunks != (p.n_values + p.chunk_values - 1) / p.chunk_values) return fail("inconsistent chunk count");
    if (p.quant.L < kLevelLo - 1 || p.quant.L > kLevelHi || p.quant.L_min < kLevelLo - 1 || p.quant.L_min > kLevelHi)
        return fail("bad grid level");
    *info = i;
    return SMESH_SUCCESS;
}

void smz_encode_index_entry(const ChunkEntry &e, u8 out[kSmzIndexEntryBytes]) {
    put_at<u64>(out, 0, e.offset);
    put_at<u64>(out, 8, e.hash);
    put_at<u32>(out, 16, e.stored_bytes);
    put_at<u32>(out, 20, e.stage1_bytes);
    put_at<u32>(out, 24, e.n_escaped);
    put_at<u8>(out, 28, u8(e.codec));
    put_at<u8>(out, 29, u8(e.byte_stage));
    put_at<u16>(out, 30, e.flags);
}

ChunkEntry smz_decode_index_entry(const u8 in[kSmzIndexEntryBytes]) {
    ChunkEntry e;
    e.offset       = get_at<u64>(in, 0);
    e.hash         = get_at<u64>(in, 8);
    e.stored_bytes = get_at<u32>(in, 16);
    e.stage1_bytes = get_at<u32>(in, 20);
    e.n_escaped    = get_at<u32>(in, 24);
    e.codec        = Codec(get_at<u8>(in, 28));
    e.byte_stage   = ByteStage(get_at<u8>(in, 29));
    e.flags        = get_at<u16>(in, 30);
    return e;
}

int smz_validate_index(const SmzInfo &info, const ChunkEntry *entries, u64 c0, u64 n, u64 expected_offset,
                       std::string *why) {
    auto fail = [why](const char *msg) {
        if (why) *why = msg;
        return SMESH_FAILURE;
    };
    u64 off = expected_offset;
    for (u64 k = 0; k < n; k++) {
        const ChunkEntry &e = entries[k];
        if (!codec_valid(u8(e.codec)) || !stage_valid(u8(e.byte_stage))) return fail("unknown codec in the index");
        if (e.offset != off) return fail("index offsets do not tile the payload");
        if (e.byte_stage == ByteStage::None && e.stage1_bytes != e.stored_bytes) return fail("inconsistent chunk sizes");
        const u64 m = info.plan.chunk_size(c0 + k);
        if (u64(e.n_escaped) > m) return fail("escape count larger than the chunk");
        off += e.stored_bytes;
    }
    if (c0 + n == info.plan.n_chunks && off != info.payload_bytes) return fail("index does not end at the payload size");
    if (off > info.payload_bytes) return fail("index runs past the payload");
    return SMESH_SUCCESS;
}

}  // namespace compression
}  // namespace smesh
