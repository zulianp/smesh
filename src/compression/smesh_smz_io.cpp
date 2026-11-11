#include "smesh_smz_io.hpp"

#include "smesh_chunked_codec.impl.hpp"
#include "smesh_smz_format.hpp"

#include "smesh_alloc.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>

namespace smesh {
namespace compression {

namespace {
using Clock = std::chrono::steady_clock;

struct FileCloser {
    void operator()(FILE *f) const {
        if (f) fclose(f);
    }
};
using FilePtr = std::unique_ptr<FILE, FileCloser>;

bool read_at(FILE *f, u64 off, void *dst, size_t n) {
    if (n == 0) return true;
    if (fseeko(f, off_t(off), SEEK_SET) != 0) return false;
    return fread(dst, 1, n, f) == n;
}

/// Open, read and validate header and the full index.
int open_smz(const Path &path, FilePtr *file, SmzInfo *info, std::vector<ChunkEntry> *index) {
    FilePtr f(fopen(path.c_str(), "rb"));
    if (!f) {
        fprintf(stderr, "smz: cannot open %s\n", path.c_str());
        return SMESH_FAILURE;
    }
    u8 hdr[kSmzHeaderBytes];
    if (!read_at(f.get(), 0, hdr, kSmzHeaderBytes)) {
        fprintf(stderr, "smz: %s is too short for a header\n", path.c_str());
        return SMESH_FAILURE;
    }
    std::string why;
    if (smz_decode_header(hdr, info, &why) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: %s: %s\n", path.c_str(), why.c_str());
        return SMESH_FAILURE;
    }
    if (fseeko(f.get(), 0, SEEK_END) != 0 || u64(ftello(f.get())) != info->file_bytes()) {
        fprintf(stderr, "smz: %s: file size does not match its header (truncated?)\n", path.c_str());
        return SMESH_FAILURE;
    }
    if (index) {
        std::vector<u8> raw(size_t(info->plan.n_chunks) * kSmzIndexEntryBytes);
        if (!read_at(f.get(), info->index_offset(), raw.data(), raw.size())) {
            fprintf(stderr, "smz: %s: cannot read the index\n", path.c_str());
            return SMESH_FAILURE;
        }
        index->resize(size_t(info->plan.n_chunks));
        for (size_t c = 0; c < index->size(); c++) (*index)[c] = smz_decode_index_entry(raw.data() + c * kSmzIndexEntryBytes);
        if (smz_validate_index(*info, index->data(), 0, info->plan.n_chunks, 0, &why) != SMESH_SUCCESS) {
            fprintf(stderr, "smz: %s: %s\n", path.c_str(), why.c_str());
            return SMESH_FAILURE;
        }
    }
    *file = std::move(f);
    return SMESH_SUCCESS;
}

template <class FT, class T>
int decode_range(const SmzInfo &info, const std::vector<ChunkEntry> &index, FILE *f, u64 begin, u64 end, T *out,
                 const Path &path) {
    u64 c0, c1;
    chunks_covering(info.plan, begin, end, &c0, &c1);
    if (c1 <= c0) return SMESH_SUCCESS;
    const u64       first = index[c0].offset;
    const u64       last  = index[c1 - 1].offset + index[c1 - 1].stored_bytes;
    std::vector<u8> payload(size_t(last - first) + kSlack);
    if (!read_at(f, info.payload_offset() + first, payload.data(), size_t(last - first))) {
        fprintf(stderr, "smz: %s: cannot read the payload\n", path.c_str());
        return SMESH_FAILURE;
    }
    if (decode_chunks<FT, T>(info.plan, payload.data(), index.data(), c0, c1, begin, end, out,
                             smz_verify_read_checksums()) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: %s: corrupt data\n", path.c_str());
        return SMESH_FAILURE;
    }
    return SMESH_SUCCESS;
}

template <class T>
int decode_range_any(const SmzInfo &info, const std::vector<ChunkEntry> &index, FILE *f, u64 begin, u64 end, T *out,
                     const Path &path) {
    switch (info.plan.dtype) {
        case DType::I8: return decode_range<i8, T>(info, index, f, begin, end, out, path);
        case DType::U8: return decode_range<u8, T>(info, index, f, begin, end, out, path);
        case DType::I16: return decode_range<i16, T>(info, index, f, begin, end, out, path);
        case DType::U16: return decode_range<u16, T>(info, index, f, begin, end, out, path);
        case DType::I32: return decode_range<i32, T>(info, index, f, begin, end, out, path);
        case DType::U32: return decode_range<u32, T>(info, index, f, begin, end, out, path);
        case DType::I64: return decode_range<i64, T>(info, index, f, begin, end, out, path);
        case DType::U64: return decode_range<u64, T>(info, index, f, begin, end, out, path);
        case DType::F16: return decode_range<f16, T>(info, index, f, begin, end, out, path);
        case DType::F32: return decode_range<f32, T>(info, index, f, begin, end, out, path);
        case DType::F64: return decode_range<f64, T>(info, index, f, begin, end, out, path);
    }
    return SMESH_FAILURE;
}
}  // namespace

bool smz_verify_read_checksums() {
    const char *v = std::getenv("SMESH_COMPRESS_VERIFY_READ");
    return !(v && *v && std::atoi(v) == 0);
}

void smz_maybe_print_report(const Path &path, const SmzReport &r) {
    const char *v = std::getenv("SMESH_COMPRESS_REPORT");
    if (!(v && *v && std::atoi(v) != 0)) return;
    printf("smz: %s  %s  %llu values  %.1f -> %.1f MB (%.2fx)  max err %.3g (%.3f of bound)  escaped %llu  "
           "lossy->lossless chunks %llu  encode %.3f s  write %.3f s\n",
           path.c_str(), r.options.c_str(), (unsigned long long)r.n_values, double(r.raw_bytes) / 1e6,
           double(r.file_bytes) / 1e6, r.ratio(), r.max_abs_err, r.max_ratio, (unsigned long long)r.n_escaped,
           (unsigned long long)r.n_lossy_fallback, r.encode_s, r.write_s);
    fflush(stdout);
}

int smz_read_info(const Path &path, SmzInfo *info) {
    FilePtr f;
    return open_smz(path, &f, info, nullptr);
}

int smz_read_index(const Path &path, SmzInfo *info, std::vector<ChunkEntry> *index) {
    FilePtr f;
    return open_smz(path, &f, info, index);
}

int smz_write_encoded(const Path &path, const ArrayPlan &plan, const EncodedChunks &enc, const ArrayStats &stats,
                      SmzReport *report) {
    const auto t0 = Clock::now();
    SmzInfo    info;
    info.plan          = plan;
    info.payload_bytes = enc.payload.size();
    info.vmin          = stats.vmin;
    info.vmax          = stats.vmax;
    info.max_abs_err   = enc.stats.max_abs_err;
    info.max_ratio     = enc.stats.max_ratio;
    info.flags         = plan.verify_write ? 1u : 0u;

    u8 hdr[kSmzHeaderBytes];
    smz_encode_header(info, hdr);
    std::vector<u8> idx(enc.index.size() * kSmzIndexEntryBytes);
    for (size_t c = 0; c < enc.index.size(); c++) smz_encode_index_entry(enc.index[c], idx.data() + c * kSmzIndexEntryBytes);

    const std::string tmp = path.to_string() + ".tmp";
    FILE             *f   = fopen(tmp.c_str(), "wb");
    if (!f) {
        fprintf(stderr, "smz: cannot write %s\n", tmp.c_str());
        return SMESH_FAILURE;
    }
    bool ok = fwrite(hdr, 1, kSmzHeaderBytes, f) == kSmzHeaderBytes;
    ok      = ok && fwrite(idx.data(), 1, idx.size(), f) == idx.size();
    ok      = ok && (enc.payload.empty() || fwrite(enc.payload.data(), 1, enc.payload.size(), f) == enc.payload.size());
    ok      = (fclose(f) == 0) && ok;
    if (!ok || std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::remove(tmp.c_str());
        fprintf(stderr, "smz: failed to write %s\n", path.c_str());
        return SMESH_FAILURE;
    }
    if (report) {
        u64 esc = 0;
        for (const auto &e : enc.index) esc += e.n_escaped;
        report->n_values         = plan.n_values;
        report->raw_bytes        = plan.n_values * dtype_bytes(plan.dtype);
        report->file_bytes       = info.file_bytes();
        report->n_chunks         = plan.n_chunks;
        report->n_escaped        = esc;
        report->n_lossy_fallback = enc.n_lossy_fallback;
        report->max_abs_err      = enc.stats.max_abs_err;
        report->max_ratio        = enc.stats.max_ratio;
        report->write_s          = std::chrono::duration<double>(Clock::now() - t0).count();
    }
    return SMESH_SUCCESS;
}

template <class T>
int smz_write(const Path &path, const T *x, ptrdiff_t n, const CompressionOptions &opts_in, const double *bound,
              SmzReport *report) {
    CompressionOptions opts = opts_in;
    if (!opts.enabled) opts = CompressionOptions::lossless();
    const auto       t0    = Clock::now();
    const ArrayStats stats = local_stats<T>(x, n, bound);
    ArrayPlan        plan;
    if (plan_array(opts, DTypeOf<T>::value, u64(n), stats, &plan) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: cannot compress %s with '%s'\n", path.c_str(), opts.to_string().c_str());
        return SMESH_FAILURE;
    }
    EncodedChunks enc;
    if (encode_chunks<T>(plan, x, 0, u64(n), bound, &enc) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: encoding failed for %s\n", path.c_str());
        return SMESH_FAILURE;
    }
    const double encode_s = std::chrono::duration<double>(Clock::now() - t0).count();
    SmzReport    r;
    if (smz_write_encoded(path, plan, enc, stats, &r) != SMESH_SUCCESS) return SMESH_FAILURE;
    r.encode_s = encode_s;
    r.options  = opts.to_string();
    smz_maybe_print_report(path, r);
    if (report) *report = r;
    return SMESH_SUCCESS;
}

template <class T>
int smz_read_convert(const Path &path, T **data, ptrdiff_t *n) {
    *data = nullptr;
    *n    = 0;
    FilePtr                 f;
    SmzInfo                 info;
    std::vector<ChunkEntry> index;
    if (open_smz(path, &f, &info, &index) != SMESH_SUCCESS) return SMESH_FAILURE;
    const u64 nv = info.plan.n_values;
    T        *out = static_cast<T *>(SMESH_ALLOC(size_t(nv ? nv : 1) * sizeof(T)));
    if (decode_range_any<T>(info, index, f.get(), 0, nv, out, path) != SMESH_SUCCESS) {
        SMESH_FREE(out);
        return SMESH_FAILURE;
    }
    *data = out;
    *n    = ptrdiff_t(nv);
    return SMESH_SUCCESS;
}

template <class T>
int smz_read_range(const Path &path, ptrdiff_t begin, ptrdiff_t end, T *out) {
    FilePtr                 f;
    SmzInfo                 info;
    std::vector<ChunkEntry> index;
    if (open_smz(path, &f, &info, &index) != SMESH_SUCCESS) return SMESH_FAILURE;
    if (begin < 0 || end < begin || u64(end) > info.plan.n_values) {
        fprintf(stderr, "smz: %s: range [%td, %td) outside [0, %llu)\n", path.c_str(), begin, end,
                (unsigned long long)info.plan.n_values);
        return SMESH_FAILURE;
    }
    return decode_range_any<T>(info, index, f.get(), u64(begin), u64(end), out, path);
}

#define SMESH_SMZ_IO_INSTANTIATE(T)                                                                                    \
    template int smz_write<T>(const Path &, const T *, ptrdiff_t, const CompressionOptions &, const double *,         \
                              SmzReport *);                                                                           \
    template int smz_read_convert<T>(const Path &, T **, ptrdiff_t *);                                               \
    template int smz_read_range<T>(const Path &, ptrdiff_t, ptrdiff_t, T *);

SMESH_SMZ_IO_INSTANTIATE(i8)
SMESH_SMZ_IO_INSTANTIATE(u8)
SMESH_SMZ_IO_INSTANTIATE(i16)
SMESH_SMZ_IO_INSTANTIATE(u16)
SMESH_SMZ_IO_INSTANTIATE(i32)
SMESH_SMZ_IO_INSTANTIATE(u32)
SMESH_SMZ_IO_INSTANTIATE(i64)
SMESH_SMZ_IO_INSTANTIATE(u64)
SMESH_SMZ_IO_INSTANTIATE(f16)
SMESH_SMZ_IO_INSTANTIATE(f32)
SMESH_SMZ_IO_INSTANTIATE(f64)

}  // namespace compression
}  // namespace smesh
