#include "smesh_smz_mpi.hpp"

#include "smesh_alloc.hpp"
#include "smesh_alltoallv.impl.hpp"
#include "smesh_chunked_codec.impl.hpp"
#include "smesh_smz_format.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>

namespace smesh {
namespace compression {

namespace {
using Clock = std::chrono::steady_clock;

static constexpr i64 kRoundBytes = i64(1) << 30;

int any_error(MPI_Comm comm, int local_err) {
    int global = 0;
    MPI_Allreduce(&local_err, &global, 1, MPI_INT, MPI_MAX, comm);
    return global ? SMESH_FAILURE : SMESH_SUCCESS;
}

/// Collective write of n bytes at `offset` (n may differ, also be 0, per rank).
int write_all_rounds(MPI_File f, MPI_Offset offset, const u8 *buf, i64 n, MPI_Comm comm) {
    i64 rounds = (n + kRoundBytes - 1) / kRoundBytes, max_rounds = 0;
    MPI_Allreduce(&rounds, &max_rounds, 1, MPI_INT64_T, MPI_MAX, comm);
    int err = 0;
    for (i64 r = 0; r < max_rounds; r++) {
        const i64  b   = r * kRoundBytes;
        const i64  len = n > b ? std::min(kRoundBytes, n - b) : 0;
        MPI_Status st;
        if (MPI_File_write_at_all(f, offset + MPI_Offset(b), len ? buf + b : nullptr, int(len), MPI_BYTE, &st) !=
            MPI_SUCCESS)
            err = 1;
    }
    return err;
}

int read_all_rounds(MPI_File f, MPI_Offset offset, u8 *buf, i64 n, MPI_Comm comm) {
    i64 rounds = (n + kRoundBytes - 1) / kRoundBytes, max_rounds = 0;
    MPI_Allreduce(&rounds, &max_rounds, 1, MPI_INT64_T, MPI_MAX, comm);
    int err = 0;
    for (i64 r = 0; r < max_rounds; r++) {
        const i64  b   = r * kRoundBytes;
        const i64  len = n > b ? std::min(kRoundBytes, n - b) : 0;
        MPI_Status st;
        if (MPI_File_read_at_all(f, offset + MPI_Offset(b), len ? buf + b : nullptr, int(len), MPI_BYTE, &st) !=
            MPI_SUCCESS)
            err = 1;
        int got = 0;
        MPI_Get_count(&st, MPI_BYTE, &got);
        if (len && got != int(len)) err = 1;
    }
    return err;
}

enum class Combine { Overwrite, Min };

/// Send values to the ranks owning their global positions: rank r owns [starts[r], starts[r+1]).
/// out has starts[rank+1] - starts[rank] values and must be initialised by the caller.
template <class V, class MapFn>
int redistribute_to_ranges(MPI_Comm comm, ptrdiff_t n_local, MapFn mapping, const V *values,
                           const std::vector<i64> &starts, V *out, Combine combine) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    auto dest_of = [&](i64 g) {
        const auto it = std::upper_bound(starts.begin(), starts.end(), g);
        int        d  = int(it - starts.begin()) - 1;
        return d < 0 ? 0 : (d >= size ? size - 1 : d);
    };
    std::vector<i64> send_count(size_t(size), 0), recv_count(size_t(size), 0);
    for (ptrdiff_t i = 0; i < n_local; i++) send_count[size_t(dest_of(i64(mapping(i))))]++;
    MPI_Alltoall(send_count.data(), 1, MPI_INT64_T, recv_count.data(), 1, MPI_INT64_T, comm);
    std::vector<i64> send_displs(size_t(size), 0), recv_displs(size_t(size), 0);
    for (int r = 1; r < size; r++) {
        send_displs[size_t(r)] = send_displs[size_t(r - 1)] + send_count[size_t(r - 1)];
        recv_displs[size_t(r)] = recv_displs[size_t(r - 1)] + recv_count[size_t(r - 1)];
    }
    const i64        n_recv = recv_displs[size_t(size - 1)] + recv_count[size_t(size - 1)];
    std::vector<i64> send_idx(size_t(std::max<ptrdiff_t>(n_local, 1)));
    std::vector<V>   send_val(size_t(std::max<ptrdiff_t>(n_local, 1)));
    std::vector<i64> fill(send_displs);
    for (ptrdiff_t i = 0; i < n_local; i++) {
        const i64    g = i64(mapping(i));
        const size_t k = size_t(fill[size_t(dest_of(g))]++);
        send_idx[k]    = g;
        send_val[k]    = values[i];
    }
    std::vector<i64> recv_idx(size_t(std::max<i64>(n_recv, 1)));
    std::vector<V>   recv_val(size_t(std::max<i64>(n_recv, 1)));
    // Values travel as bytes: counts and displacements in bytes.
    std::vector<i64> sc_b(send_count), sd_b(send_displs), rc_b(recv_count), rd_b(recv_displs);
    for (int r = 0; r < size; r++) {
        sc_b[size_t(r)] *= i64(sizeof(V));
        sd_b[size_t(r)] *= i64(sizeof(V));
        rc_b[size_t(r)] *= i64(sizeof(V));
        rd_b[size_t(r)] *= i64(sizeof(V));
    }
    const i64 max_chunk = std::numeric_limits<i32>::max() / 2;
    all_to_allv_64_b(send_idx.data(), send_count.data(), send_displs.data(), MPI_INT64_T, recv_idx.data(),
                     recv_count.data(), recv_displs.data(), MPI_INT64_T, comm, max_chunk);
    all_to_allv_64_b(send_val.data(), sc_b.data(), sd_b.data(), MPI_BYTE, recv_val.data(), rc_b.data(), rd_b.data(),
                     MPI_BYTE, comm, max_chunk);
    const i64 begin = starts[size_t(rank)], m = starts[size_t(rank) + 1] - begin;
    int       err   = 0;
    for (i64 k = 0; k < n_recv; k++) {
        const i64 d = recv_idx[size_t(k)] - begin;
        if (d < 0 || d >= m) {
            err = 1;
            continue;
        }
        if (combine == Combine::Min) {
            out[d] = recv_val[size_t(k)] < out[d] ? recv_val[size_t(k)] : out[d];
        } else {
            out[d] = recv_val[size_t(k)];
        }
    }
    return err;
}

struct CollectiveFile {
    MPI_File f{MPI_FILE_NULL};
    ~CollectiveFile() {
        if (f != MPI_FILE_NULL) MPI_File_close(&f);
    }
};

/// Open an .smz file collectively and validate its header (read by rank 0).
int open_collective(MPI_Comm comm, const Path &path, CollectiveFile *file, SmzInfo *info) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    u8  hdr[kSmzHeaderBytes] = {0};
    int err                  = 0;
    if (MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file->f) != MPI_SUCCESS) {
        file->f = MPI_FILE_NULL;
        if (!rank) fprintf(stderr, "smz: cannot open %s\n", path.c_str());
        return SMESH_FAILURE;  // MPI_File_open is collective: all ranks failed
    }
    MPI_Offset size = 0;
    if (!rank) {
        MPI_Status st;
        int        got = 0;
        if (MPI_File_read_at(file->f, 0, hdr, int(kSmzHeaderBytes), MPI_BYTE, &st) != MPI_SUCCESS ||
            MPI_Get_count(&st, MPI_BYTE, &got) != MPI_SUCCESS || got != int(kSmzHeaderBytes)) {
            err = 1;
        }
        MPI_File_get_size(file->f, &size);
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, comm);
    if (err) {
        if (!rank) fprintf(stderr, "smz: %s is too short for a header\n", path.c_str());
        return SMESH_FAILURE;
    }
    MPI_Bcast(hdr, int(kSmzHeaderBytes), MPI_BYTE, 0, comm);
    i64 size64 = i64(size);
    MPI_Bcast(&size64, 1, MPI_INT64_T, 0, comm);
    std::string why;
    if (smz_decode_header(hdr, info, &why) != SMESH_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: %s: %s\n", path.c_str(), why.c_str());
        return SMESH_FAILURE;
    }
    if (u64(size64) != info->file_bytes()) {
        if (!rank) fprintf(stderr, "smz: %s: file size does not match its header (truncated?)\n", path.c_str());
        return SMESH_FAILURE;
    }
    return SMESH_SUCCESS;
}

template <class FT, class T>
int decode_collective(MPI_Comm comm, MPI_File f, const SmzInfo &info, u64 begin, u64 end, T *out, const Path &path) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    u64 c0 = 0, c1 = 0;
    if (end > begin) chunks_covering(info.plan, begin, end, &c0, &c1);
    const u64       nc = c1 - c0;
    std::vector<u8> raw(size_t(nc) * kSmzIndexEntryBytes);
    int             err = read_all_rounds(f, MPI_Offset(info.index_offset() + kSmzIndexEntryBytes * c0), raw.data(),
                                          i64(raw.size()), comm);
    // decode_chunks indexes entries by absolute chunk number: keep a full-size index, filled in [c0, c1).
    std::vector<ChunkEntry> index(static_cast<size_t>(c1));
    for (u64 k = 0; k < nc; k++) index[size_t(c0 + k)] = smz_decode_index_entry(raw.data() + k * kSmzIndexEntryBytes);
    std::string why;
    if (!err && nc &&
        smz_validate_index(info, index.data() + c0, c0, nc, index[size_t(c0)].offset, &why) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: %s (rank %d): %s\n", path.c_str(), rank, why.c_str());
        err = 1;
    }
    u64 first = 0, last = 0;
    if (!err && nc) {
        first = index[size_t(c0)].offset;
        last  = index[size_t(c1 - 1)].offset + index[size_t(c1 - 1)].stored_bytes;
    }
    std::vector<u8> payload(size_t(last - first) + kSlack);
    err |= read_all_rounds(f, MPI_Offset(info.payload_offset() + first), payload.data(), i64(last - first), comm);
    if (!err && nc &&
        decode_chunks<FT, T>(info.plan, payload.data(), index.data(), c0, c1, begin, end, out,
                             smz_verify_read_checksums()) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: %s (rank %d): corrupt data\n", path.c_str(), rank);
        err = 1;
    }
    return any_error(comm, err);
}

template <class T>
int decode_collective_any(MPI_Comm comm, MPI_File f, const SmzInfo &info, u64 begin, u64 end, T *out,
                          const Path &path) {
    switch (info.plan.dtype) {
        case DType::I8: return decode_collective<i8, T>(comm, f, info, begin, end, out, path);
        case DType::U8: return decode_collective<u8, T>(comm, f, info, begin, end, out, path);
        case DType::I16: return decode_collective<i16, T>(comm, f, info, begin, end, out, path);
        case DType::U16: return decode_collective<u16, T>(comm, f, info, begin, end, out, path);
        case DType::I32: return decode_collective<i32, T>(comm, f, info, begin, end, out, path);
        case DType::U32: return decode_collective<u32, T>(comm, f, info, begin, end, out, path);
        case DType::I64: return decode_collective<i64, T>(comm, f, info, begin, end, out, path);
        case DType::U64: return decode_collective<u64, T>(comm, f, info, begin, end, out, path);
        case DType::F16: return decode_collective<f16, T>(comm, f, info, begin, end, out, path);
        case DType::F32: return decode_collective<f32, T>(comm, f, info, begin, end, out, path);
        case DType::F64: return decode_collective<f64, T>(comm, f, info, begin, end, out, path);
    }
    return SMESH_FAILURE;
}
}  // namespace

std::vector<i64> chunk_aligned_starts(i64 n_global, i64 chunk_values, int size) {
    std::vector<i64> s(size_t(size) + 1);
    const i64        nc = (n_global + chunk_values - 1) / chunk_values;
    for (int r = 0; r <= size; r++) {
        const i64 c       = nc * r / size;
        s[size_t(r)]      = std::min(n_global, c * chunk_values);
    }
    return s;
}

void uniform_split(i64 n, int size, int rank, i64 *begin, i64 *count) {
    const i64 u = n / size, rem = n - u * size;
    *count      = u + (rank < rem ? 1 : 0);
    *begin      = rank * u + std::min<i64>(rank, rem);
}

template <class T>
int smz_write_slice(MPI_Comm comm, const Path &path, const T *x, ptrdiff_t n_local, ptrdiff_t n_global,
                    const CompressionOptions &opts_in, const double *bound, SmzReport *report) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    CompressionOptions opts = opts_in.enabled ? opts_in : CompressionOptions::lossless();
    const auto         t0   = Clock::now();

    const std::vector<i64> starts = chunk_aligned_starts(n_global, opts.chunk_values, size);
    int                    err    = (starts[size_t(rank) + 1] - starts[size_t(rank)] != n_local) ? 1 : 0;
    if (any_error(comm, err) != SMESH_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: %s: the local slices are not the chunk-aligned split\n", path.c_str());
        return SMESH_FAILURE;
    }

    // Global statistics -> the same plan on every rank.
    ArrayStats s = local_stats<T>(x, n_local, bound);
    double     mins[3] = {s.vmin, s.amin_nz, s.bmin}, maxs[2] = {s.vmax, s.amax};
    double     gmins[3], gmaxs[2];
    u64        nnf = 0;
    MPI_Allreduce(mins, gmins, 3, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(maxs, gmaxs, 2, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&s.n_nonfinite, &nnf, 1, MPI_UINT64_T, MPI_SUM, comm);
    ArrayStats g;
    g.vmin = gmins[0], g.amin_nz = gmins[1], g.bmin = gmins[2], g.vmax = gmaxs[0], g.amax = gmaxs[1];
    g.n_nonfinite = nnf;
    ArrayPlan plan;
    if (any_error(comm, plan_array(opts, DTypeOf<T>::value, u64(n_global), g, &plan) != SMESH_SUCCESS) !=
        SMESH_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: cannot compress %s with '%s'\n", path.c_str(), opts.to_string().c_str());
        return SMESH_FAILURE;
    }

    EncodedChunks enc;
    err = encode_chunks<T>(plan, x, u64(starts[size_t(rank)]), u64(n_local), bound, &enc) != SMESH_SUCCESS;
    if (any_error(comm, err) != SMESH_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: encoding failed for %s\n", path.c_str());
        return SMESH_FAILURE;
    }
    const double encode_s = std::chrono::duration<double>(Clock::now() - t0).count();
    const auto   t1       = Clock::now();

    // Payload offsets and global results.
    u64 bytes = enc.payload.size(), base = 0, total = 0, esc = 0, n_esc = 0, fb = 0, n_fb = 0;
    MPI_Exscan(&bytes, &base, 1, MPI_UINT64_T, MPI_SUM, comm);
    if (!rank) base = 0;
    MPI_Allreduce(&bytes, &total, 1, MPI_UINT64_T, MPI_SUM, comm);
    for (const auto &e : enc.index) esc += e.n_escaped;
    fb = enc.n_lossy_fallback;
    MPI_Allreduce(&esc, &n_esc, 1, MPI_UINT64_T, MPI_SUM, comm);
    MPI_Allreduce(&fb, &n_fb, 1, MPI_UINT64_T, MPI_SUM, comm);
    double lmax[2] = {enc.stats.max_abs_err, enc.stats.max_ratio}, gmax[2];
    MPI_Allreduce(lmax, gmax, 2, MPI_DOUBLE, MPI_MAX, comm);

    SmzInfo info;
    info.plan          = plan;
    info.payload_bytes = total;
    info.vmin          = g.vmin;
    info.vmax          = g.vmax;
    info.max_abs_err   = gmax[0];
    info.max_ratio     = gmax[1];
    info.flags         = plan.verify_write ? 1u : 0u;

    const u64       c0 = u64(starts[size_t(rank)]) / plan.chunk_values;
    std::vector<u8> idx(enc.index.size() * kSmzIndexEntryBytes);
    for (size_t k = 0; k < enc.index.size(); k++) {
        ChunkEntry e = enc.index[k];
        e.offset += base;
        smz_encode_index_entry(e, idx.data() + k * kSmzIndexEntryBytes);
    }

    const std::string tmp = path.to_string() + ".tmp";
    MPI_File          f;
    if (MPI_File_open(comm, tmp.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f) != MPI_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: cannot write %s\n", tmp.c_str());
        return SMESH_FAILURE;
    }
    err = MPI_File_set_size(f, MPI_Offset(info.file_bytes())) != MPI_SUCCESS;
    if (!rank) {
        u8 hdr[kSmzHeaderBytes];
        smz_encode_header(info, hdr);
        MPI_Status st;
        err |= MPI_File_write_at(f, 0, hdr, int(kSmzHeaderBytes), MPI_BYTE, &st) != MPI_SUCCESS;
    }
    err |= write_all_rounds(f, MPI_Offset(info.index_offset() + kSmzIndexEntryBytes * c0), idx.data(), i64(idx.size()),
                            comm);
    err |= write_all_rounds(f, MPI_Offset(info.payload_offset() + base), enc.payload.data(), i64(enc.payload.size()),
                            comm);
    err |= MPI_File_close(&f) != MPI_SUCCESS;
    if (any_error(comm, err) != SMESH_SUCCESS) {
        if (!rank) {
            std::remove(tmp.c_str());
            fprintf(stderr, "smz: failed to write %s\n", path.c_str());
        }
        return SMESH_FAILURE;
    }
    if (!rank) err = std::rename(tmp.c_str(), path.c_str()) != 0;
    MPI_Bcast(&err, 1, MPI_INT, 0, comm);
    if (err) {
        if (!rank) fprintf(stderr, "smz: cannot rename %s\n", tmp.c_str());
        return SMESH_FAILURE;
    }

    SmzReport r;
    r.n_values         = plan.n_values;
    r.raw_bytes        = plan.n_values * dtype_bytes(plan.dtype);
    r.file_bytes       = info.file_bytes();
    r.n_chunks         = plan.n_chunks;
    r.n_escaped        = n_esc;
    r.n_lossy_fallback = n_fb;
    r.max_abs_err      = gmax[0];
    r.max_ratio        = gmax[1];
    r.encode_s         = encode_s;
    r.write_s          = std::chrono::duration<double>(Clock::now() - t1).count();
    r.options          = opts.to_string();
    if (!rank) smz_maybe_print_report(path, r);
    if (report) *report = r;
    return SMESH_SUCCESS;
}

template <class T, class G>
int smz_write_mapped(MPI_Comm comm, const Path &path, ptrdiff_t n_local, ptrdiff_t n_global, const G *mapping,
                     const T *data, const CompressionOptions &opts_in, ptrdiff_t n_contrib, const G *contrib_map,
                     const double *bound_contrib, SmzReport *report) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const CompressionOptions opts   = opts_in.enabled ? opts_in : CompressionOptions::lossless();
    const std::vector<i64>   starts = chunk_aligned_starts(n_global, opts.chunk_values, size);
    const i64                m      = starts[size_t(rank) + 1] - starts[size_t(rank)];
    std::vector<T>           slice(size_t(m), T(0));
    int err = redistribute_to_ranges<T>(comm, n_local, [mapping](ptrdiff_t i) { return mapping[i]; }, data, starts,
                                        slice.data(), Combine::Overwrite);
    int any_bound = bound_contrib ? 1 : 0, all_bound = 0;
    MPI_Allreduce(&any_bound, &all_bound, 1, MPI_INT, MPI_MAX, comm);
    std::vector<double> bound;
    if (all_bound) {
        bound.assign(size_t(m), INFINITY);
        err |= redistribute_to_ranges<double>(
                comm, bound_contrib ? n_contrib : 0, [contrib_map](ptrdiff_t i) { return contrib_map[i]; },
                bound_contrib, starts, bound.data(), Combine::Min);
    }
    if (any_error(comm, err) != SMESH_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: %s: global indices out of range\n", path.c_str());
        return SMESH_FAILURE;
    }
    return smz_write_slice<T>(comm, path, slice.data(), ptrdiff_t(m), n_global, opts,
                              all_bound ? bound.data() : nullptr, report);
}

template <class T>
int smz_write_contiguous(MPI_Comm comm, const Path &path, const T *x, ptrdiff_t n_local, ptrdiff_t n_global,
                         const CompressionOptions &opts_in, SmzReport *report) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const CompressionOptions opts   = opts_in.enabled ? opts_in : CompressionOptions::lossless();
    const std::vector<i64>   starts = chunk_aligned_starts(n_global, opts.chunk_values, size);
    i64                      n = n_local, offset = 0;
    MPI_Exscan(&n, &offset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (!rank) offset = 0;
    const i64      m = starts[size_t(rank) + 1] - starts[size_t(rank)];
    std::vector<T> slice(size_t(m), T(0));
    const int err = redistribute_to_ranges<T>(comm, n_local, [offset](ptrdiff_t i) { return offset + i64(i); }, x,
                                              starts, slice.data(), Combine::Overwrite);
    if (any_error(comm, err) != SMESH_SUCCESS) {
        if (!rank) fprintf(stderr, "smz: %s: local sizes do not add up to %td\n", path.c_str(), n_global);
        return SMESH_FAILURE;
    }
    return smz_write_slice<T>(comm, path, slice.data(), ptrdiff_t(m), n_global, opts, nullptr, report);
}

template <class G>
int smz_write_mapped_mpi_type(MPI_Comm comm, const Path &path, ptrdiff_t n_local, ptrdiff_t n_global,
                              const G *mapping, MPI_Datatype type, const void *data, const CompressionOptions &opts) {
#define SMESH_SMZ_DISPATCH(T)                                                                                          \
    return smz_write_mapped<T, G>(comm, path, n_local, n_global, mapping, static_cast<const T *>(data), opts)
    if (type == MPI_FLOAT) SMESH_SMZ_DISPATCH(f32);
    if (type == MPI_DOUBLE) SMESH_SMZ_DISPATCH(f64);
    if (type == MPI_INT || type == MPI_INT32_T) SMESH_SMZ_DISPATCH(i32);
    if (type == MPI_LONG_LONG || type == MPI_INT64_T || (sizeof(long) == 8 && type == MPI_LONG))
        SMESH_SMZ_DISPATCH(i64);
    if (type == MPI_SHORT || type == MPI_INT16_T) SMESH_SMZ_DISPATCH(i16);
    if (type == MPI_UNSIGNED_CHAR || type == MPI_UINT8_T) SMESH_SMZ_DISPATCH(u8);
    if (type == MPI_UNSIGNED_SHORT || type == MPI_UINT16_T) SMESH_SMZ_DISPATCH(u16);
    if (type == MPI_UNSIGNED || type == MPI_UINT32_T) SMESH_SMZ_DISPATCH(u32);
    if (type == MPI_UINT64_T || type == MPI_UNSIGNED_LONG_LONG) SMESH_SMZ_DISPATCH(u64);
    if (type == mpi_type<f16>()) SMESH_SMZ_DISPATCH(f16);
#undef SMESH_SMZ_DISPATCH
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (!rank) fprintf(stderr, "smz: %s: unsupported MPI datatype for compression\n", path.c_str());
    return SMESH_FAILURE;
}

template <class T>
int smz_read_range(MPI_Comm comm, const Path &path, ptrdiff_t begin, ptrdiff_t end, T *out) {
    CollectiveFile f;
    SmzInfo        info;
    if (open_collective(comm, path, &f, &info) != SMESH_SUCCESS) return SMESH_FAILURE;
    const int bad = begin < 0 || end < begin || u64(end) > info.plan.n_values;
    if (any_error(comm, bad) != SMESH_SUCCESS) {
        fprintf(stderr, "smz: %s: range outside the array\n", path.c_str());
        return SMESH_FAILURE;
    }
    return decode_collective_any<T>(comm, f.f, info, u64(begin), u64(end), out, path);
}

template <class T>
int smz_read_uniform(MPI_Comm comm, const Path &path, T **data, ptrdiff_t *n_local, ptrdiff_t *n_global) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    *data = nullptr;
    CollectiveFile f;
    SmzInfo        info;
    if (open_collective(comm, path, &f, &info) != SMESH_SUCCESS) return SMESH_FAILURE;
    i64 begin, count;
    uniform_split(i64(info.plan.n_values), size, rank, &begin, &count);
    T *out = static_cast<T *>(SMESH_ALLOC(size_t(count > 0 ? count : 1) * sizeof(T)));
    if (decode_collective_any<T>(comm, f.f, info, u64(begin), u64(begin + count), out, path) != SMESH_SUCCESS) {
        SMESH_FREE(out);
        return SMESH_FAILURE;
    }
    *data     = out;
    *n_local  = ptrdiff_t(count);
    *n_global = ptrdiff_t(info.plan.n_values);
    return SMESH_SUCCESS;
}

#define SMESH_SMZ_MPI_INSTANTIATE(T)                                                                                   \
    template int smz_write_slice<T>(MPI_Comm, const Path &, const T *, ptrdiff_t, ptrdiff_t,                          \
                                    const CompressionOptions &, const double *, SmzReport *);                         \
    template int smz_write_mapped<T, i32>(MPI_Comm, const Path &, ptrdiff_t, ptrdiff_t, const i32 *, const T *,       \
                                          const CompressionOptions &, ptrdiff_t, const i32 *, const double *,         \
                                          SmzReport *);                                                               \
    template int smz_write_mapped<T, i64>(MPI_Comm, const Path &, ptrdiff_t, ptrdiff_t, const i64 *, const T *,       \
                                          const CompressionOptions &, ptrdiff_t, const i64 *, const double *,         \
                                          SmzReport *);                                                               \
    template int smz_write_contiguous<T>(MPI_Comm, const Path &, const T *, ptrdiff_t, ptrdiff_t,                     \
                                         const CompressionOptions &, SmzReport *);                                    \
    template int smz_read_uniform<T>(MPI_Comm, const Path &, T **, ptrdiff_t *, ptrdiff_t *);                         \
    template int smz_read_range<T>(MPI_Comm, const Path &, ptrdiff_t, ptrdiff_t, T *);

SMESH_SMZ_MPI_INSTANTIATE(i8)
SMESH_SMZ_MPI_INSTANTIATE(u8)
SMESH_SMZ_MPI_INSTANTIATE(i16)
SMESH_SMZ_MPI_INSTANTIATE(u16)
SMESH_SMZ_MPI_INSTANTIATE(i32)
SMESH_SMZ_MPI_INSTANTIATE(u32)
SMESH_SMZ_MPI_INSTANTIATE(i64)
SMESH_SMZ_MPI_INSTANTIATE(u64)
SMESH_SMZ_MPI_INSTANTIATE(f16)
SMESH_SMZ_MPI_INSTANTIATE(f32)
SMESH_SMZ_MPI_INSTANTIATE(f64)

template int smz_write_mapped_mpi_type<i32>(MPI_Comm, const Path &, ptrdiff_t, ptrdiff_t, const i32 *, MPI_Datatype,
                                            const void *, const CompressionOptions &);
template int smz_write_mapped_mpi_type<i64>(MPI_Comm, const Path &, ptrdiff_t, ptrdiff_t, const i64 *, MPI_Datatype,
                                            const void *, const CompressionOptions &);

}  // namespace compression
}  // namespace smesh
