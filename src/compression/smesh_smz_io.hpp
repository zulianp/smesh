#ifndef SMESH_SMZ_IO_HPP
#define SMESH_SMZ_IO_HPP

// Serial reading and writing of .smz files (one compressed array each, see smesh_smz_format.hpp).
// The MPI versions are in src/distributed/smesh_smz_mpi.hpp. Errors return SMESH_FAILURE with
// a message on stderr; nothing here aborts.
//
// Environment: SMESH_COMPRESS_VERIFY_READ=0 skips the chunk checksums on read (default: check);
// SMESH_COMPRESS_REPORT=1 prints ratio, error and time per file written.

#include "smesh_compression_options.hpp"
#include "smesh_path.hpp"

#include <string>
#include <vector>

namespace smesh {
namespace compression {

struct SmzInfo;
struct ArrayPlan;
struct ArrayStats;
struct EncodedChunks;
struct ChunkEntry;

/// Value types the .smz functions are instantiated for.
template <class T>
static constexpr bool is_smz_type_v =
        std::is_same_v<T, i8> || std::is_same_v<T, u8> || std::is_same_v<T, i16> || std::is_same_v<T, u16> ||
        std::is_same_v<T, i32> || std::is_same_v<T, u32> || std::is_same_v<T, i64> || std::is_same_v<T, u64> ||
        std::is_same_v<T, f16> || std::is_same_v<T, f32> || std::is_same_v<T, f64>;

/// What a write did.
struct SmzReport {
    u64         n_values{0};
    u64         raw_bytes{0};   // n_values * value size
    u64         file_bytes{0};  // header + index + payload
    u64         n_chunks{0};
    u64         n_escaped{0};
    u64         n_lossy_fallback{0};
    double      max_abs_err{0};
    double      max_ratio{0};
    double      encode_s{0};
    double      write_s{0};
    std::string options;

    double ratio() const { return file_bytes ? double(raw_bytes) / double(file_bytes) : 0; }
};

/// Header of an .smz file.
int smz_read_info(const Path &path, SmzInfo *info);

/// Header and validated chunk index.
int smz_read_index(const Path &path, SmzInfo *info, std::vector<ChunkEntry> *index);

/// Compress x[0, n) into `path` (written to path.tmp, then renamed). `bound`: per-value absolute
/// bounds for BoundMode::NodeHeight, else null. Options that are off mean lossless.
template <class T>
int smz_write(const Path &path, const T *x, ptrdiff_t n, const CompressionOptions &opts, const double *bound = nullptr,
              SmzReport *report = nullptr);

/// Read a whole .smz file into a new SMESH_ALLOC'd array, converted to T.
template <class T>
int smz_read_convert(const Path &path, T **data, ptrdiff_t *n);

/// Read values [begin, end) into out (converted to T); only the chunks covering them are read.
template <class T>
int smz_read_range(const Path &path, ptrdiff_t begin, ptrdiff_t end, T *out);

/// Write an array that is already encoded (header, index, payload), atomically.
int smz_write_encoded(const Path &path, const ArrayPlan &plan, const EncodedChunks &enc, const ArrayStats &stats,
                      SmzReport *report);

/// Print a report line (rank 0 callers) if SMESH_COMPRESS_REPORT is set.
void smz_maybe_print_report(const Path &path, const SmzReport &r);

bool smz_verify_read_checksums();

/// Compress values of type T stored as type FT (e.g. idx_t written as int32).
template <class FT, class T>
int smz_write_convert(const Path &path, const T *x, ptrdiff_t n, const CompressionOptions &opts,
                      const double *bound = nullptr, SmzReport *report = nullptr) {
    if constexpr (std::is_same_v<FT, T>) {
        return smz_write<T>(path, x, n, opts, bound, report);
    } else {
        std::vector<FT> tmp(static_cast<size_t>(n));
        for (ptrdiff_t i = 0; i < n; i++) tmp[size_t(i)] = static_cast<FT>(x[i]);
        return smz_write<FT>(path, tmp.data(), n, opts, bound, report);
    }
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_SMZ_IO_HPP
