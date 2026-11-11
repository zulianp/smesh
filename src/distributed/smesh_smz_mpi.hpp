#ifndef SMESH_SMZ_MPI_HPP
#define SMESH_SMZ_MPI_HPP

// Collective reading and writing of .smz files (src/compression), producing the same file as the
// serial writer whatever the number of ranks:
//
//  * chunk boundaries are global (multiples of chunk_values), and each rank encodes whole chunks;
//  * the plan (codec, grid) comes from statistics reduced over all ranks;
//  * a chunk's bytes depend only on its values and the plan;
//  * payload offsets come from an exscan of the encoded sizes; index slices and payload slices
//    are written with MPI_File_write_at_all (rounds of at most 1 GiB), the header by rank 0.
//
// Readers take the same split as matrix.io's array_create_from_file (n/size per rank, the
// remainder to the first ranks), so the distributed mesh readers need no change, and read only
// the chunks covering their range.
//
// All functions are collective; every rank returns the same status.

#include "smesh_compression_options.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_path.hpp"
#include "smesh_smz_io.hpp"

#include <vector>

namespace smesh {
namespace compression {

/// First global value of each rank's chunk-aligned slice: rank r gets chunks
/// [r * n_chunks / size, (r + 1) * n_chunks / size). Size size + 1.
std::vector<i64> chunk_aligned_starts(i64 n_global, i64 chunk_values, int size);

/// matrix.io's uniform split of n values: rank r reads [*begin, *begin + *count).
void uniform_split(i64 n, int size, int rank, i64 *begin, i64 *count);

/// Write values that are already this rank's chunk-aligned slice (chunk_aligned_starts with
/// opts.chunk_values). `bound`: NodeHeight bounds aligned with x, else null.
template <class T>
int smz_write_slice(MPI_Comm comm, const Path &path, const T *x, ptrdiff_t n_local, ptrdiff_t n_global,
                    const CompressionOptions &opts, const double *bound = nullptr, SmzReport *report = nullptr);

/// Write a global array given as (global index, value) pairs on each rank, in any order, with
/// duplicates allowed (shared nodes: the values must agree). Optionally bound contributions
/// (contrib_map[k], bound_contrib[k]): each global value's bound is the minimum of its
/// contributions (node heights from the elements of every rank).
template <class T, class G>
int smz_write_mapped(MPI_Comm comm, const Path &path, ptrdiff_t n_local, ptrdiff_t n_global, const G *mapping,
                     const T *data, const CompressionOptions &opts, ptrdiff_t n_contrib = 0,
                     const G *contrib_map = nullptr, const double *bound_contrib = nullptr,
                     SmzReport *report = nullptr);

/// Write a global array given as contiguous blocks in rank order (as matrix.io's array_write).
template <class T>
int smz_write_contiguous(MPI_Comm comm, const Path &path, const T *x, ptrdiff_t n_local, ptrdiff_t n_global,
                         const CompressionOptions &opts, SmzReport *report = nullptr);

/// write_mapped_field for .smz paths: the value type given as an MPI datatype.
template <class G>
int smz_write_mapped_mpi_type(MPI_Comm comm, const Path &path, ptrdiff_t n_local, ptrdiff_t n_global,
                              const G *mapping, MPI_Datatype type, const void *data, const CompressionOptions &opts);

/// Read this rank's part of the uniform split into a new SMESH_ALLOC'd array (converted to T).
template <class T>
int smz_read_uniform(MPI_Comm comm, const Path &path, T **data, ptrdiff_t *n_local, ptrdiff_t *n_global);

/// Read values [begin, end) (any range per rank, possibly empty) into out.
template <class T>
int smz_read_range(MPI_Comm comm, const Path &path, ptrdiff_t begin, ptrdiff_t end, T *out);

/// Contiguous write with conversion to the file type FT.
template <class FT, class T>
int smz_write_contiguous_convert(MPI_Comm comm, const Path &path, const T *x, ptrdiff_t n_local, ptrdiff_t n_global,
                                 const CompressionOptions &opts) {
    if constexpr (std::is_same_v<FT, T>) {
        return smz_write_contiguous<T>(comm, path, x, n_local, n_global, opts);
    } else {
        std::vector<FT> tmp(static_cast<size_t>(n_local));
        for (ptrdiff_t i = 0; i < n_local; i++) tmp[size_t(i)] = static_cast<FT>(x[i]);
        return smz_write_contiguous<FT>(comm, path, tmp.data(), n_local, n_global, opts);
    }
}

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_SMZ_MPI_HPP
