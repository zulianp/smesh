#ifndef SMESH_COMPRESSION_OPTIONS_HPP
#define SMESH_COMPRESSION_OPTIONS_HPP

// How one array is compressed. Text form (environment variables, command lines):
//
//   off | lossless | abs:E | range:E | pwrel:E | h:E   [,zstd[:LEVEL]] [,lz4] [,shuffle]
//                                                      [,chunk:N] [,nochecksum] [,verify]
//
//   lossless   bit-exact (integers and floats)
//   abs:E      |x' - x| <= E
//   range:E    |x' - x| <= E * (max - min) of the array
//   pwrel:E    |x' - x| <= E * |x|, 0 < E < 1
//   h:E        coordinates: |x' - x| <= E * h_n, h_n the smallest element height at node n
//   zstd, lz4  general-purpose stage after the codec, kept per chunk when it saves >= 5%
//   shuffle    lossless floats as byte planes (use with zstd / lz4)
//   chunk:N    values per chunk (default 65536, or SMESH_COMPRESS_CHUNK)
//   verify     decode every chunk again after encoding it and check it
//
// Defaults: off; "on" without a mode is lossless with smesh's own codecs and no byte stage
// (a file written with zstd needs zstd to be read back).

#include "smesh_compression_types.hpp"

#include <string>
#include <string_view>

namespace smesh {
namespace compression {

struct CompressionOptions {
    bool      enabled{false};
    BoundMode mode{BoundMode::Lossless};
    double    eb{0};
    ByteStage byte_stage{ByteStage::None};
    int       zstd_level{3};
    bool      shuffle{false};
    i64       chunk_values{i64(1) << 16};
    bool      checksums{true};
    bool      verify_write{false};
    double    max_escape_fraction{0.5};  // lossy chunks with more escapes are written losslessly

    static CompressionOptions lossless();

    /// Parse the text form; on error returns SMESH_FAILURE and prints why.
    static int parse(std::string_view spec, CompressionOptions *out);

    /// Options from environment variable `var` (unset: `dflt`), with SMESH_COMPRESS_CHUNK and
    /// SMESH_COMPRESS_VERIFY_WRITE applied on top. A malformed value is an error (SMESH_ERROR).
    static CompressionOptions from_env(const char *var, const CompressionOptions &dflt);
    static CompressionOptions from_env(const char *var);

    std::string to_string() const;
    bool        is_lossy() const { return enabled && mode != BoundMode::Lossless; }
};

/// Compression of each kind of array a mesh folder holds.
///
/// Environment (MeshWriteOptions::from_env, used by Mesh::write(path)):
///   SMESH_COMPRESS                 all arrays of the mesh (e.g. "lossless,zstd")
///   SMESH_COMPRESS_CONNECTIVITY    i0, i1, ... (lossless only)
///   SMESH_COMPRESS_COORDS          x, y, z (e.g. "h:1e-3" for bounds relative to the element height)
///   SMESH_COMPRESS_SETS            sidesets, nodesets, edgesets (lossless only)
/// A specific variable overrides SMESH_COMPRESS for its arrays.
struct MeshWriteOptions {
    CompressionOptions connectivity;
    CompressionOptions coordinates;
    CompressionOptions sets;      // and node_mapping
    bool               any() const { return connectivity.enabled || coordinates.enabled || sets.enabled; }

    static MeshWriteOptions all(const CompressionOptions &o);
    static MeshWriteOptions from_env();
    /// SMESH_FAILURE (with a message) for requests that cannot work, e.g. lossy connectivity.
    int validate() const;
};

/// ".smz" if the options compress, else "".
static inline const char *compressed_suffix(const CompressionOptions *o) { return (o && o->enabled) ? ".smz" : ""; }

}  // namespace compression
}  // namespace smesh

#endif  // SMESH_COMPRESSION_OPTIONS_HPP
