# Compressed arrays (`.smz`)

smesh can store any array of a mesh folder, or any field, compressed: losslessly, or with a
guaranteed error bound. Serial and MPI runs read and write the same files. A rank count does not
change a file, and raw folders keep working.

## Quick start

```bash
# Compress every array of a mesh folder written by Mesh::write (environment, off by default)
SMESH_COMPRESS=lossless,zstd ./my_app ...
# Coordinates with an error below 1e-3 of the local element height, the rest lossless
SMESH_COMPRESS=lossless SMESH_COMPRESS_COORDS=h:1e-3 ./my_app ...
# Fields written through smesh::Output: error below 1e-4 of each field's range
SMESH_OUTPUT_COMPRESS=range:1e-4,zstd ./my_app ...

# Command line
smz mesh in_folder out_folder lossless,zstd h:1e-3,zstd   # MPI-parallel (mpiexec -n P smz mesh ...)
smz compress u.float64 range:1e-5,zstd                    # -> u.float64.smz
smz info u.float64.smz                                    # ratio, bound, achieved error, chunks
smz decompress u.float64.smz                              # -> u.float64 (bit-exact if lossless)
smz verify u.float64.smz                                  # decode everything, check checksums
```

```cpp
#include "smesh_compression_options.hpp"
using namespace smesh::compression;

CompressionOptions coords;
CompressionOptions::parse("h:1e-3,zstd", &coords);
MeshWriteOptions opts = MeshWriteOptions::all(CompressionOptions::lossless());
opts.coordinates = coords;
mesh->write(path, opts);                      // Mesh::create_from_file reads it as usual

output->set_compression(range_opts);          // smesh::Output
buffer->to_file(path, opts);                  // Buffer<T>: writes path + ".smz"
smz_write<f64>(path, x, n, opts);             // single arrays; smz_read_convert<T>, smz_read_range<T>
```

Readers need no option: `array_read_convert_from_extension`, `Buffer::from_file`, the mesh and
set readers, serial and distributed, recognise `<stem>.<dtype>.smz` and decode it. A folder may
mix raw and compressed files. If both exist for one array, the one named in `meta.yaml` is read,
otherwise the raw one.

Python: the converters in `python/smesh` (`common/raw_io.py`, `raw_to_db.py`) read `.smz` files
too; `common/smz.py` parses headers in pure Python and decodes through the `smz` tool (on `PATH`,
or `SMESH_SMZ=/path/to/smz`). `Mesh::write_with_xdmf` always writes raw files, since XDMF points
readers at raw binary arrays.

## Specs

`off | lossless | abs:E | range:E | pwrel:E | h:E` followed by options
`[,zstd[:LEVEL]] [,lz4] [,shuffle] [,chunk:N] [,nochecksum] [,verify]`.

| mode | guarantee, for every finite value | use |
|---|---|---|
| `lossless` | bit-exact (also NaN payloads, -0, denormals) | connectivity, sets, anything exact |
| `abs:E` | \|x' - x\| <= E | fields with a known tolerance |
| `range:E` | \|x' - x\| <= E (max - min), over the whole array (all ranks) | simulation output |
| `pwrel:E` | \|x' - x\| <= E \|x\|, same sign, zeros exact (0 < E < 1) | fields spanning many magnitudes |
| `h:E` | \|x' - x\| <= E h_n per component, h_n = smallest element height at node n (E < 0.25) | mesh coordinates |

NaN and Inf are always stored exactly. Lossy modes on integer arrays are rejected; connectivity
and sets are always lossless.

Why `h` for coordinates: a lossy element geometry must not break the operator. Quantising
per-element geometric factors (e.g. sfem's FFF) independently breaks the patch test; moving the
nodes instead gives the exact operator of a slightly perturbed mesh. The bound must be relative
to the element height, not its edges: the geometric factors change by about |dx| / height, and
graded or flat elements have edges far longer than their heights. With `h:E`, every simplex and
every hex corner tet keeps its orientation, and its volume changes by at most 4 sqrt(3) E
(measured: below 1.6 E on strongly graded meshes, `smz_node_height_test`).

Options: `zstd` / `lz4` add a general-purpose stage after smesh's codecs (kept per chunk only
when it saves 5%; reading such a file needs a build with that library). `shuffle` stores
lossless floats as byte planes, which zstd compresses well. `chunk:N` sets the values per
chunk (default 65536, or `SMESH_COMPRESS_CHUNK`). `verify` decodes every chunk after encoding it
(`SMESH_COMPRESS_VERIFY_WRITE=1` for all writes).

Environment: `SMESH_COMPRESS` (every mesh array), `SMESH_COMPRESS_CONNECTIVITY`,
`SMESH_COMPRESS_COORDS`, `SMESH_COMPRESS_SETS`, `SMESH_OUTPUT_COMPRESS` (Output fields),
`SMESH_COMPRESS_REPORT=1` (ratio, error and time per file written),
`SMESH_COMPRESS_VERIFY_READ=0` (skip the checksums when reading).

## How it works

* Values map to order-preserving unsigned words, so one integer codec serves every type:
  blocks of 128 values as the smallest of frame-of-reference, delta (zig-zag), and their
  patched forms (a few outliers stored apart), bit-packed.
* Lossy floats are quantised on power-of-two grids, x' = origin + q 2^L, and coded as integer
  residuals: one step for `abs`/`range`, a step per value for `pwrel`/`h`. The decoded value
  is exact up to one rounding, so CPUs, GPUs and compilers decode the same bits. The encoder
  checks every value with the decoder's own function and stores any value that misses its
  bound exactly; a chunk with too many such values is written losslessly.
* Arrays are cut into chunks; a chunk depends only on its values and on parameters fixed for
  the whole array (reduced over ranks). Chunks are encoded in parallel (OpenMP) and joined in
  order, so files do not depend on the number of threads or ranks.
* MPI: writers redistribute values to chunk-aligned slices and write with collective MPI-IO at
  offsets from an exscan of the encoded sizes; readers take matrix.io's uniform split and read
  only the chunks they need. Node heights for `h` are the minimum over every rank's elements.

## File format (version 1.0)

```
[header 128 B][index: 32 B per chunk][payload: the chunks, in order]     (little-endian)
```

Header: magic `\x89SMZ\r\n\x1a\n`, version, value type, bound mode, planned codec and byte stage,
counts (values, chunk size, chunks, payload bytes), the bound requested and used, the grid
(origin, levels), the finite range, the achieved maximum error and error/bound, and an XXH64
of the header. Index entry: offset, XXH64 of the stored bytes, stored and pre-stage sizes,
number of values stored exactly, codec and byte stage. The full layout is in
`src/compression/smesh_smz_format.hpp`. Readers check the header hash, the version, that the
chunks tile the payload and the file size, and every chunk checksum; a damaged file is
reported, never decoded silently.

## Code

| layer | files |
|---|---|
| kernels (header-only, no IO/MPI, `SMESH_HOST_DEVICE` decoders) | `src/compression/smesh_codec_{bits,int,float,quant,hash}.hpp` |
| zstd / lz4 stage | `src/compression/smesh_byte_codec.*` |
| chunked arrays | `src/compression/smesh_chunked_codec.*` |
| format, serial IO, options | `src/compression/smesh_smz_{format,io}.*`, `smesh_compression_options.*` |
| MPI IO | `src/distributed/smesh_smz_mpi.*` |
| element heights | `src/mesh/geometry/smesh_node_height.*` |
| tools | `src/drivers/io/smz.exe.cpp`, `src/benchmark/drivers/smz_bench.exe.cpp` |
| tests | `src/compression/tests/smz_*`, `src/distributed/tests/smz_*_mpi_test.cpp` (`ctest -L smz`, `ctest -L mpi`) |

The kernels are meant to be reused for operator data held in memory (sfem's packed operators):
`dec_packed<T, false>` and `dec_quant_var<T, false>` decode without bounds checks, and
`enc_quant_var` takes a gather functor to encode a pack's node list from global arrays with the
global grid, so shared nodes decode identically in every pack.

## Build

`SMESH_ENABLE_ZSTD` (default ON, used if found) and `SMESH_ENABLE_LZ4` (default OFF, used if
found); set `ZSTD_ROOT` / `LZ4_ROOT` if they are not in a standard place. Without them, everything
but the byte stages works.

## Measured (Mac M1 Max, 8 threads; small arrays from page cache)

| data | spec | ratio | encode / decode GB/s |
|---|---|---|---|
| float64 field (smooth) | `lossless,shuffle,zstd` | 3.9x | 2.8 / 6.0 |
| float64 field | `range:1e-5,zstd` | 28.7x | 2.7 / 15.6 |
| int32 connectivity, structured HEX8 | `lossless` / `lossless,zstd` | 23x / 234x | 3.5 / 17.6 |
| HEX8 64^3 folder (smesh cube) | `lossless` / `lossless,zstd` | 3.3x / 34x | |
| TET4 40^3 folder | `lossless` / `lossless,zstd` | 2.9x / 6.1x | |
| graded HEX8 / TET4 coordinates (h ratio 1e5) | `h:1e-3` | 1.9x / 1.7x | |

Decode of the packed integer codec: 11.5 GB/s per core (`smz_lossless_kernels_test`).
