// smz: inspect, create and unpack compressed array files (.smz), see src/compression.
//
//   smz info <file.smz>                         header, chunk codecs, ratio, error, escapes
//   smz compress <file.<dtype>> [spec] [out]    raw array -> <file.<dtype>>.smz (spec: lossless,
//                                               abs:E, range:E, pwrel:E, [,zstd][,chunk:N]...)
//   smz decompress <file.<dtype>.smz> [out]     -> the raw array, bit-exact for lossless files
//   smz verify <file.smz>                       decode everything, check every chunk checksum
//   smz mesh <in_folder> <out_folder> [spec] [coords_spec]
//                                               rewrite a mesh folder compressed (MPI-parallel;
//                                               e.g. "lossless,zstd" "h:1e-3,zstd"); "off" writes raw

#include "smesh_alloc.hpp"
#include "smesh_base.hpp"
#include "smesh_chunked_codec.hpp"
#include "smesh_communicator.hpp"
#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_read.hpp"
#include "smesh_read.impl.hpp"
#include "smesh_smz_format.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_write.hpp"
#include "smesh_write.impl.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace smesh;
using namespace smesh::compression;

static int usage() {
    fprintf(stderr,
            "usage: smz info <file.smz>\n"
            "       smz compress <file.<dtype>> [spec] [out.smz]\n"
            "       smz decompress <file.<dtype>.smz> [out]\n"
            "       smz verify <file.smz>\n"
            "       smz mesh <in_folder> <out_folder> [spec] [coords_spec]\n"
            "spec: off|lossless|abs:E|range:E|pwrel:E [,zstd[:L]] [,lz4] [,shuffle] [,chunk:N] [,verify]\n");
    return SMESH_FAILURE;
}

static int cmd_info(const Path &path) {
    SmzInfo                 info;
    std::vector<ChunkEntry> index;
    if (smz_read_index(path, &info, &index) != SMESH_SUCCESS) return SMESH_FAILURE;
    const ArrayPlan &p   = info.plan;
    const double     raw = double(p.n_values) * double(dtype_bytes(p.dtype));
    std::map<std::string, u64> codecs;
    u64                        esc = 0;
    for (const auto &e : index) {
        std::string k = to_string(e.codec);
        if (e.byte_stage != ByteStage::None) k += std::string("+") + to_string(e.byte_stage);
        codecs[k]++;
        esc += e.n_escaped;
    }
    printf("file        %s\n", path.c_str());
    printf("format      smz %u.%u\n", info.version_major, info.version_minor);
    printf("values      %llu x %s (%.3f MB raw)\n", (unsigned long long)p.n_values, to_string(p.dtype), raw / 1e6);
    printf("size        %.3f MB (%.2fx), %llu chunks of %llu values\n", double(info.file_bytes()) / 1e6,
           info.file_bytes() ? raw / double(info.file_bytes()) : 0.0, (unsigned long long)p.n_chunks,
           (unsigned long long)p.chunk_values);
    printf("bound       %s", to_string(p.requested));
    if (p.requested != BoundMode::Lossless) printf(" %g", p.quant.eb);
    if (p.codec == Codec::QuantFixed) printf(" (absolute %g, step 2^%d)", p.quant.eb_abs, p.quant.L);
    if (p.codec == Codec::QuantVar) printf(" (finest step 2^%d)", p.quant.L_min);
    printf("\nachieved    max |x'-x| %.6g, %.4f of the bound, %llu values stored exactly\n", info.max_abs_err,
           info.max_ratio, (unsigned long long)esc);
    printf("range       [%.17g, %.17g]\n", info.vmin, info.vmax);
    printf("checksums   %s%s\n", p.checksums ? "xxh64" : "none", (info.flags & 1u) ? ", verified on write" : "");
    printf("chunks     ");
    for (const auto &kv : codecs) printf(" %s:%llu", kv.first.c_str(), (unsigned long long)kv.second);
    printf("\n");
    return SMESH_SUCCESS;
}

template <class T>
static int compress_as(const Path &in, const Path &out, const CompressionOptions &opts) {
    T        *x = nullptr;
    ptrdiff_t n = 0;
    if (array_read_convert_from_extension<T>(in, &x, &n) != SMESH_SUCCESS) return SMESH_FAILURE;
    SmzReport r;
    const int ret = smz_write<T>(out, x, n, opts, nullptr, &r);
    SMESH_FREE(x);
    if (ret == SMESH_SUCCESS) {
        printf("%s -> %s: %.3f -> %.3f MB (%.2fx), max err %.6g (%.4f of the bound), %llu exact, encode %.3f s\n",
               in.c_str(), out.c_str(), double(r.raw_bytes) / 1e6, double(r.file_bytes) / 1e6, r.ratio(),
               r.max_abs_err, r.max_ratio, (unsigned long long)r.n_escaped, r.encode_s);
    }
    return ret;
}

static int cmd_compress(const Path &in, const std::string &spec, const Path &out) {
    CompressionOptions opts;
    if (CompressionOptions::parse(spec, &opts) != SMESH_SUCCESS) return SMESH_FAILURE;
    if (!opts.enabled) opts = CompressionOptions::lossless();
    const std::string ext = typed_extension(in);
    if (ext == "float16") return compress_as<f16>(in, out, opts);
    if (ext == "float32") return compress_as<f32>(in, out, opts);
    if (ext == "float64") return compress_as<f64>(in, out, opts);
    if (ext == "int16") return compress_as<i16>(in, out, opts);
    if (ext == "int32") return compress_as<i32>(in, out, opts);
    if (ext == "int64") return compress_as<i64>(in, out, opts);
    fprintf(stderr, "smz: cannot tell the value type of %s (expected a .float32/.int32/... name)\n", in.c_str());
    return SMESH_FAILURE;
}

template <class T>
static int decompress_as(const Path &in, const Path &out) {
    T        *x = nullptr;
    ptrdiff_t n = 0;
    if (smz_read_convert<T>(in, &x, &n) != SMESH_SUCCESS) return SMESH_FAILURE;
    const int ret = array_write<T>(out, x, n);
    SMESH_FREE(x);
    return ret;
}

static int cmd_decompress(const Path &in, const Path &out) {
    SmzInfo info;
    if (smz_read_info(in, &info) != SMESH_SUCCESS) return SMESH_FAILURE;
    switch (info.plan.dtype) {
        case DType::I8: return decompress_as<i8>(in, out);
        case DType::U8: return decompress_as<u8>(in, out);
        case DType::I16: return decompress_as<i16>(in, out);
        case DType::U16: return decompress_as<u16>(in, out);
        case DType::I32: return decompress_as<i32>(in, out);
        case DType::U32: return decompress_as<u32>(in, out);
        case DType::I64: return decompress_as<i64>(in, out);
        case DType::U64: return decompress_as<u64>(in, out);
        case DType::F16: return decompress_as<f16>(in, out);
        case DType::F32: return decompress_as<f32>(in, out);
        case DType::F64: return decompress_as<f64>(in, out);
    }
    return SMESH_FAILURE;
}

static int cmd_verify(const Path &in) {
    const auto t0 = std::chrono::steady_clock::now();
    f64       *x  = nullptr;
    ptrdiff_t  n  = 0;
    if (smz_read_convert<f64>(in, &x, &n) != SMESH_SUCCESS) {
        printf("%s: CORRUPT\n", in.c_str());
        return SMESH_FAILURE;
    }
    SMESH_FREE(x);
    printf("%s: OK, %td values, decoded in %.3f s\n", in.c_str(), n,
           std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return SMESH_SUCCESS;
}

static int cmd_mesh(const Path &in, const Path &out, const std::string &spec, const std::string &coords) {
    CompressionOptions a, b;
    if (CompressionOptions::parse(spec, &a) != SMESH_SUCCESS) return SMESH_FAILURE;
    b = a;
    if (!coords.empty() && CompressionOptions::parse(coords, &b) != SMESH_SUCCESS) return SMESH_FAILURE;
    MeshWriteOptions o = MeshWriteOptions::all(a);
    o.coordinates      = b;
    auto comm          = Communicator::world();
    auto mesh          = Mesh::create_from_file(comm, in);
    if (!mesh) return SMESH_FAILURE;
    return mesh->write(out, o);
}

int main(int argc, char *argv[]) {
    smesh::Context context(argc, argv);
    if (argc < 3) return usage();
    const std::string cmd = argv[1];
    const Path        in(argv[2]);
    int               ret = SMESH_FAILURE;
    if (cmd == "info") {
        ret = cmd_info(in);
    } else if (cmd == "compress") {
        const std::string spec = argc > 3 ? argv[3] : "lossless";
        const Path        out  = argc > 4 ? Path(argv[4]) : Path(in.to_string() + ".smz");
        ret                    = cmd_compress(in, spec, out);
    } else if (cmd == "decompress") {
        Path out = argc > 3 ? Path(argv[3]) : in.without_extension();
        ret      = cmd_decompress(in, out);
    } else if (cmd == "verify") {
        ret = cmd_verify(in);
    } else if (cmd == "mesh" && argc > 3) {
        ret = cmd_mesh(in, Path(argv[3]), argc > 4 ? argv[4] : "lossless", argc > 5 ? argv[5] : "");
    } else {
        return usage();
    }
    return ret;
}
