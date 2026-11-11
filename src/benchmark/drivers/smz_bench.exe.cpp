// smz_bench: compression ratio and speed of smesh's compressed IO (src/compression).
//
//   smz_bench array <file.<dtype>> <scratch_dir> [spec ...]
//       serial (OpenMP): ratio, encode and decode GB/s, achieved error, per spec
//   mpiexec -n P smz_bench mesh <mesh_folder> <scratch_dir> [spec ...]
//       Mesh::write and Mesh::create_from_file of the whole mesh folder, raw and per spec
//   mpiexec -n P smz_bench read <folder> [folder ...]
//       Mesh::create_from_file only, e.g. of folders written by `smz mesh` in an earlier job, so
//       that the files are not in this node's page cache (the first read is reported)
//
// A spec is a compression spec for every array ("lossless,zstd"), or "<connectivity>|<coords>"
// ("lossless,zstd|h:1e-3,zstd") for meshes. Without specs a default sweep runs. Rows are also
// appended as CSV to $SMESH_BENCH_CSV if set. Repetitions: $SMESH_BENCH_REPEAT (default 3, the
// best time is reported).

#include "smesh_alloc.hpp"
#include "smesh_communicator.hpp"
#include "smesh_compression_options.hpp"
#include "smesh_context.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_mesh.hpp"
#include "smesh_read.hpp"
#include "smesh_read.impl.hpp"
#include "smesh_smz_io.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace smesh;
using namespace smesh::compression;
namespace fs = std::filesystem;

namespace {
using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

int repeats() {
    const char *v = std::getenv("SMESH_BENCH_REPEAT");
    return v && std::atoi(v) > 0 ? std::atoi(v) : 3;
}

int threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

u64 tree_bytes(const fs::path &p) {
    u64 b = 0;
    if (fs::is_regular_file(p)) return fs::file_size(p);
    for (const auto &e : fs::recursive_directory_iterator(p))
        if (e.is_regular_file()) b += e.file_size();
    return b;
}

FILE *csv_file(const char *header) {
    const char *path = std::getenv("SMESH_BENCH_CSV");
    if (!path || !*path) return nullptr;
    const bool fresh = !fs::exists(path) || fs::file_size(path) == 0;
    FILE      *f     = fopen(path, "a");
    if (f && fresh) fprintf(f, "%s\n", header);
    return f;
}

// ------------------------------------------------------------------------------ array
template <class T>
int bench_array(const Path &in, const fs::path &scratch, const std::vector<std::string> &specs) {
    T        *x = nullptr;
    ptrdiff_t n = 0;
    if (array_read_convert_from_extension<T>(in, &x, &n) != SMESH_SUCCESS) return SMESH_FAILURE;
    const double raw = double(n) * sizeof(T);
    FILE *csv = csv_file("mode,file,spec,threads,n_values,raw_MB,file_MB,ratio,encode_GBps,decode_GBps,max_err,"
                         "err_over_bound,escaped");
    printf("%-34s %9s %9s %8s %10s %10s %11s %9s\n", "spec", "raw MB", "file MB", "ratio", "enc GB/s", "dec GB/s",
           "max err", "err/bound");
    const Path out((scratch / ("bench." + std::string(to_string(DTypeOf<T>::value)) + ".smz")).string());
    for (const auto &spec : specs) {
        CompressionOptions o;
        if (CompressionOptions::parse(spec, &o) != SMESH_SUCCESS) continue;
        if (!o.enabled) o = CompressionOptions::lossless();
        double    enc = 1e30, dec = 1e30;
        SmzReport r;
        for (int k = 0; k < repeats(); k++) {
            auto t0 = Clock::now();
            if (smz_write<T>(out, x, n, o, nullptr, &r) != SMESH_SUCCESS) {
                SMESH_FREE(x);
                return SMESH_FAILURE;
            }
            enc = std::min(enc, r.encode_s);
            (void)t0;
            T        *y  = nullptr;
            ptrdiff_t ny = 0;
            t0           = Clock::now();
            if (smz_read_convert<T>(out, &y, &ny) != SMESH_SUCCESS) {
                SMESH_FREE(x);
                return SMESH_FAILURE;
            }
            dec = std::min(dec, seconds_since(t0));
            SMESH_FREE(y);
        }
        const double file = double(fs::file_size(out.to_string()));
        printf("%-34s %9.2f %9.2f %8.2f %10.2f %10.2f %11.3g %9.3f\n", spec.c_str(), raw / 1e6, file / 1e6, raw / file,
               raw / enc / 1e9, raw / dec / 1e9, r.max_abs_err, r.max_ratio);
        if (csv)
            fprintf(csv, "array,%s,\"%s\",%d,%td,%.3f,%.3f,%.4f,%.3f,%.3f,%.6g,%.4f,%llu\n", in.c_str(), spec.c_str(),
                    threads(), n, raw / 1e6, file / 1e6, raw / file, raw / enc / 1e9, raw / dec / 1e9, r.max_abs_err,
                    r.max_ratio, (unsigned long long)r.n_escaped);
    }
    if (csv) fclose(csv);
    fs::remove(out.to_string());
    SMESH_FREE(x);
    return SMESH_SUCCESS;
}

int cmd_array(const Path &in, const fs::path &scratch, std::vector<std::string> specs) {
    const std::string ext = typed_extension(in);
    const bool        flt = ext.rfind("float", 0) == 0;
    if (specs.empty()) {
        specs = {"lossless", "lossless,zstd", "lossless,shuffle,zstd"};
        if (flt) specs.insert(specs.end(), {"range:1e-3", "range:1e-5", "range:1e-5,zstd", "pwrel:1e-3"});
    }
    if (ext == "float32") return bench_array<f32>(in, scratch, specs);
    if (ext == "float64") return bench_array<f64>(in, scratch, specs);
    if (ext == "int32") return bench_array<i32>(in, scratch, specs);
    if (ext == "int64") return bench_array<i64>(in, scratch, specs);
    if (ext == "int16") return bench_array<i16>(in, scratch, specs);
    fprintf(stderr, "smz_bench: unsupported array type '%s'\n", ext.c_str());
    return SMESH_FAILURE;
}

// ------------------------------------------------------------------------------ mesh
int cmd_mesh(const Path &in, const fs::path &scratch, std::vector<std::string> specs) {
    auto comm = Communicator::world();
    if (specs.empty())
        specs = {"off", "lossless", "lossless,zstd", "lossless|h:1e-2", "lossless|h:1e-3", "lossless,zstd|h:1e-3,zstd"};
    FILE *csv = comm->rank() == 0 ? csv_file("mode,mesh,spec,ranks,threads,n_elements,n_nodes,folder_MB,ratio,write_s,"
                                             "read_s,write_GBps_raw,read_GBps_raw")
                                  : nullptr;
    double raw_bytes = 0;
    {
        auto t0   = Clock::now();
        auto mesh = Mesh::create_from_file(comm, in);
        if (!mesh) return SMESH_FAILURE;
        if (comm->rank() == 0) printf("read %s in %.3f s\n", in.c_str(), seconds_since(t0));
    }
    auto mesh = Mesh::create_from_file(comm, in);
    if (comm->rank() == 0)
        printf("%-30s %10s %8s %9s %9s %12s %12s\n", "spec", "folder MB", "ratio", "write s", "read s", "write GB/s*",
               "read GB/s*");
    for (size_t k = 0; k < specs.size(); k++) {
        const std::string &spec = specs[k];
        MeshWriteOptions   o;
        const size_t       bar = spec.find('|');
        CompressionOptions a, b;
        if (CompressionOptions::parse(spec.substr(0, bar), &a) != SMESH_SUCCESS) continue;
        b = a;
        if (bar != std::string::npos && CompressionOptions::parse(spec.substr(bar + 1), &b) != SMESH_SUCCESS) continue;
        o              = MeshWriteOptions::all(a);
        o.coordinates  = b;
        const fs::path out = scratch / ("mesh_" + std::to_string(k));
        if (comm->rank() == 0) fs::remove_all(out);
        comm->barrier();
        double wbest = 1e30, rbest = 1e30;
        for (int r = 0; r < repeats(); r++) {
            comm->barrier();
            auto t0 = Clock::now();
            if (mesh->write(Path(out.string()), o) != SMESH_SUCCESS) return SMESH_FAILURE;
            comm->barrier();
            wbest = std::min(wbest, seconds_since(t0));
            t0    = Clock::now();
            auto back = Mesh::create_from_file(comm, Path(out.string()));
            comm->barrier();
            rbest = std::min(rbest, seconds_since(t0));
            if (!back) return SMESH_FAILURE;
        }
        if (comm->rank() == 0) {
            const double bytes = double(tree_bytes(out));
            if (k == 0 || raw_bytes == 0) raw_bytes = spec == "off" ? bytes : raw_bytes;
            const double ratio = raw_bytes > 0 ? raw_bytes / bytes : 0;
            printf("%-30s %10.2f %8.2f %9.3f %9.3f %12.2f %12.2f\n", spec.c_str(), bytes / 1e6, ratio, wbest, rbest,
                   raw_bytes / wbest / 1e9, raw_bytes / rbest / 1e9);
            if (csv)
                fprintf(csv, "mesh,%s,\"%s\",%d,%d,%td,%td,%.3f,%.4f,%.4f,%.4f,%.3f,%.3f\n", in.c_str(), spec.c_str(),
                        comm->size(), threads(), mesh->n_elements(), mesh->n_nodes(), bytes / 1e6, ratio, wbest,
                        rbest, raw_bytes / wbest / 1e9, raw_bytes / rbest / 1e9);
            fs::remove_all(out);
        }
    }
    if (comm->rank() == 0) printf("* raw folder bytes / time\n");
    if (csv) fclose(csv);
    return SMESH_SUCCESS;
}
int cmd_read(const std::vector<std::string> &folders) {
    auto  comm = Communicator::world();
    FILE *csv  = comm->rank() == 0 ? csv_file("mode,folder,ranks,threads,n_elements,n_nodes,folder_MB,first_read_s,"
                                              "best_read_s,GBps_first")
                                   : nullptr;
    if (comm->rank() == 0)
        printf("%-50s %10s %12s %11s %10s\n", "folder", "MB", "1st read s", "best read s", "GB/s 1st");
    for (const auto &f : folders) {
        double first = 0, best = 1e30;
        std::shared_ptr<Mesh> m;
        for (int r = 0; r < repeats(); r++) {
            comm->barrier();
            const auto t0 = Clock::now();
            m             = Mesh::create_from_file(comm, Path(f));
            comm->barrier();
            const double t = seconds_since(t0);
            if (!m) return SMESH_FAILURE;
            if (r == 0) first = t;
            best = std::min(best, t);
        }
        if (comm->rank() == 0) {
            const double bytes = double(tree_bytes(f));
            printf("%-50s %10.2f %12.3f %11.3f %10.2f\n", f.c_str(), bytes / 1e6, first, best, bytes / first / 1e9);
            if (csv)
                fprintf(csv, "read,%s,%d,%d,%td,%td,%.3f,%.4f,%.4f,%.3f\n", f.c_str(), comm->size(), threads(),
                        m->n_elements(), m->n_nodes(), bytes / 1e6, first, best, bytes / first / 1e9);
        }
    }
    if (csv) fclose(csv);
    return SMESH_SUCCESS;
}
}  // namespace

int main(int argc, char *argv[]) {
    smesh::Context context(argc, argv);
    if (argc >= 3 && std::string(argv[1]) == "read") {
        return cmd_read(std::vector<std::string>(argv + 2, argv + argc));
    }
    if (argc < 4) {
        fprintf(stderr,
                "usage: smz_bench array <file.<dtype>> <scratch_dir> [spec ...]\n"
                "       smz_bench mesh <mesh_folder> <scratch_dir> [spec|conn_spec|coords_spec ...]\n"
                "       smz_bench read <folder> [folder ...]\n");
        return SMESH_FAILURE;
    }
    const std::string        mode = argv[1];
    std::vector<std::string> specs;
    for (int i = 4; i < argc; i++) specs.push_back(argv[i]);
    fs::create_directories(argv[3]);
    if (mode == "array") return cmd_array(Path(argv[2]), argv[3], specs);
    if (mode == "mesh") return cmd_mesh(Path(argv[2]), argv[3], specs);
    fprintf(stderr, "smz_bench: unknown mode %s\n", mode.c_str());
    return SMESH_FAILURE;
}
