#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <unistd.h>

using namespace smesh;

namespace {

struct MeshStats {
    ptrdiff_t n_elements{0};
    ptrdiff_t n_nodes{0};
    ptrdiff_t n_elements_out{0};
    ptrdiff_t n_nodes_out{0};
    long long bytes{0};
};

int omp_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

int slurm_nodes() {
    const char *v = std::getenv("SLURM_JOB_NUM_NODES");
    return v ? std::atoi(v) : 0;
}

const char *job_id() {
    const char *v = std::getenv("SLURM_JOB_ID");
    return v ? v : "";
}

void hostname_buf(char *buf, const size_t n) {
    if (gethostname(buf, n) != 0) {
        std::snprintf(buf, n, "unknown");
    }
    buf[n - 1] = 0;
}

MeshStats stats_of(const std::shared_ptr<Mesh> &mesh) {
    MeshStats s;
    if (!mesh) {
        return s;
    }
    if (mesh->is_distributed()) {
        auto d         = mesh->distributed();
        s.n_elements   = d->n_elements_global();
        s.n_nodes      = d->n_nodes_global();
    } else {
        s.n_elements = mesh->n_elements();
        s.n_nodes    = mesh->n_nodes();
    }
    s.n_elements_out = s.n_elements;
    s.n_nodes_out    = s.n_nodes;
    const int nxe    = elem_num_nodes(mesh->element_type(0));
    const int sdim   = mesh->spatial_dimension();
    s.bytes = (long long)s.n_nodes * (long long)sdim * (long long)sizeof(geom_t) +
              (long long)s.n_elements * (long long)nxe * (long long)sizeof(idx_t);
    return s;
}

double max_time(const std::shared_ptr<Communicator> &comm, double dt) {
#ifdef SMESH_ENABLE_MPI
    if (comm && comm->size() > 1) {
        double tmax = 0;
        MPI_Allreduce(&dt, &tmax, 1, MPI_DOUBLE, MPI_MAX, comm->get());
        return tmax;
    }
#else
    SMESH_UNUSED(comm);
#endif
    return dt;
}

double timed_barrier(const std::shared_ptr<Communicator> &comm) {
    comm->barrier();
    const double t0 = time_seconds();
    return t0;
}

double elapsed_max(const std::shared_ptr<Communicator> &comm, const double t0) {
    const double dt = time_seconds() - t0;
    comm->barrier();
    return max_time(comm, dt);
}

std::shared_ptr<Mesh> make_cube(const std::shared_ptr<Communicator> &comm,
                                const enum ElemType                  et,
                                const ptrdiff_t                      n) {
    return Mesh::create_cube(comm, et, n, n, n, 0, 0, 0, 1, 1, 1);
}

int ensure_csv_header(const Path &csv) {
    const bool write_header = !csv.exists();
    std::ofstream os(csv.to_string(), std::ios::app);
    if (!os.good()) {
        SMESH_ERROR("smesh_bench: cannot open CSV %s\n", csv.c_str());
        return SMESH_FAILURE;
    }
    if (write_header) {
        os << "kernel,elem_type,nx,n_elements,n_nodes,n_elements_out,n_nodes_out,"
              "mpi_ranks,slurm_nodes,omp_threads,repeat,time_s,melems_s,bytes,gib_s,"
              "job_id,host\n";
    }
    return SMESH_SUCCESS;
}

int append_row(const Path                             &csv,
               const char                             *kernel,
               const enum ElemType                     et,
               const ptrdiff_t                         nx,
               const MeshStats                        &st,
               const int                               mpi_ranks,
               const int                               repeat,
               const double                            time_s) {
    const double melems = time_s > 0 ? (double)st.n_elements / time_s / 1e6 : 0;
    const double gib    = time_s > 0 ? (double)st.bytes / time_s / (1024.0 * 1024.0 * 1024.0) : 0;
    char         host[256];
    hostname_buf(host, sizeof(host));
    std::ofstream os(csv.to_string(), std::ios::app);
    if (!os.good()) {
        SMESH_ERROR("smesh_bench: cannot append CSV %s\n", csv.c_str());
        return SMESH_FAILURE;
    }
    os << kernel << ',' << type_to_string(et) << ',' << nx << ',' << st.n_elements << ',' << st.n_nodes << ','
       << st.n_elements_out << ',' << st.n_nodes_out << ',' << mpi_ranks << ',' << slurm_nodes() << ','
       << omp_threads() << ',' << repeat << ',' << time_s << ',' << melems << ',' << st.bytes << ',' << gib << ','
       << job_id() << ',' << host << '\n';
    return SMESH_SUCCESS;
}

int run_repeats(const std::shared_ptr<Communicator>                 &comm,
                const Path                                          &csv,
                const char                                          *kernel,
                const enum ElemType                                  et,
                const ptrdiff_t                                      nx,
                const int                                            nrep,
                const std::function<int(int, MeshStats *, double *)> &once) {
    if (comm->rank() == 0) {
        if (ensure_csv_header(csv) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
    }
    comm->barrier();
    for (int r = 0; r < nrep; ++r) {
        MeshStats st;
        double    dt = 0;
        if (once(r, &st, &dt) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        if (comm->rank() == 0) {
            if (append_row(csv, kernel, et, nx, st, comm->size(), r, dt) != SMESH_SUCCESS) {
                return SMESH_FAILURE;
            }
            std::printf("smesh_bench kernel=%s elem=%s N=%ld ranks=%d repeat=%d time_s=%.6f melems_s=%.6f gib_s=%.6f\n",
                        kernel,
                        type_to_string(et),
                        (long)nx,
                        comm->size(),
                        r,
                        dt,
                        dt > 0 ? (double)st.n_elements / dt / 1e6 : 0,
                        dt > 0 ? (double)st.bytes / dt / (1024.0 * 1024.0 * 1024.0) : 0);
            std::fflush(stdout);
        }
    }
    return SMESH_SUCCESS;
}

}  // namespace

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("smesh_bench.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc < 5 || argc > 6) {
        if (comm->rank() == 0) {
            std::fprintf(stderr,
                         "Usage: %s <generate|read|write|promote|refine|io> <HEX8|TET4> <N> <csv> [mesh_dir]\n",
                         argv[0]);
        }
        return SMESH_FAILURE;
    }

    const std::string   kernel = argv[1];
    const enum ElemType et     = type_from_string(argv[2]);
    const ptrdiff_t     n      = (ptrdiff_t)std::atoll(argv[3]);
    const Path          csv(argv[4]);
    const Path          mesh_dir(argc == 6 ? argv[5] : "");
    const int           nrep   = Env::read<int32_t>("SMESH_BENCH_REPEAT", 3);
    const int           levels = Env::read<int32_t>("SMESH_REFINEMENT_LEVELS", 1);
    const bool          keep   = Env::read<bool>("SMESH_BENCH_KEEP", false);

    if (et != HEX8 && et != TET4) {
        if (comm->rank() == 0) {
            SMESH_ERROR("smesh_bench: elem must be HEX8 or TET4\n");
        }
        return SMESH_FAILURE;
    }
    if (n <= 0 || nrep <= 0) {
        if (comm->rank() == 0) {
            SMESH_ERROR("smesh_bench: N and SMESH_BENCH_REPEAT must be > 0\n");
        }
        return SMESH_FAILURE;
    }
    if ((kernel == "read" || kernel == "write" || kernel == "io") && mesh_dir.empty()) {
        if (comm->rank() == 0) {
            SMESH_ERROR("smesh_bench: %s requires mesh_dir\n", kernel.c_str());
        }
        return SMESH_FAILURE;
    }
    if (kernel == "promote" && et != TET4) {
        if (comm->rank() == 0) {
            SMESH_ERROR("smesh_bench: promote requires TET4 (TET4→TET10)\n");
        }
        return SMESH_FAILURE;
    }

    auto time_generate = [&](int, MeshStats *st, double *dt) -> int {
        comm->barrier();
        const double t0 = timed_barrier(comm);
        auto         mesh = make_cube(comm, et, n);
        *dt               = elapsed_max(comm, t0);
        if (!mesh) {
            return SMESH_FAILURE;
        }
        *st = stats_of(mesh);
        return SMESH_SUCCESS;
    };

    auto time_write = [&](int, MeshStats *st, double *dt) -> int {
        std::shared_ptr<Mesh> mesh;
        if (mesh_dir.exists()) {
            mesh = Mesh::create_from_file(comm, mesh_dir);
        } else {
            mesh = make_cube(comm, et, n);
        }
        if (!mesh) {
            return SMESH_FAILURE;
        }
        *st = stats_of(mesh);
        comm->barrier();
        const double t0 = timed_barrier(comm);
        const int    err = mesh->write(mesh_dir);
        *dt              = elapsed_max(comm, t0);
        return err;
    };

    auto time_read = [&](int, MeshStats *st, double *dt) -> int {
        if (!mesh_dir.exists()) {
            if (comm->rank() == 0) {
                SMESH_ERROR("smesh_bench: mesh_dir does not exist: %s\n", mesh_dir.c_str());
            }
            return SMESH_FAILURE;
        }
        comm->barrier();
        const double t0   = timed_barrier(comm);
        auto         mesh = Mesh::create_from_file(comm, mesh_dir);
        *dt               = elapsed_max(comm, t0);
        if (!mesh) {
            return SMESH_FAILURE;
        }
        *st = stats_of(mesh);
        return SMESH_SUCCESS;
    };

    auto time_promote = [&](int, MeshStats *st, double *dt) -> int {
        auto coarse = make_cube(comm, TET4, n);
        if (!coarse) {
            return SMESH_FAILURE;
        }
        comm->barrier();
        const double t0     = timed_barrier(comm);
        auto         fine   = promote_to(TET10, coarse);
        *dt                 = elapsed_max(comm, t0);
        if (!fine) {
            return SMESH_FAILURE;
        }
        *st                 = stats_of(coarse);
        const MeshStats out = stats_of(fine);
        st->n_elements_out  = out.n_elements;
        st->n_nodes_out     = out.n_nodes;
        st->bytes           = out.bytes;
        if (keep && !mesh_dir.empty()) {
            fine->write(mesh_dir);
        }
        return SMESH_SUCCESS;
    };

    auto time_refine = [&](int, MeshStats *st, double *dt) -> int {
        auto coarse = make_cube(comm, et, n);
        if (!coarse) {
            return SMESH_FAILURE;
        }
        comm->barrier();
        const double t0     = timed_barrier(comm);
        auto         fine   = refine(coarse, levels);
        *dt                 = elapsed_max(comm, t0);
        if (!fine) {
            return SMESH_FAILURE;
        }
        *st                 = stats_of(coarse);
        const MeshStats out = stats_of(fine);
        st->n_elements_out  = out.n_elements;
        st->n_nodes_out     = out.n_nodes;
        st->bytes           = out.bytes;
        if (keep && !mesh_dir.empty()) {
            fine->write(mesh_dir);
        }
        return SMESH_SUCCESS;
    };

    int err = SMESH_FAILURE;
    if (kernel == "generate") {
        err = run_repeats(comm, csv, "generate", et, n, nrep, time_generate);
    } else if (kernel == "write") {
        err = run_repeats(comm, csv, "write", et, n, nrep, time_write);
    } else if (kernel == "read") {
        err = run_repeats(comm, csv, "read", et, n, nrep, time_read);
    } else if (kernel == "promote") {
        err = run_repeats(comm, csv, "promote", TET4, n, nrep, time_promote);
    } else if (kernel == "refine") {
        err = run_repeats(comm, csv, "refine", et, n, nrep, time_refine);
    } else if (kernel == "io") {
        err = run_repeats(comm, csv, "write", et, n, nrep, time_write);
        if (err == SMESH_SUCCESS) {
            err = run_repeats(comm, csv, "read", et, n, nrep, time_read);
        }
    } else {
        if (comm->rank() == 0) {
            SMESH_ERROR("smesh_bench: unknown kernel %s\n", kernel.c_str());
        }
        return SMESH_FAILURE;
    }
    return err;
}
