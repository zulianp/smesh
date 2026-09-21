#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("adapt_refine.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 3) {
        if (!comm->rank()) {
            fprintf(stderr, "Usage: %s <mesh_folder> <output_folder>\n", argv[0]);
        }
        return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_from_file(comm, Path(argv[1]));
    if (!mesh) {
        return SMESH_FAILURE;
    }

    AdaptRefineOptions opt;
    opt.cells_per_radius    = Env::read<geom_t>("SMESH_CELLS_PER_RADIUS", static_cast<geom_t>(8));
    opt.sharp_cos_threshold = Env::read<geom_t>("SMESH_SHARP_EDGES_THRESHOLD", static_cast<geom_t>(0.15));
    opt.max_levels          = Env::read<int>("SMESH_ADAPT_LEVELS", 8);
    opt.smooth_iters        = Env::read<int>("SMESH_SMOOTH_ITERS", 10);
    opt.use_parametrization = true;

    auto fine = adapt_refine(mesh, opt);
    if (!fine) {
        return SMESH_FAILURE;
    }
    return fine->write(Path(argv[2]));
}
