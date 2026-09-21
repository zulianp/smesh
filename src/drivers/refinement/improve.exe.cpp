#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("improve.exe");
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

    ImproveOptions opt;
    opt.q_min               = Env::read<geom_t>("SMESH_Q_MIN", static_cast<geom_t>(0.3));
    opt.max_abs_dev         = Env::read<geom_t>("SMESH_MAX_ABS_DEV", static_cast<geom_t>(0));
    opt.max_normal_dev      = Env::read<geom_t>("SMESH_MAX_NORMAL_DEV", static_cast<geom_t>(0));
    opt.sharp_cos_threshold = Env::read<geom_t>("SMESH_SHARP_EDGES_THRESHOLD", static_cast<geom_t>(0.15));
    opt.max_passes          = Env::read<int>("SMESH_IMPROVE_PASSES", 8);
    opt.smooth_iters        = Env::read<int>("SMESH_SMOOTH_ITERS", 8);
    opt.use_parametrization = true;

    if (improve(*mesh, opt) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[2]));
}
