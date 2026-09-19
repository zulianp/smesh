#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cstdlib>
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("lshape.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 10) {
        if (!comm->rank()) {
            fprintf(stderr,
                    "Usage: %s <nx> <ny> <nz> <xmax> <ymax> <zmax> <step_x> <step_y> "
                    "<output_folder>\n"
                    "  HEX8 backward-facing step. Notch [0,step_x] x [0,step_y] x [0,zmax] "
                    "removed. The step must fall on a grid line.\n",
                    argv[0]);
        }
        return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_hex8_lshape(comm,
                                         std::atoi(argv[1]),
                                         std::atoi(argv[2]),
                                         std::atoi(argv[3]),
                                         std::atof(argv[4]),
                                         std::atof(argv[5]),
                                         std::atof(argv[6]),
                                         std::atof(argv[7]),
                                         std::atof(argv[8]));
    if (!mesh) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[9]));
}
