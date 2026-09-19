#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cstdlib>
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("ssquad_square.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 9) {
        if (!comm->rank()) {
            fprintf(stderr,
                    "Usage: %s <micro_per_dim> <nx> <ny> <xmin> <ymin> "
                    "<xmax> <ymax> <output_folder>\n"
                    "  PROTEUS_QUAD square (semistructured lattice).\n",
                    argv[0]);
        }
        return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_semistructured_quad_square(comm,
                                                        std::atoi(argv[1]),
                                                        std::atoi(argv[2]),
                                                        std::atoi(argv[3]),
                                                        std::atof(argv[4]),
                                                        std::atof(argv[5]),
                                                        std::atof(argv[6]),
                                                        std::atof(argv[7]));
    if (!mesh) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[8]));
}
