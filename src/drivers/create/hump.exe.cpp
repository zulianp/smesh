#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cstdlib>
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("hump.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 12) {
        if (!comm->rank()) {
            fprintf(stderr,
                    "Usage: %s <element_type> <nx> <ny> <nz> <length> <height> <width> "
                    "<hump_start> <hump_length> <hump_height> <output_folder>\n"
                    "  Wall-mounted hump. HEX8 or TET4.\n",
                    argv[0]);
        }
        return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_wall_mounted_hump(comm,
                                               type_from_string(argv[1]),
                                               std::atoi(argv[2]),
                                               std::atoi(argv[3]),
                                               std::atoi(argv[4]),
                                               std::atof(argv[5]),
                                               std::atof(argv[6]),
                                               std::atof(argv[7]),
                                               std::atof(argv[8]),
                                               std::atof(argv[9]),
                                               std::atof(argv[10]));
    if (!mesh) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[11]));
}
