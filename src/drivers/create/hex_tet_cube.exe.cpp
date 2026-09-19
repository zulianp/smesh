#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cstdlib>
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("hex_tet_cube.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 11) {
        if (!comm->rank()) {
            fprintf(stderr,
                    "Usage: %s <nx> <ny> <nz> <xmin> <ymin> <zmin> "
                    "<xmax> <ymax> <zmax> <output_folder>\n"
                    "  HEX8 cube; the second half of hexes is converted to 6 TET4 each.\n",
                    argv[0]);
        }
        return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_hex8_tet4_cube(comm,
                                            std::atoi(argv[1]),
                                            std::atoi(argv[2]),
                                            std::atoi(argv[3]),
                                            std::atof(argv[4]),
                                            std::atof(argv[5]),
                                            std::atof(argv[6]),
                                            std::atof(argv[7]),
                                            std::atof(argv[8]),
                                            std::atof(argv[9]));
    if (!mesh) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[10]));
}
