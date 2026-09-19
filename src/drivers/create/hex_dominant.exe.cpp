#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("hex_dominant.exe");
    auto ctx  = initialize_serial(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 2) {
        if (!comm->rank()) {
            fprintf(stderr,
                    "Usage: %s <output_folder>\n"
                    "  Serial HEX-dominant unit: HEX8 + PYRAMID5 + WEDGE6 + TET4.\n",
                    argv[0]);
        }
        return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_hex_dominant_serial(comm);
    if (!mesh) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[1]));
}
