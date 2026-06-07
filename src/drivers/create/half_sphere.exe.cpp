#include <stdio.h>
#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("half_sphere.exe");

    auto ctx = smesh::initialize_serial(argc, argv);

    if (argc != 6) {
        fprintf(stderr,
                "usage: %s <radius> <nx> <ny> <nz> "
                "<output_folder>\n",
                argv[0]);
        return SMESH_FAILURE;
    }

    {
        auto mesh = Mesh::create_hex8_half_sphere(
                ctx->communicator(), std::atof(argv[1]), std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
        mesh->write(Path(argv[5]));
    }

    return SMESH_SUCCESS;
}
