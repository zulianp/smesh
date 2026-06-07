#include <stdio.h>
#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("half_sphere.exe");

    auto ctx = smesh::initialize_serial(argc, argv);

    if (argc != 7) {
        fprintf(stderr,
                "usage: %s <element_type> <radius> <nx> <ny> <nz> "
                "<output_folder>\n",
                argv[0]);
        return SMESH_FAILURE;
    }

    {
        const enum ElemType element_type = type_from_string(argv[1]);
        auto mesh = Mesh::create_half_sphere(
                ctx->communicator(), element_type, std::atof(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]));
        mesh->write(Path(argv[6]));
    }

    return SMESH_SUCCESS;
}
