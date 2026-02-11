#include <stdio.h>
#include "smesh_path.hpp"
#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"

using namespace smesh;

int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <block_size> <input_folder> <output_folder>\n", argv[0]);
        return 1;
    }

    auto ctx = smesh::initialize_serial(argc, argv);
    {
        auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));
        auto packed = PackedMesh<i16>::create(mesh);
        packed->write(Path(argv[2]));
    }
    // int block_size = std::atoi(argv[1]);
    // Path input_folder(argv[2]);
    // Path output_folder(argv[3]);

    return SMESH_SUCCESS;
}