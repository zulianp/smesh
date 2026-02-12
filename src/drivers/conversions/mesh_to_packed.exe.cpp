#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("mesh_to_packed.exe");

  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 4) {
    fprintf(stderr, "Usage: %s <block_size> <input_folder> <output_folder>\n",
            argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[2]));
    auto packed = PackedMesh<i16>::create(mesh, {}, false, std::atoi(argv[1]));
    packed->write(Path(argv[3]));
  }

  return SMESH_SUCCESS;
}