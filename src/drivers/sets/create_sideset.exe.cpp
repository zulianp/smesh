#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("create_sideset.exe");
  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 4) {
    fprintf(stderr,
            "Usage: %s <mesh_folder> <x> <y> <z> <angle_threshold> "
            "<output_folder>\n",
            argv[0]);
    return SMESH_FAILURE;
  }

  int ret = SMESH_SUCCESS;
  {
    auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));

  }

  return ret;
}