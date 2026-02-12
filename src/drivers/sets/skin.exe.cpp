#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_sideset.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("create_sideset.exe");
  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 3) {
    fprintf(stderr,
            "Usage: %s <mesh_folder> "
            "<output_folder>\n",
            argv[0]);
    return SMESH_FAILURE;
  }

  int ret = SMESH_SUCCESS;
  {
    auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));
    auto sideset = skin_sideset(mesh);
    auto surf = mesh_from_sideset(mesh, sideset);
    surf->write(Path(argv[2]));
    sideset->write(Path(argv[2]) / "parent_sideset");
  }

  return ret;
}