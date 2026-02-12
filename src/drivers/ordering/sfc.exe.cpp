#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh_reorder.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("sfc.exe");
  auto ctx = initialize_serial(argc, argv);

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
    auto sfc = SFC::create_from_env();
    ret = sfc->reorder(*mesh);

    if (ret != SMESH_SUCCESS) {
      mesh->write(Path(argv[2]));
    }
  }

  return ret;
}