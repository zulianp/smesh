#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_env.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("refine.exe");
  auto ctx  = initialize(argc, argv);
  auto comm = ctx->communicator();

  if (argc != 3) {
    if (!comm->rank()) {
      fprintf(stderr, "Usage: %s <mesh_folder> <output_folder>\n", argv[0]);
    }
    return SMESH_FAILURE;
  }

  auto mesh = Mesh::create_from_file(comm, Path(argv[1]));
  if (!mesh) {
    return SMESH_FAILURE;
  }
  auto refined = refine(mesh, Env::read<int>("SMESH_REFINEMENT_LEVELS", 1));
  if (!refined) {
    return SMESH_FAILURE;
  }
  return refined->write(Path(argv[2]));
}
