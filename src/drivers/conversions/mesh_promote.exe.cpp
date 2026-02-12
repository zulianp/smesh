#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("mesh_promote.exe");

  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 4) {
    fprintf(stderr, "Usage: %s <to_element_type> <input_mesh> <output_mesh>\n", argv[0]);
    return SMESH_FAILURE;
  }

  int ret = SMESH_SUCCESS;
  {
    auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));
    auto new_mesh = promote_to(type_from_string(argv[2]), mesh);
    ret = new_mesh->write(Path(argv[3]));
  }

  return ret;
}