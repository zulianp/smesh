#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("square.exe");

  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 9) {
    fprintf(stderr,
            "Usage: %s <element_type> <nx> <ny> <xmin> <ymin> "
            "<xmax> <ymax> <output_folder>\n",
            argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto mesh = Mesh::create_square(
        ctx->communicator(), type_from_string(argv[1]), std::atoi(argv[2]),
        std::atoi(argv[3]), std::atof(argv[4]), std::atof(argv[5]),
        std::atof(argv[6]), std::atof(argv[7]));
    mesh->write(Path(argv[8]));
  }

  return SMESH_SUCCESS;
}