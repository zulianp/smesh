#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("bidomain_cube.exe");

  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 11) {
    fprintf(stderr,
            "Usage: %s <nx> <ny> <nz> <xmin> <ymin> <zmin> "
            "<xmax> <ymax> <zmax> <output_folder>n",
            argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto mesh = Mesh::create_hex8_bidomain_cube(
        ctx->communicator(), std::atoi(argv[1]),
        std::atoi(argv[2]), std::atoi(argv[3]), std::atof(argv[4]),
        std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]),
        std::atof(argv[8]), std::atof(argv[9]));
    mesh->write(Path(argv[10]));
  }

  return SMESH_SUCCESS;
}