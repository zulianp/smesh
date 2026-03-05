#include "smesh_context.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_env.hpp"
#include "smesh_gencube.hpp"


#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  using namespace smesh;

  SMESH_TRACE_SCOPE("cube.exe");

  auto ctx = initialize_serial(argc, argv);

  if (argc != 12) {
    fprintf(stderr,
            "Usage: %s <element_type> <nx> <ny> <nz> <xmin> <ymin> <zmin> "
            "<xmax> <ymax> <zmax> <output_folder>\n",
            argv[0]);
    return SMESH_FAILURE;
  }

  const enum ElemType element_type = type_from_string(argv[1]);
  const ptrdiff_t nx = std::atoi(argv[2]);
  const ptrdiff_t ny = std::atoi(argv[3]);
  const ptrdiff_t nz = std::atoi(argv[4]);
  const f32 xmin = std::atof(argv[5]);
  const f32 ymin = std::atof(argv[6]);
  const f32 zmin = std::atof(argv[7]);
  const f32 xmax = std::atof(argv[8]);
  const f32 ymax = std::atof(argv[9]);
  const f32 zmax = std::atof(argv[10]);
  const Path output_folder = Path(argv[11]);

  i64 chunk_size = Env::read<i64>("SMESH_CHUNK_SIZE", (i64)(2000));
  const i64 n_elements = nx * ny * nz;
  if(n_elements > chunk_size && element_type == HEX8) {
    mesh_hex8_cube_to_folder<i64, f32>(
      output_folder, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax, chunk_size);

  } else {
    auto mesh = Mesh::create_cube(
        ctx->communicator(), element_type, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax);
    mesh->write(output_folder);
  }

  return SMESH_SUCCESS;
}