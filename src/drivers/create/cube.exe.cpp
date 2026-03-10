#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_gencube.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh_reorder.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_sideset.hpp"
#include "smesh_tracer.hpp"

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

  i64 z_chunk_size = Env::read<i64>("SMESH_Z_CHUNK_SIZE", (i64)(1000));
  const i64 n_elements = nx * ny * nz;
  if ((n_elements > z_chunk_size * nx * ny ||
       n_elements > std::numeric_limits<int>::max()) &&
      element_type == HEX8) {
    // Huge meshes for testing!
    mesh_hex8_cube_to_folder<i64, f32>(output_folder, nx, ny, nz, xmin, ymin,
                                       zmin, xmax, ymax, zmax, z_chunk_size);
  } else {
    auto mesh = Mesh::create_cube(ctx->communicator(), element_type, nx, ny, nz,
                                  xmin, ymin, zmin, xmax, ymax, zmax);

    if (Env::read<bool>("SMESH_REORDER", true)) {
      auto reordering = SFC::create_from_env();
      reordering->reorder(*mesh);
    }

    mesh->write(output_folder);

    if (Env::read<bool>("SMESH_CREATE_SIDESETS", true)) {
      create_directory(output_folder / "sidesets");

      for (const auto &sideset :
           Sideset::create_from_plane(mesh, 0, 0, 1, zmin)) {
        sideset->write(output_folder / "sidesets/back");
      }

      for (const auto &sideset :
           Sideset::create_from_plane(mesh, 0, 0, 1, zmax)) {
        sideset->write(output_folder / "sidesets/front");
      }

      for (const auto &sideset :
           Sideset::create_from_plane(mesh, 0, 1, 0, ymin)) {
        sideset->write(output_folder / "sidesets/bottom");
      }

      for (const auto &sideset :
           Sideset::create_from_plane(mesh, 0, 1, 0, ymax)) {
        sideset->write(output_folder / "sidesets/top");
      }

      for (const auto &sideset :
           Sideset::create_from_plane(mesh, 1, 0, 0, xmin)) {
        sideset->write(output_folder / "sidesets/left");
      }

      for (const auto &sideset :
           Sideset::create_from_plane(mesh, 1, 0, 0, xmax)) {
        sideset->write(output_folder / "sidesets/right");
      }
    }
  }

  return SMESH_SUCCESS;
}