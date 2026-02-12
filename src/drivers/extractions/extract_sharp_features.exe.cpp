#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_extractions.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("extract_sharp_features.exe");
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

    auto sharp_edges = extract_sharp_edges(
        *mesh, Env::read("SMESH_SHARP_EDGES_THRESHOLD", 0.15));

    auto sharp_corners =
        extract_sharp_corners(mesh->n_nodes(), sharp_edges, true);

    auto disconnected_faces = extract_disconnected_faces(*mesh, sharp_edges);

    auto output_path = Path(argv[2]);

    create_directory(output_path);
    create_directory(output_path / ("edges"));
    create_directory(output_path / ("corners"));
    create_directory(output_path / ("disconnected_faces"));

    {
      std::string ext(TypeToString<idx_t>::value());
      sharp_edges->to_files(output_path / ("edges/i%d." + ext));
      sharp_corners->to_file(output_path / ("corners/i0." + ext));
    }

    {
      std::string ext(TypeToString<element_idx_t>::value());
      disconnected_faces->to_file(output_path /
                                  ("disconnected_faces/i0." + ext));
    }
  }

  return ret;
}