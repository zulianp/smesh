#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_extractions.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_file_extensions.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("extract_sharp_features.exe");
    auto ctx = initialize_serial(argc, argv);

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <mesh_folder> <output_folder>\n", argv[0]);
        return SMESH_FAILURE;
    }

    int ret = SMESH_SUCCESS;
    {
        auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));

        auto sharp_edges = extract_sharp_edges(*mesh, Env::read("SMESH_SHARP_EDGES_THRESHOLD", 0.15));
        if (!sharp_edges) {
            return SMESH_FAILURE;
        }

        auto sharp_corners = extract_sharp_corners(*mesh, sharp_edges, true);
        if (!sharp_corners) {
            return SMESH_FAILURE;
        }

        auto disconnected_faces = extract_disconnected_faces(*mesh, *sharp_edges);

        mesh->add_edgeset("sharp_edges", sharp_edges);
        mesh->add_nodeset("sharp_corners", sharp_corners);

        auto output_path = Path(argv[2]);
        if (mesh->write(output_path) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }

        if (disconnected_faces) {
            create_directory(output_path / "disconnected_faces");
            std::string ext(TypeToString<element_idx_t>::value());
            disconnected_faces->to_file(output_path / ("disconnected_faces/i0." + ext));
        }
    }

    return ret;
}
