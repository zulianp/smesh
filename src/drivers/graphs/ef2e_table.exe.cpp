#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_common.hpp"

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc < 2) {
        fprintf(stderr, "usage: %s <folder> <output_folder>\n", argv[0]);
        return EXIT_FAILURE;
    }

    auto output_folder = Path(argv[2]);
    create_directory(output_folder);

    double tick = time_seconds();

    ///////////////////////////////////////////////////////////////////////////////
    // Read data
    ///////////////////////////////////////////////////////////////////////////////

    auto            folder       = Path(argv[1]);
    auto            mesh         = Mesh::create_from_file(ctx->communicator(), folder);
    smesh::ElemType element_type = mesh->element_type(0);

    if (is_semistructured_type(element_type)) {
        mesh = derefine(mesh, 1);
    }

    auto table = mesh->half_face_table();
    table->to_file(output_folder / ("ef2e_table." + str(TypeToString<element_idx_t>::value())));

    double tock = time_seconds();

    printf("----------------------------------------\n");
    printf("TTS:\t\t\t%g seconds\n", tock - tick);

    return SMESH_SUCCESS;
}
