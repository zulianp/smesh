#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_prolongation.hpp"
#include "smesh_tracer.hpp"

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);
    SMESH_TRACE_SCOPE("hierarchical_prolongation");

    if (argc != 6) {
        fprintf(stderr,
                "usage: %s <mesh> <from_element> <to_element> <input.float64> "
                "<output.float64>\n",
                argv[0]);
        return SMESH_FAILURE;
    }

    auto folder       = Path(argv[1]);
    auto from_element = type_from_string(argv[2]);
    auto to_element   = type_from_string(argv[3]);
    auto path_input   = Path(argv[4]);
    auto path_output  = Path(argv[5]);

    auto            mesh       = Mesh::create_from_file(ctx->communicator(), folder);
    const ptrdiff_t n_elements = mesh->n_elements();
    const ptrdiff_t n_nodes    = mesh->n_nodes();

    auto from   = Buffer<real_t>::from_file(path_input);
    auto to     = create_host_buffer<real_t>(n_nodes);
    auto b_from = from->data();
    auto b_to   = to->data();

    hierarchical_prolongation(from_element, to_element, n_elements, mesh->elements(0)->data(), 1, from->data(), to->data()) ||
            to->to_file(path_output);

    return SMESH_SUCCESS;
}
