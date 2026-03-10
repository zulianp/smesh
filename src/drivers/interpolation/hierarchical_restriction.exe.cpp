#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_restriction.hpp"
#include "smesh_tracer.hpp"

using namespace smesh;
int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);
    SMESH_TRACE_SCOPE("hierarchical_restriction");

    if (argc != 6) {
        fprintf(stderr,
                "usage: %s <mesh> <from_element> <to_element> <input.float64> "
                "<output.float64>\n",
                argv[0]);
        return SMESH_FAILURE;
    }

    auto            folder       = Path(argv[1]);
    smesh::ElemType from_element = type_from_string(argv[2]);
    smesh::ElemType to_element   = type_from_string(argv[3]);
    auto            path_input   = Path(argv[4]);
    auto            path_output  = Path(argv[5]);

    bool SFEM_USE_CRS_GRAPH_RESTRICT = Env::read("SFEM_USE_CRS_GRAPH_RESTRICT", false);

    auto            mesh       = Mesh::create_from_file(ctx->communicator(), folder);
    const ptrdiff_t n_elements = mesh->n_elements();
    const ptrdiff_t n_nodes    = mesh->n_nodes();

    ptrdiff_t n_coarse_nodes = max_node_id(to_element, n_elements, mesh->elements(0)->data()) + 1;

    auto from   = Buffer<real_t>::from_file(path_input);
    auto to     = create_host_buffer<real_t>(n_coarse_nodes);
    auto b_from = from->data();
    auto b_to   = to->data();

    if (!SFEM_USE_CRS_GRAPH_RESTRICT) {
        int  nxe                               = elem_num_nodes(mesh->element_type(0));
        auto element_to_node_incidence_count   = create_host_buffer<uint16_t>(n_nodes);
        auto b_element_to_node_incidence_count = element_to_node_incidence_count->data();

        auto b_elements = mesh->elements(0)->data();
        for (int d = 0; d < nxe; d++) {
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < n_elements; ++i) {
#pragma omp atomic update
                b_element_to_node_incidence_count[b_elements[d][i]]++;
            }
        }

        hierarchical_restriction(from_element,
                                 to_element,
                                 n_elements,
                                 mesh->elements(0)->data(),
                                 element_to_node_incidence_count->data(),
                                 1,
                                 from->data(),
                                 to->data());
    } else {
        auto crs_graph = mesh->node_to_node_graph();
        auto rowptr    = crs_graph->rowptr()->data();
        auto colidx    = crs_graph->colidx()->data();
        hierarchical_restriction(n_coarse_nodes, rowptr, colidx, 1, from->data(), to->data());
    }

    to->to_file(path_output);

    return SMESH_SUCCESS;
}
