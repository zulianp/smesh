#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "smesh_context.hpp"
#include "smesh_types.hpp"
#include "smesh_buffer.hpp"
#include "smesh_write.hpp"
#include "smesh_env.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_graph.hpp"
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

    auto folder = Path(argv[1]);
    auto        mesh = Mesh::create_from_file(ctx->communicator(), folder);
    ElemType element_type          = mesh->element_type(0);
    ElemType element_type_for_algo = element_type;
    if (element_type == TET10) {
        element_type_for_algo = TET4;
    } else if (element_type == TRI6) {
        element_type_for_algo = TRI3;
    }

    // Read only the data we need
    const int nnxe = elem_num_nodes(element_type_for_algo);
    idx_t **const elems = mesh->elements(0)->data();
    const ptrdiff_t n_elements = mesh->n_elements();
    const ptrdiff_t n_nodes = mesh->n_nodes();

    count_t       *adj_ptr = 0;
    element_idx_t *adj_idx = 0;
    create_dual_graph(n_elements, n_nodes, element_type_for_algo, elems, &adj_ptr, &adj_idx);
    auto adj_ptr_buffer = manage_host_buffer<count_t>(n_elements + 1, adj_ptr);
    auto adj_idx_buffer = manage_host_buffer<element_idx_t>(adj_ptr[n_elements], adj_idx);
    adj_ptr_buffer->to_file(output_folder / ("adj_ptr." + str(TypeToString<count_t>::value())));
    adj_idx_buffer->to_file(output_folder / ("adj_idx." + str(TypeToString<element_idx_t>::value())));

    double tock = time_seconds();
    printf("Dual graph %ld elements, %ld nnz\n", (long)n_elements, (long)adj_ptr[n_elements]);
    printf("TTS:\t\t\t%g seconds\n", tock - tick);
    return SMESH_SUCCESS;
}
