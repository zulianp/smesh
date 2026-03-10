#include "smesh_base.hpp"

#include "smesh_buffer.hpp"
#include "smesh_common.hpp"
#include "smesh_context.hpp"
#include "smesh_glob.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"
#include "smesh_env.hpp"  
#include "smesh_mesh.hpp"

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc != 3) {
        fprintf(stderr, "usage: %s <folder> <output>\n", argv[0]);
        return EXIT_FAILURE;
    }

    auto path_count = Path(argv[2]);
    auto folder = Path(argv[1]);
    auto m = Mesh::create_from_file(ctx->communicator(), folder);

    auto count = create_host_buffer<int>(m->n_nodes());
    int nxe = elem_num_nodes(m->element_type(0));

    const ptrdiff_t nelements = m->n_elements();
    const auto elements = m->elements(0)->data();
    auto b_count = count->data();
    for (int d = 0; d < nxe; d++) {
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < nelements; ++i) {
#pragma omp atomic update
            b_count[elements[d][i]]++;
        }
    }

    count->to_file(path_count);

    return SMESH_SUCCESS;
}
