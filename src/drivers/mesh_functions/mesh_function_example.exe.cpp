#include "smesh_base.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

using namespace smesh;

void eval_function(const ptrdiff_t nnodes, geom_t **points, real_t *f) {
    SMESH_TRACE_SCOPE("eval_function");

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < nnodes; i++) {
        geom_t x = points[0][i];
        geom_t y = points[1][i];
        geom_t z = points[2][i];
        f[i]     = x * y * z;
    }
}

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc != 3) {
        fprintf(stderr, "usage: %s <mesh> <function.raw>\n", argv[0]);
        return EXIT_FAILURE;
    }

    auto        mesh   = Mesh::create_from_file(ctx->communicator(), smesh::Path(argv[1]));
    const char *output = argv[2];

    const ptrdiff_t n_nodes    = mesh->n_nodes();

    auto points = mesh->points()->data();
    auto f      = create_host_buffer<real_t>(n_nodes);
    eval_function(n_nodes, points, f->data());
    return f->to_file(smesh::Path(output));
}
