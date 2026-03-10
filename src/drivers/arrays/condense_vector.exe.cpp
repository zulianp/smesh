// #include "smesh_base.hpp"
// #include "smesh_context.hpp"
// #include "smesh_types.hpp"
// #include "smesh_buffer.hpp"
// #include "smesh_write.hpp"
// #include "smesh_env.hpp"


// using namespace smesh;

// int main(int argc, char *argv[]) {
//     auto ctx = initialize_serial(argc, argv);

//     auto input_path = Path(argv[1]);
//     auto dirichlet_nodes_path = Path(argv[2]);
//     auto output_path = Path(argv[3]);

//     if (argc != 3) {
//         fprintf(stderr, "usage: %s <array.raw> <dirichlet_nodes.raw> [out=condensed.raw]\n", argv[0]);
//         return EXIT_FAILURE;
//     }

//     if (input_path.to_string() == output_path.to_string()) {
//         fprintf(stderr, "Input and output are the same! Quitting!\n");
//         fprintf(stderr, help, argv[0]);
//         return EXIT_FAILURE;
//     }

//     auto values = Buffer<real_t>::from_file(input_path);
//     auto dirichlet_nodes = Buffer<idx_t>::from_file(dirichlet_nodes_path);
//     const ptrdiff_t nnodes = values->size();
//     const ptrdiff_t ndirichlet = dirichlet_nodes->size();
//     const ptrdiff_t new_nnodes = nnodes - ndirichlet; 

//     auto is_dirichlet = create_host_buffer<idx_t>(nnodes);
//     auto new_values = create_host_buffer<real_t>(new_nnodes);

//     auto b_is_dirichlet = is_dirichlet->data();
//     auto b_values = values->data();
//     auto b_new_values = new_values->data();

//     for (ptrdiff_t node = 0; node < nnodes; ++node) {
//         if (!b_is_dirichlet[node]) {
//             b_new_values[node] = b_values[node];
//         }
//     }

//     new_values->to_file(output_path);
//     return EXIT_SUCCESS;
// }


int main() { return 1; }