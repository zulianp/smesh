// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>

// #include "smesh_base.hpp"
#include "smesh_alloc.hpp"
// #include "smesh_context.hpp"
// #include "smesh_types.hpp"
// #include "smesh_buffer.hpp"
// #include "smesh_write.hpp"
// #include "smesh_env.hpp"

// using namespace smesh;

// int main(int argc, char *argv[]) {
//     auto ctx = initialize_serial(argc, argv);

//     auto condensed_path = Path(argv[1]);
//     auto dirichlet_nodes_path = Path(argv[2]);
//     auto output_path = Path(argv[3]);

//     if (argc != 3) {
//         fprintf(stderr, "usage: %s <condensed.raw> <dirichlet_nodes.raw> [out=full.raw]\n", argv[0]);
//         fprintf(stderr, help, argv[0]);
//         return EXIT_FAILURE;
//     }

//     const char *output_path = "./full.raw";

//     if (argc > 3) {
//         output_path = argv[3];
//     }

//     if (strcmp(output_path, argv[1]) == 0) {
//         fprintf(stderr, "Input and output are the same! Quitting!\n");
//         fprintf(stderr, help, argv[0]);
//         return EXIT_FAILURE;
//     }

//     double tick = MPI_Wtime();

//     MPI_Datatype values_mpi_t = MPI_DOUBLE;

//     real_t *values;
//     ptrdiff_t nlocal_, nnodes;
//     array_create_from_file(comm, argv[1], values_mpi_t, (void **)&values, &nlocal_, &nnodes);

//     idx_t *is_dirichlet = 0;
//     ptrdiff_t new_nnodes = 0;
//     {
//         idx_t *dirichlet_nodes = 0;
//         ptrdiff_t nlocal_, ndirichlet;
//         array_create_from_file(
//             comm, argv[2], MPI_INT, (void **)&dirichlet_nodes, &nlocal_, &ndirichlet);

//         new_nnodes = nnodes + ndirichlet;

//         is_dirichlet = (idx_t *)SMESH_ALLOC(new_nnodes * sizeof(idx_t));
//         memset(is_dirichlet, 0, new_nnodes * sizeof(idx_t));

//         for (ptrdiff_t node = 0; node < ndirichlet; ++node) {
//             idx_t i = dirichlet_nodes[node];
//             is_dirichlet[i] = 1;
//         }

//         SMESH_FREE(dirichlet_nodes);
//     }

//     real_t *new_values = (real_t *)SMESH_ALLOC(new_nnodes * sizeof(real_t));
//     memset(new_values, 0, new_nnodes * sizeof(real_t));

//     for (ptrdiff_t node = 0, old_node_idx = 0; node < new_nnodes; ++node) {
//         if (!is_dirichlet[node]) {
//             new_values[node] = values[old_node_idx++];
//         }
//     }

//     array_write(comm, output_path, values_mpi_t, (void *)new_values, new_nnodes, new_nnodes);

//     SMESH_FREE(is_dirichlet);
//     SMESH_FREE(values);
//     SMESH_FREE(new_values);

//     double tock = MPI_Wtime();

//     if (!rank) {
//         printf("Remapped dofs: from %ld to %ld\n", (long)nnodes, (long)new_nnodes);
//         printf("TTS: %g seconds\n", tock - tick);
//     }

//     MPI_Finalize();
//     return EXIT_SUCCESS;
// }


int main() { return 1; }