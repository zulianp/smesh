#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "smesh_context.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_common.hpp"
#include "smesh_graph.hpp"
#include "smesh_tracer.hpp"
#include "smesh_env.hpp"

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);
    SMESH_TRACE_SCOPE("adjaciency_matrix.exe");

    if (argc < 2) {
        fprintf(stderr, "usage: %s <folder> [output_folder=./]\n", argv[0]);
        return EXIT_FAILURE;
    }

    auto output_folder = Path("./");
    if (argc > 2) {
        output_folder = Path(argv[2]);
    }

    bool SFEM_EXPORT_FP32 = Env::read("SFEM_EXPORT_FP32", false);
    bool SFEM_REMOVE_DIAGONAL = Env::read("SFEM_REMOVE_DIAGONAL", true);
    bool SFEM_GRAPH_LAPLACIAN = Env::read("SFEM_GRAPH_LAPLACIAN", false);
    bool SFEM_NORMALIZE_ROWS = Env::read("SFEM_NORMALIZE_ROWS", false);
    bool SFEM_REMOVE_LOWER_TRIANGULAR = Env::read("SFEM_REMOVE_LOWER_TRIANGULAR", false);
    bool SFEM_NEGATE_LOWER_TRIANGULAR = Env::read("SFEM_NEGATE_LOWER_TRIANGULAR", false);
    bool SFEM_DIRECTED = Env::read("SFEM_DIRECTED", false);

    std::string ext = SFEM_EXPORT_FP32? "float32" : "float64";


    auto folder = Path(argv[1]);
    auto mesh               = Mesh::create_from_file(ctx->communicator(), folder);
    const ptrdiff_t n_nodes = mesh->n_nodes();

    auto crs_graph = mesh->node_to_node_graph();
    const ptrdiff_t nnz = crs_graph->nnz();
    auto rowptr = crs_graph->rowptr()->data();
    auto colidx = crs_graph->colidx()->data();
    auto values = create_host_buffer<real_t>(crs_graph->nnz());
    auto b_values = values->data();

    if (SFEM_DIRECTED) {
        printf("SFEM_DIRECTED\n");
        for (ptrdiff_t i = 0; i < n_nodes; i++) {
            for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                idx_t col = colidx[k];
                if (col != i) {
                    if (col < i) {
                        b_values[k] = 1;
                    }
                    // else {
                    //     values[k] = -1;
                    // }
                }
            }
        }
    } else {
        if (SFEM_GRAPH_LAPLACIAN) {
            for (count_t i = 0; i < nnz; i++) {
                b_values[i] = 1;
            }

            for (ptrdiff_t i = 0; i < n_nodes; i++) {
                real_t row_sum = 0;
                for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                    row_sum += b_values[k];
                }

                for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                    idx_t col = colidx[k];
                    if (col == i) {
                        b_values[k] = 1 - row_sum;
                    }
                }
            }

        } else {
            for (count_t i = 0; i < nnz; i++) {
                b_values[i] = 1;
            }

            if (SFEM_REMOVE_DIAGONAL) {
                for (ptrdiff_t i = 0; i < n_nodes; i++) {
                    for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                        idx_t col = colidx[k];

                        if (col == i) {
                            b_values[k] = 0;
                        }
                    }
                }
            }

            if (SFEM_REMOVE_LOWER_TRIANGULAR) {
                for (ptrdiff_t i = 0; i < n_nodes; i++) {
                    for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                        idx_t col = colidx[k];
                        if (col < i) {
                            b_values[k] = 0;
                        }
                    }
                }
            }

            if (SFEM_NEGATE_LOWER_TRIANGULAR) {
                for (ptrdiff_t i = 0; i < n_nodes; i++) {
                    for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                        idx_t col = colidx[k];
                        if (col < i) {
                            b_values[k] = -b_values[k];
                        }
                    }
                }
            }

            if (SFEM_NORMALIZE_ROWS) {
                for (ptrdiff_t i = 0; i < n_nodes; i++) {
                    real_t row_sum = 0;
                    for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                        row_sum += b_values[k];
                    }

                    for (count_t k = rowptr[i]; k < rowptr[i + 1]; k++) {
                        b_values[k] /= row_sum;
                    }
                }
            }
        }
    }

    

    crs_graph->rowptr()->to_file(output_folder / ("rowptr." + str(TypeToString<count_t>::value())));
    crs_graph->colidx()->to_file(output_folder / ("colidx." + str(TypeToString<idx_t>::value())));
    values->to_file(output_folder / ("values." + ext));

    return SMESH_SUCCESS;
}
