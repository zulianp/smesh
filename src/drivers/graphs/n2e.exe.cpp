#include "smesh_context.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("n2e.exe");
  auto ctx = initialize(argc, argv);

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <mesh_folder> <output_folder>\n", argv[0]);
    return SMESH_FAILURE;
  }

  auto mesh_folder = Path(argv[1]);
  auto output_folder = Path(argv[2]);

  auto comm = ctx->communicator();

  int nnodesxelem;
  idx_t **elems;
  ptrdiff_t n_local_elements;
  ptrdiff_t n_global_elements;
  mesh_block_from_folder(comm->get(), mesh_folder, &nnodesxelem, &elems,
                         &n_local_elements, &n_global_elements);

  int spatial_dim;
  geom_t **points;
  ptrdiff_t n_local_nodes;
  ptrdiff_t n_global_nodes;
  mesh_coordinates_from_folder(comm->get(), mesh_folder, &spatial_dim, &points,
                               &n_local_nodes, &n_global_nodes);

  if (comm->size() == 1) {
    count_t *n2eptr;
    element_idx_t *n2e_idx;
    create_n2e(n_local_elements, n_local_nodes, nnodesxelem, elems, &n2eptr,
               &n2e_idx);

    // for (ptrdiff_t i = 0; i < n_local_nodes; i++) {
    //   count_t begin = n2eptr[i];
    //   count_t end = n2eptr[i + 1];
    //   for (count_t j = begin; j < end; j++) {
    //     std::cout << n2e_idx[j] << " ";
    //   }
    //   std::cout << "\n";
    // }
    // std::cout << "\n";
    // std::cout << "------------\n";
  }

//   comm->print_callback([&](std::ostream &os) {
//     os << "N nodesxelem: " << nnodesxelem << "\n";
//     os << "N elements: " << n_local_elements << "\n";
//     os << "N global elements: " << n_global_elements << "\n";
//   });

  count_t *n2eptr;
  element_idx_t *n2e_idx;
  create_n2e<idx_t, count_t, element_idx_t>(
      comm->get(), n_local_elements, n_global_elements, n_local_nodes, n_global_nodes, nnodesxelem,
      elems, &n2eptr, &n2e_idx);

      if (comm->rank() == 0) {
        printf("#elements: %ld  #nodes: %ld\n", n_global_elements, n_global_nodes);
      }

//   comm->print_callback([&](std::ostream &os) {
//     for (ptrdiff_t i = 0; i < n_local_nodes; i++) {
//       count_t begin = n2eptr[i];
//       count_t end = n2eptr[i + 1];
//       for (count_t j = begin; j < end; j++) {
//         os << n2e_idx[j] << " ";
//       }
//       os << "\n";
//     }
//     os << "\n";
//   });

  free(n2eptr);
  free(n2e_idx);
  free(elems);
  free(points);

  return SMESH_SUCCESS;
}
