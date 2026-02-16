#include "smesh_context.hpp"
#include "smesh_glob.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_distributed_write.hpp"
#endif

using namespace smesh;

int main(int argc, char **argv) {

  SMESH_TRACE_SCOPE("n2e.exe");
  auto ctx = initialize(argc, argv);

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <mesh_folder> <output_folder>\n", argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto mesh_folder = Path(argv[1]);
    auto output_folder = Path(argv[2]);

    auto comm = ctx->communicator();

    int nnodesxelem;
    idx_t **elems;
    ptrdiff_t n_local_elements;

#ifdef SMESH_ENABLE_MPI
    ptrdiff_t n_global_elements;
    mesh_block_from_folder<idx_t>(comm->get(), mesh_folder, &nnodesxelem,
                                  &elems, &n_local_elements,
                                  &n_global_elements);
#else
    mesh_block_from_folder<idx_t>(mesh_folder, &nnodesxelem, &elems,
                                  &n_local_elements);
#endif

    int spatial_dim;
    geom_t **points;
    ptrdiff_t n_local_nodes;

#ifdef SMESH_ENABLE_MPI
    ptrdiff_t n_global_nodes;
    mesh_coordinates_from_folder(comm->get(), mesh_folder, &spatial_dim,
                                 &points, &n_local_nodes, &n_global_nodes);
#else
    mesh_coordinates_from_folder<geom_t>(mesh_folder, &spatial_dim, &points,
                                         &n_local_nodes);
#endif

    printf("Memory (elements): %g [GB]\n",
           (n_local_nodes * nnodesxelem) * sizeof(element_idx_t) * 1e-9);

    count_t *n2eptr;
    element_idx_t *n2e_idx;

    if (!comm->rank()) {
      create_directory(output_folder);
    }

#ifdef SMESH_ENABLE_MPI
    comm->barrier();

    if (comm->size() == 1)
#endif
    {
      SMESH_TRACE_SCOPE("n2e+output (serial)");

      create_n2e<idx_t, count_t, element_idx_t>(n_local_elements, n_local_nodes,
                                                nnodesxelem, elems, &n2eptr,
                                                &n2e_idx);
      // Ensure it is always the same
      sort_n2e<count_t, element_idx_t>(n_local_nodes, n2eptr, n2e_idx);

      printf("Memory (n2e): %g [GB]\n",
             (n_local_nodes + 1) * sizeof(count_t) * 1e-9 +
                 n2eptr[n_local_nodes] * sizeof(element_idx_t) * 1e-9);

      count_t *n2n_ptr;
      idx_t *n2n_idx;
      create_n2n_from_n2e(n_local_elements, n_local_nodes, nnodesxelem, elems,
                          n2eptr, n2e_idx, &n2n_ptr, &n2n_idx);

      printf("Memory (n2n): %g [GB]\n",
             (n_local_nodes + 1) * sizeof(count_t) * 1e-9 +
                 n2n_ptr[n_local_nodes] * sizeof(idx_t) * 1e-9);

      array_write_convert_from_extension(output_folder / Path("n2n_ptr.int32"),
                                         n2n_ptr, n_local_nodes);

      array_write_convert_from_extension(output_folder / Path("n2n_idx.int32"),
                                         n2n_idx, n2n_ptr[n_local_nodes]);
    }
#ifdef SMESH_ENABLE_MPI
    else {
      SMESH_TRACE_SCOPE("n2e+output (distributed)");

      create_n2e<idx_t, count_t, element_idx_t>(
          comm->get(), n_local_elements, n_global_elements, n_local_nodes,
          n_global_nodes, nnodesxelem, elems, &n2eptr, &n2e_idx);

      // Ensure it is always the same
      sort_n2e<count_t, element_idx_t>(n_local_nodes, n2eptr, n2e_idx);


      SMESH_ERROR("Not implemented");
    }

    if (!comm->rank()) {
      printf("#elements: %ld  #nodes: %ld\n", n_global_elements,
             n_global_nodes);
    }

#endif // SMESH_ENABLE_MPI
    free(n2eptr);
    free(n2e_idx);
    free(elems);
    free(points);
  }

  return SMESH_SUCCESS;
}
