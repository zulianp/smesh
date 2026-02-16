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

      array_write_convert_from_extension(output_folder / Path("n2e_ptr.int32"),
                                         n2eptr, n_local_nodes);

      array_write_convert_from_extension(output_folder / Path("n2e_idx.int32"),
                                         n2e_idx, n2eptr[n_local_nodes]);

      printf("n2e: #elements: %ld  #nodes: %ld\n", n_local_elements,
             n_local_nodes);

    }
#ifdef SMESH_ENABLE_MPI
    else {
      SMESH_TRACE_SCOPE("n2e+output (distributed)");
      create_n2e<idx_t, count_t, element_idx_t>(
          comm->get(), n_local_elements, n_global_elements, n_local_nodes,
          n_global_nodes, nnodesxelem, elems, &n2eptr, &n2e_idx);

      // Ensure it is always the same
      sort_n2e<count_t, element_idx_t>(n_local_nodes, n2eptr, n2e_idx);

      ptrdiff_t n2e_nidx = n2eptr[n_local_nodes];
      ptrdiff_t n2e_offset = 0;
      MPI_Exscan(&n2e_nidx, &n2e_offset, 1, smesh::mpi_type<ptrdiff_t>(),
                 MPI_SUM, comm->get());

      ptrdiff_t n_local_n2e_idx = n2eptr[n_local_nodes];
      ptrdiff_t n_global_n2e_idx = n2e_offset + n_local_n2e_idx;
      MPI_Bcast(&n_global_n2e_idx, 1, smesh::mpi_type<ptrdiff_t>(),
                comm->size() - 1, comm->get());

      for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
        n2eptr[i] = n2eptr[i] + n2e_offset;
      }

      array_write_convert_from_extension(comm->get(),
                                         output_folder / Path("n2e_ptr.int32"),
                                         n2eptr, n_local_nodes, n_global_nodes);

      array_write_convert_from_extension(
          comm->get(), output_folder / Path("n2e_idx.int32"), n2e_idx,
          n_local_n2e_idx, n_global_n2e_idx);
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
