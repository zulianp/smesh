#include "smesh_config.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_context.hpp"
#include "smesh_glob.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_sort.hpp"
#include "smesh_tracer.hpp"

#include "smesh_communicator.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_distributed_write.hpp"

using namespace smesh;

int main(int argc, char **argv) {

  SMESH_TRACE_SCOPE("pe2n.exe");
  auto ctx = initialize(argc, argv);

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <mesh_folder> <output_folder>\n", argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto mesh_folder = Path(argv[1]);
    auto output_folder = Path(argv[2]);

    auto comm = ctx->communicator();
    int comm_size = comm->size();
    int comm_rank = comm->rank();

    int nnodesxelem;
    idx_t **elems;
    ptrdiff_t n_local_elements;
    ptrdiff_t n_global_elements;
    mesh_block_from_folder<idx_t>(comm->get(), mesh_folder, &nnodesxelem,
                                  &elems, &n_local_elements,
                                  &n_global_elements);

    int spatial_dim;
    geom_t **points;
    ptrdiff_t n_local2global;
    ptrdiff_t n_global_nodes;
    mesh_coordinates_from_folder(comm->get(), mesh_folder, &spatial_dim,
                                 &points, &n_local2global, &n_global_nodes);

    count_t *n2eptr;
    element_idx_t *n2e_idx;

    if (!comm_rank) {
      create_directory(output_folder);
    }

    comm->barrier();

    SMESH_TRACE_SCOPE("pe2n+output (distributed)");
    create_n2e<idx_t, count_t, element_idx_t>(
        comm->get(), n_local_elements, n_global_elements, n_local2global,
        n_global_nodes, nnodesxelem, elems, &n2eptr, &n2e_idx);

    // Ensure it is always the same
    sort_n2e<count_t, element_idx_t>(n_local2global, n2eptr, n2e_idx);

    ptrdiff_t local2global_size = 0;
    idx_t *local2global = nullptr;
    count_t *local_n2e_ptr = nullptr;
    element_idx_t *local_n2e_idx = nullptr;
    redistribute_n2e(comm->get(), comm_size, comm_rank, n_local2global,
                     n_global_nodes, n_global_elements, n2eptr, n2e_idx,
                     &local2global_size, &local2global, &local_n2e_ptr,
                     &local_n2e_idx);

    // We do not need them anymore
    free(n2eptr);
    free(n2e_idx);

    idx_t **local_elements = (idx_t **)malloc(nnodesxelem * sizeof(idx_t *));
    for (int d = 0; d < nnodesxelem; ++d) {
      local_elements[d] = (idx_t *)malloc(n_local_elements * sizeof(idx_t));
    }

    localize_element_indices(comm_size, comm_rank, n_global_elements,
                             n_local_elements, nnodesxelem, elems,
                             local2global_size, local_n2e_ptr, local_n2e_idx,
                             local2global, local_elements);

    ptrdiff_t n_owned  = 0;
    ptrdiff_t n_shared = 0;
    ptrdiff_t n_ghosts = 0;
    rearrange_local_nodes(comm_size, comm_rank, n_global_elements,
                          n_local_elements, nnodesxelem, local2global_size,
                          local_n2e_ptr, local_n2e_idx, local2global,
                          local_elements, &n_owned, &n_shared, &n_ghosts);

    ptrdiff_t n_owned_not_shared = 0;
rearrange_local_elements(comm_size, comm_rank, n_global_elements,
                          n_local_elements, nnodesxelem, local2global_size,
                          local_n2e_ptr, local_n2e_idx, local_elements, n_owned,
                          &n_owned_not_shared);

    for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
      for (int d = 0; d < nnodesxelem; ++d) {
        const idx_t li = local_elements[d][i];
        if (li == invalid_idx<idx_t>() || li < 0 ||
            static_cast<ptrdiff_t>(li) >= local2global_size) {
          SMESH_ERROR(
              "pe2n.exe: unmapped node for local element %ld (d=%d) on rank "
              "%d (li=%d recv_nodes_size=%ld)\n",
              (long)i, d, comm_rank, (int)li, (long)local2global_size);
        }
      }
    }

    
    

    if (!comm_rank) {
      printf("#elements: %ld  #nodes: %ld\n", n_global_elements,
             n_global_nodes);
    }

    free(elems);
    free(points);

    free(local2global);
    free(local_n2e_ptr);
    free(local_n2e_idx);
    for (int d = 0; d < nnodesxelem; ++d) {
      free(local_elements[d]);
    }
    free(local_elements);
  }

  return SMESH_SUCCESS;
}

#else
int main() { return 1; }
#endif
