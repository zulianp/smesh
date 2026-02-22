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
#include "smesh_distributed_aura.hpp"
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

    ptrdiff_t n_owned = 0;
    ptrdiff_t n_shared = 0;
    ptrdiff_t n_ghosts = 0;
    rearrange_local_nodes(comm_size, comm_rank, n_global_elements,
                          n_local_elements, nnodesxelem, local2global_size,
                          local_n2e_ptr, local_n2e_idx, local2global,
                          local_elements, &n_owned, &n_shared, &n_ghosts);

    element_idx_t *element_mapping =
        (element_idx_t *)malloc(n_local_elements * sizeof(element_idx_t));
    ptrdiff_t n_owned_not_shared = 0;
    rearrange_local_elements(comm_size, comm_rank, n_global_elements,
                             n_local_elements, nnodesxelem, local2global_size,
                             local_n2e_ptr, local_n2e_idx, local_elements,
                             n_owned, &n_owned_not_shared, element_mapping);

    idx_t *aura_elements = nullptr;
    idx_t **aura_element_nodes =
        (idx_t **)malloc(nnodesxelem * sizeof(idx_t *));
    for (int d = 0; d < nnodesxelem; ++d) {
      aura_element_nodes[d] = nullptr;
    }
    ptrdiff_t n_aura_elements = 0;
    expand_aura_elements_inconsistent(
        comm->get(), n_global_elements, n_local_elements, nnodesxelem,
        local_n2e_ptr, local_n2e_idx, local2global, local_elements,
        element_mapping, n_owned, n_ghosts, &aura_elements, aura_element_nodes,
        &n_aura_elements);

    long long owned_nodes_start_ll = 0;
    long long n_owned_ll = (long long)n_owned;
    SMESH_MPI_CATCH(MPI_Exscan(&n_owned_ll, &owned_nodes_start_ll, 1,
                               MPI_LONG_LONG, MPI_SUM, comm->get()));
    if (!comm_rank) {
      owned_nodes_start_ll = 0;
    }
    const ptrdiff_t owned_nodes_start =
        static_cast<ptrdiff_t>(owned_nodes_start_ll);

    idx_t *global2owned = (idx_t *)calloc(
        rank_split(n_global_nodes, comm_size, comm_rank), sizeof(idx_t));
    prepare_node_renumbering(comm->get(), n_global_nodes, owned_nodes_start,
                             n_owned, local2global, global2owned);

    ptrdiff_t *owned_node_ranges =
        (ptrdiff_t *)malloc((comm_size + 1) * sizeof(ptrdiff_t));
    node_ownership_ranges(comm->get(), n_owned, owned_node_ranges);

    idx_t *local2global_with_aura = nullptr;
    ptrdiff_t n_aura_nodes = 0;
    stitch_aura_elements(comm->get(), n_owned, n_shared, n_ghosts, local2global,
                         nnodesxelem, n_aura_elements, aura_element_nodes,
                         n_local_elements, local_elements,
                         &local2global_with_aura, &n_aura_nodes);
    free(local2global);
    local2global = local2global_with_aura;
    local2global_size = n_owned + n_ghosts + n_aura_nodes;

    SMESH_ASSERT(n_ghosts + n_aura_nodes > 0 || comm_size == 1);
    idx_t *ghost_and_aura_to_owned =
        (idx_t *)malloc((n_ghosts + n_aura_nodes) * sizeof(idx_t));
    collect_ghost_and_aura_import_indices(
        comm->get(), n_owned, n_ghosts, n_aura_nodes, n_global_nodes,
        local2global, global2owned, owned_node_ranges, ghost_and_aura_to_owned);

    node_ownership_ranges(comm->get(), n_owned, owned_node_ranges);
    int *owner =
        (int *)malloc((n_owned + n_ghosts + n_aura_nodes) * sizeof(int));
    determine_ownership(comm_size, comm_rank, n_owned, n_ghosts, n_aura_nodes,
                        ghost_and_aura_to_owned, owned_node_ranges, owner);

    group_ghost_and_aura_by_rank(comm_size, n_owned, n_ghosts, n_aura_nodes,
                                 local2global, ghost_and_aura_to_owned, owner,
                                 nnodesxelem, n_local_elements, n_aura_elements,
                                 local_elements);

    const ptrdiff_t n_local_nodes = n_owned + n_ghosts + n_aura_nodes;
    geom_t **local_points = (geom_t **)malloc(spatial_dim * sizeof(geom_t *));
    for (int d = 0; d < spatial_dim; ++d) {
      local_points[d] = (geom_t *)malloc(n_local_nodes * sizeof(geom_t));
      gather_mapped_field(comm->get(), n_local_nodes, n_global_nodes,
                          local2global, smesh::mpi_type<geom_t>(), points[d],
                          local_points[d]);
    }

    Path path_block = output_folder / std::to_string(comm_rank);
    create_directory(path_block);
    mesh_to_folder(path_block,
                   nnodesxelem == 8    ? HEX8
                   : nnodesxelem == 10 ? TET10
                                       : TET4,
                   n_local_elements + n_aura_elements, local_elements,
                   spatial_dim, n_local_nodes, local_points);

    array_write(path_block / "owner.int32", owner, n_local_nodes);

    { // FIXME find better solution for coding this
      int *send_count = (int *)malloc(comm_size * sizeof(int));
      int *send_displs = (int *)malloc((comm_size + 1) * sizeof(int));
      int *recv_count = (int *)malloc(comm_size * sizeof(int));
      int *recv_displs = (int *)malloc((comm_size + 1) * sizeof(int));
      idx_t *scatter_idx = nullptr;

      exchange_create<idx_t>(comm->get(), n_local_nodes, n_owned, owner,
                             owned_node_ranges, ghost_and_aura_to_owned,
                             send_count, send_displs, recv_count, recv_displs,
                             &scatter_idx);

      idx_t *owner_global = (idx_t *)malloc(n_local_nodes * sizeof(idx_t));
      for (ptrdiff_t i = 0; i < n_owned; ++i) {
        owner_global[i] = owned_node_ranges[comm_rank] + i;
      }

      idx_t *gather_buffer = (idx_t *)malloc(
          (recv_count[comm_size - 1] + recv_displs[comm_size - 1]) *
          sizeof(idx_t));

      exchange_gather(comm->get(), n_owned, recv_count, recv_displs, send_count,
                      send_displs, scatter_idx, owner_global, gather_buffer);

      array_write(path_block / "owner_global.int32", owner_global,
                  n_local_nodes);

      free(send_count);
      free(send_displs);
      free(recv_count);
      free(recv_displs);
      free(scatter_idx);
      free(gather_buffer);
      free(owner_global);
    }

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

    free(ghost_and_aura_to_owned);
    free(global2owned);
    free(owned_node_ranges);
    free(aura_elements);
    for (int d = 0; d < nnodesxelem; ++d) {
      free(aura_element_nodes[d]);
    }
    free(aura_element_nodes);

    free(elems);
    free(points);

    free(local2global);
    free(local_n2e_ptr);
    free(local_n2e_idx);
    free(element_mapping);
    for (int d = 0; d < nnodesxelem; ++d) {
      free(local_elements[d]);
    }
    free(local_elements);
    free(owner);

    for (int d = 0; d < spatial_dim; ++d) {
      free(local_points[d]);
    }
    free(local_points);
  }

  return SMESH_SUCCESS;
}

#else
int main() { return 1; }
#endif
