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
    ptrdiff_t n_local_nodes;
    ptrdiff_t n_global_nodes;
    mesh_coordinates_from_folder(comm->get(), mesh_folder, &spatial_dim,
                                 &points, &n_local_nodes, &n_global_nodes);

    count_t *n2eptr;
    element_idx_t *n2e_idx;

    if (!comm_rank) {
      create_directory(output_folder);
    }

    comm->barrier();

    SMESH_TRACE_SCOPE("pe2n+output (distributed)");
    create_n2e<idx_t, count_t, element_idx_t>(
        comm->get(), n_local_elements, n_global_elements, n_local_nodes,
        n_global_nodes, nnodesxelem, elems, &n2eptr, &n2e_idx);

    // Ensure it is always the same
    sort_n2e<count_t, element_idx_t>(n_local_nodes, n2eptr, n2e_idx);

    int *send_elements_displs = (int *)calloc(comm_size + 1, sizeof(int));
    int *send_elements_count = (int *)calloc(comm_size, sizeof(int));

    int *send_nodes_count = (int *)calloc(comm_size, sizeof(int));
    int *send_nodes_displs = (int *)calloc(comm_size + 1, sizeof(int));

    ptrdiff_t max_adj_count = 0;
    for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
      const count_t e_begin = n2eptr[i];
      const count_t e_end = n2eptr[i + 1];
      max_adj_count = std::max(max_adj_count, (ptrdiff_t)(e_end - e_begin));
    }

    const ptrdiff_t connected_ranks_capacity = std::max<ptrdiff_t>(1, max_adj_count);
    int *connected_ranks =
        (int *)malloc(static_cast<size_t>(connected_ranks_capacity) * sizeof(int));
    for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
      const count_t e_begin = n2eptr[i];
      const count_t e_end = n2eptr[i + 1];
      // if (e_end == e_begin) {
      //   continue;
      // }
      for (ptrdiff_t e = e_begin; e < e_end; ++e) {
        const element_idx_t element_idx = n2e_idx[e];
        const int element_owner =
            rank_owner(n_global_elements, element_idx, comm_size);
        connected_ranks[e - e_begin] = element_owner;
      }

      const size_t n_connected_ranks =
          sort_and_unique(connected_ranks, static_cast<size_t>(e_end - e_begin));

      for (size_t r = 0; r < n_connected_ranks; ++r) {
        send_nodes_displs[connected_ranks[r] + 1]++;
        send_elements_displs[connected_ranks[r] + 1] += e_end - e_begin;
      }
    }

    send_elements_displs[0] = 0;
    for (int r = 0; r < comm_size; r++) {
      send_elements_displs[r + 1] += send_elements_displs[r];
    }

    send_nodes_displs[0] = 0;
    for (int r = 0; r < comm_size; r++) {
      send_nodes_displs[r + 1] += send_nodes_displs[r];
    }

    const ptrdiff_t send_elements_size = send_elements_displs[comm_size];
    element_idx_t *send_elements =
        (element_idx_t *)malloc(send_elements_size * sizeof(element_idx_t));

    const ptrdiff_t send_nodes_size = send_nodes_displs[comm_size];
    idx_t *send_nodes = (idx_t *)malloc(send_nodes_size * sizeof(idx_t));
    count_t *send_n2e_count =
        (count_t *)malloc(send_nodes_size * sizeof(count_t));

    ptrdiff_t node_start = rank_start(n_global_nodes, comm_size, comm_rank);
    for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
      const count_t e_begin = n2eptr[i];
      const count_t e_end = n2eptr[i + 1];
      // if (e_end == e_begin) {
      //   continue;
      // }
      for (ptrdiff_t e = e_begin; e < e_end; ++e) {
        const element_idx_t element_idx = n2e_idx[e];
        const int element_owner =
            rank_owner(n_global_elements, element_idx, comm_size);
        connected_ranks[e - e_begin] = element_owner;
      }

      const size_t n_connected_ranks =
          sort_and_unique(connected_ranks, static_cast<size_t>(e_end - e_begin));

      for (size_t r = 0; r < n_connected_ranks; ++r) {
        const int cr = connected_ranks[r];
        const int node_pos = send_nodes_displs[cr] + send_nodes_count[cr];
        send_nodes[node_pos] = node_start + i;
        send_n2e_count[node_pos] = e_end - e_begin;

        send_nodes_count[cr]++;
        for (ptrdiff_t e = e_begin; e < e_end; ++e) {
          const element_idx_t element_idx = n2e_idx[e];
          send_elements[send_elements_displs[cr] + send_elements_count[cr]++] =
              element_idx;
        }
      }
    }

    // comm->print_callback([&](std::ostream &os) {
    //   for (int r = 0; r < comm_size; r++) {
    //     os << "send_nodes_count[" << r << "]: " << send_nodes_count[r] << "\n";
    //     os << "send_elements_count[" << r << "]: " << send_elements_count[r] << "\n";
    //   }

    //   for(ptrdiff_t i = 0; i < send_nodes_size; ++i) {
    //     os << "send_nodes[" << i << "]: " << send_nodes[i] << " (" << send_n2e_count[i] << ")\n";
    //   }

    //   for(ptrdiff_t i = 0; i < send_elements_size; ++i) {
    //     os << "send_elements[" << i << "]: " << send_elements[i] << "\n";
    //   }
    // });

    int *recv_nodes_count = (int *)calloc(comm_size, sizeof(int));
    int *recv_elements_count = (int *)calloc(comm_size, sizeof(int));
    MPI_Alltoall(send_nodes_count, 1, MPI_INT, recv_nodes_count, 1, MPI_INT, comm->get());

    MPI_Alltoall(send_elements_count, 1, MPI_INT, recv_elements_count, 1,
                 MPI_INT, comm->get());

    int *recv_nodes_displs = (int *)malloc((comm_size + 1) * sizeof(int));
    int *recv_elements_displs = (int *)malloc((comm_size + 1) * sizeof(int));
    recv_nodes_displs[0] = 0;
    recv_elements_displs[0] = 0;
    for (int r = 0; r < comm_size; r++) {
      recv_nodes_displs[r + 1] = recv_nodes_displs[r] + recv_nodes_count[r];
      recv_elements_displs[r + 1] =
          recv_elements_displs[r] + recv_elements_count[r];
    }

    const ptrdiff_t recv_nodes_size = recv_nodes_displs[comm_size];
    idx_t *recv_nodes = (idx_t *)malloc(recv_nodes_size * sizeof(idx_t));

    count_t *recv_n2e_ptr =
        (count_t *)malloc((recv_nodes_size + 1) * sizeof(count_t));

    const ptrdiff_t recv_elements_size = recv_elements_displs[comm_size];
    element_idx_t *recv_elements =
        (element_idx_t *)malloc(recv_elements_size * sizeof(element_idx_t));

    MPI_Alltoallv(send_nodes, send_nodes_count, send_nodes_displs,
                  smesh::mpi_type<idx_t>(), recv_nodes, recv_nodes_count,
                  recv_nodes_displs, smesh::mpi_type<idx_t>(), comm->get());

    recv_n2e_ptr[0] = 0;
    MPI_Alltoallv(send_n2e_count, send_nodes_count, send_nodes_displs,
                  smesh::mpi_type<count_t>(), &recv_n2e_ptr[1], recv_nodes_count,
                  recv_nodes_displs, smesh::mpi_type<count_t>(), comm->get());

    for (int r = 0; r < recv_nodes_size; r++) {
      recv_n2e_ptr[r + 1] += recv_n2e_ptr[r];
    }

    MPI_Alltoallv(send_elements, send_elements_count, send_elements_displs,
                  smesh::mpi_type<element_idx_t>(), recv_elements,
                  recv_elements_count, recv_elements_displs,
                  smesh::mpi_type<element_idx_t>(), comm->get());

    ptrdiff_t element_start =
        rank_start(n_global_elements, comm_size, comm_rank);
    idx_t **local_elements = (idx_t **)malloc(nnodesxelem * sizeof(idx_t *));
    for (int d = 0; d < nnodesxelem; ++d) {
      local_elements[d] = (idx_t *)malloc(n_local_elements * sizeof(idx_t));

      for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
        local_elements[d][i] = invalid_idx<idx_t>();
      }
    }

    for (ptrdiff_t i = 0; i < recv_nodes_size; ++i) {
      const count_t e_begin = recv_n2e_ptr[i];
      const count_t e_end = recv_n2e_ptr[i + 1];
      const idx_t node = recv_nodes[i];

      for (ptrdiff_t e = e_begin; e < e_end; ++e) {
        const element_idx_t element_idx = recv_elements[e];
        const int element_owner =
            rank_owner(n_global_elements, element_idx, comm_size);
        // printf("element_idx: %ld, element_owner: %d\n", (long)element_idx,
              //  element_owner);
        if (comm_rank == element_owner) {
          for (int d = 0; d < nnodesxelem; ++d) {
            if (node == elems[d][element_idx - element_start]) {
              local_elements[d][element_idx - element_start] =
                  i; // local node index
              break;
            }
          }
        }
      }
    }

    for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
      for (int d = 0; d < nnodesxelem; ++d) {
        const idx_t li = local_elements[d][i];
        if (li == invalid_idx<idx_t>() || li < 0 ||
            static_cast<ptrdiff_t>(li) >= recv_nodes_size) {
          SMESH_ERROR(
              "pe2n.exe: unmapped node for local element %ld (d=%d) on rank "
              "%d (li=%d recv_nodes_size=%ld)\n",
              (long)i, d, comm_rank, (int)li, (long)recv_nodes_size);
        }
      }
    }

    if (!comm_rank) {
      printf("#elements: %ld  #nodes: %ld\n", n_global_elements,
             n_global_nodes);
    }

    free(n2eptr);
    free(n2e_idx);
    free(elems);
    free(points);

    free(connected_ranks);
    free(send_elements_count);
    free(send_elements_displs);
    free(send_nodes_count);
    free(send_nodes_displs);
    free(send_elements);
    free(send_nodes);
    free(send_n2e_count);
    free(recv_nodes_count);
    free(recv_elements_count);
    free(recv_nodes_displs);
    free(recv_elements_displs);
    free(recv_nodes);
    free(recv_n2e_ptr);
    free(recv_elements);
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
