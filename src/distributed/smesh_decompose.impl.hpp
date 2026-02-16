#include "smesh_base.hpp"
#include "smesh_communicator.hpp"
#include "smesh_distributed_base.hpp"

#include <mpi.h>
#include <stddef.h>

namespace smesh {

inline ptrdiff_t rank_split(const ptrdiff_t n, const int comm_size,
                            const int comm_rank) {
  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t nlocal = uniform_split;
  ptrdiff_t remainder = n - nlocal * comm_size;

  if (remainder > comm_rank) {
    nlocal += 1;
  }

  return nlocal;
}

inline ptrdiff_t rank_start(const ptrdiff_t n, const int comm_size,
                            const int comm_rank) {
  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t remainder = n - uniform_split * comm_size;

  ptrdiff_t rank = comm_rank;
  ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);
  return rank_start;
}

inline ptrdiff_t rank_owner(const ptrdiff_t n, const ptrdiff_t gidx,
                            const int comm_size) {
  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t remainder = n - uniform_split * comm_size;

  ptrdiff_t rank = gidx / uniform_split;
  ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);

  if (gidx >= rank_start) {
#ifndef NDEBUG
    ptrdiff_t rank_end =
        rank_start + uniform_split + (ptrdiff_t)(rank < remainder);
    assert(gidx < rank_end);
#endif
    return rank;
  } else {
    rank -= 1;
#ifndef NDEBUG
    ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);
    ptrdiff_t rank_end =
        rank_start + uniform_split + (ptrdiff_t)(rank < remainder);

    assert(gidx >= rank_start);
    assert(gidx < rank_end);
#endif
    return rank;
  }
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(MPI_Comm comm, const ptrdiff_t n_local_elements,
               const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
               const ptrdiff_t n_global_nodes, const int nnodesxelem,
               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
               count_t **out_n2eptr, element_idx_t **out_n2e_idx) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int *send_displs = (int *)calloc(size + 1, sizeof(int));
  int *send_count = (int *)calloc(size, sizeof(int));
  int *recv_displs = (int *)calloc((size + 1), sizeof(int));
  int *recv_count = (int *)calloc(size, sizeof(int));

  const ptrdiff_t nodes_start = rank_start(n_global_nodes, size, rank);
  const ptrdiff_t elements_start = rank_start(n_global_elements, size, rank);

  assert(n_local_nodes == rank_split(n_global_nodes, size, rank));

  for (int d = 0; d < nnodesxelem; d++) {
    for (ptrdiff_t e = 0; e < n_local_elements; e++) {
      ptrdiff_t node = elements[d][e];
      int p = rank_owner(n_global_nodes, node, size);
      send_displs[p + 1]++;
    }
  }

  //   Get recv_displs values
  SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<int>(), recv_count,
                               1, mpi_type<int>(), comm));

  recv_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    recv_displs[r + 1] += recv_count[r] + recv_displs[r];
  }

  send_displs[0] = 0;
  for (int r = 0; r < size; r++) {
    send_displs[r + 1] += send_displs[r];
  }

  const ptrdiff_t send_size = send_displs[size];
  const ptrdiff_t recv_size = recv_displs[size];

  element_idx_t *send_elements =
      (element_idx_t *)malloc(send_size * sizeof(element_idx_t));

  idx_t *send_nodes = (idx_t *)malloc(send_size * sizeof(idx_t));

  element_idx_t *recv_elements =
      (element_idx_t *)malloc(recv_size * sizeof(element_idx_t));

  idx_t *recv_nodes = (idx_t *)malloc(recv_size * sizeof(idx_t));

  for (int d = 0; d < nnodesxelem; d++) {
    for (ptrdiff_t e = 0; e < n_local_elements; e++) {

      const ptrdiff_t node = elements[d][e];
      const int p = rank_owner(n_global_nodes, node, size);
      send_elements[send_displs[p] + send_count[p]] = elements_start + e;
      send_nodes[send_displs[p] + send_count[p]] = node;
      send_count[p]++;
    }
  }

  SMESH_MPI_CATCH(MPI_Alltoallv(
      send_elements, send_count, send_displs, mpi_type<element_idx_t>(),
      recv_elements, recv_count, recv_displs, mpi_type<element_idx_t>(), comm));

  SMESH_MPI_CATCH(MPI_Alltoallv(send_nodes, send_count, send_displs,
                                mpi_type<idx_t>(), recv_nodes, recv_count,
                                recv_displs, mpi_type<idx_t>(), comm));

  count_t *n2e_ptr = (count_t *)calloc((n_local_nodes + 1), sizeof(count_t));

  for (ptrdiff_t r = 0; r < size; r++) {
    int begin = recv_displs[r];
    int end = recv_displs[r + 1];
    for (int i = begin; i < end; i++) {
      recv_nodes[i] -= nodes_start;
      assert(recv_nodes[i] >= 0 && recv_nodes[i] < n_local_nodes);
      n2e_ptr[recv_nodes[i] + 1]++;
    }
  }

  for (ptrdiff_t i = 0; i < n_local_nodes; i++) {
    n2e_ptr[i + 1] += n2e_ptr[i];
  }

  element_idx_t *n2e_idx =
      (element_idx_t *)malloc(n2e_ptr[n_local_nodes] * sizeof(element_idx_t));
  count_t *book_keeping = (count_t *)calloc(n_local_nodes, sizeof(count_t));

  for (ptrdiff_t r = 0; r < size; r++) {
    int begin = recv_displs[r];
    int end = recv_displs[r + 1];
    for (int i = begin; i < end; i++) {
      n2e_idx[n2e_ptr[recv_nodes[i]] + book_keeping[recv_nodes[i]]++] =
          recv_elements[i];
    }
  }

  *out_n2eptr = n2e_ptr;
  *out_n2e_idx = n2e_idx;

  free(send_displs);
  free(send_count);
  free(recv_displs);
  free(recv_count);
  free(send_elements);
  free(send_nodes);
  free(recv_elements);
  free(recv_nodes);
  free(book_keeping);
  return SMESH_SUCCESS;
}

} // namespace smesh
