#include "smesh_decompose.hpp"

#include "smesh_base.hpp"
#include "smesh_communicator.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_sort.hpp"

#include <mpi.h>
#include <stddef.h>

namespace smesh {

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

template <typename idx_t, typename count_t, typename element_idx_t>
int redistribute_n2e(MPI_Comm comm, const int comm_size, const int comm_rank,
                     const ptrdiff_t n_local2global,
                     const ptrdiff_t n_global_nodes,
                     const ptrdiff_t n_global_elements,
                     const count_t *const SMESH_RESTRICT n2eptr,
                     const element_idx_t *const SMESH_RESTRICT n2e_idx,
                     ptrdiff_t *const SMESH_RESTRICT out_local2global_size,
                     idx_t **const SMESH_RESTRICT out_local2global,
                     count_t **const SMESH_RESTRICT out_local_n2e_ptr,
                     element_idx_t **const SMESH_RESTRICT out_local_n2e_idx) {
  int *send_elements_displs = (int *)calloc(comm_size + 1, sizeof(int));
  int *send_elements_count = (int *)calloc(comm_size, sizeof(int));

  int *send_nodes_count = (int *)calloc(comm_size, sizeof(int));
  int *send_nodes_displs = (int *)calloc(comm_size + 1, sizeof(int));

  ptrdiff_t max_adj_count = 0;
  for (ptrdiff_t i = 0; i < n_local2global; ++i) {
    const count_t e_begin = n2eptr[i];
    const count_t e_end = n2eptr[i + 1];
    max_adj_count = std::max(max_adj_count, (ptrdiff_t)(e_end - e_begin));
  }

  const ptrdiff_t connected_ranks_capacity =
      std::max<ptrdiff_t>(1, max_adj_count);
  int *connected_ranks = (int *)malloc(
      static_cast<size_t>(connected_ranks_capacity) * sizeof(int));
  for (ptrdiff_t i = 0; i < n_local2global; ++i) {
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
  for (ptrdiff_t i = 0; i < n_local2global; ++i) {
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

  int *recv_nodes_count = (int *)calloc(comm_size, sizeof(int));
  int *recv_elements_count = (int *)calloc(comm_size, sizeof(int));
  MPI_Alltoall(send_nodes_count, 1, MPI_INT, recv_nodes_count, 1, MPI_INT,
               comm);

  MPI_Alltoall(send_elements_count, 1, MPI_INT, recv_elements_count, 1, MPI_INT,
               comm);

  int *recv_nodes_displs = (int *)malloc((comm_size + 1) * sizeof(int));
  int *recv_elements_displs = (int *)malloc((comm_size + 1) * sizeof(int));
  recv_nodes_displs[0] = 0;
  recv_elements_displs[0] = 0;
  for (int r = 0; r < comm_size; r++) {
    recv_nodes_displs[r + 1] = recv_nodes_displs[r] + recv_nodes_count[r];
    recv_elements_displs[r + 1] =
        recv_elements_displs[r] + recv_elements_count[r];
  }

  const ptrdiff_t local2global_size = recv_nodes_displs[comm_size];
  idx_t *local2global = (idx_t *)malloc(local2global_size * sizeof(idx_t));

  count_t *local_n2e_ptr =
      (count_t *)malloc((local2global_size + 1) * sizeof(count_t));

  const ptrdiff_t local_n2e_size = recv_elements_displs[comm_size];
  element_idx_t *local_n2e_idx =
      (element_idx_t *)malloc(local_n2e_size * sizeof(element_idx_t));

  MPI_Alltoallv(send_nodes, send_nodes_count, send_nodes_displs,
                smesh::mpi_type<idx_t>(), local2global, recv_nodes_count,
                recv_nodes_displs, smesh::mpi_type<idx_t>(), comm);

  local_n2e_ptr[0] = 0;
  MPI_Alltoallv(send_n2e_count, send_nodes_count, send_nodes_displs,
                smesh::mpi_type<count_t>(), &local_n2e_ptr[1], recv_nodes_count,
                recv_nodes_displs, smesh::mpi_type<count_t>(), comm);

  for (int r = 0; r < local2global_size; r++) {
    local_n2e_ptr[r + 1] += local_n2e_ptr[r];
  }

  MPI_Alltoallv(send_elements, send_elements_count, send_elements_displs,
                smesh::mpi_type<element_idx_t>(), local_n2e_idx,
                recv_elements_count, recv_elements_displs,
                smesh::mpi_type<element_idx_t>(), comm);

  // Free internal send/recv buffers; keep only outputs.
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

  *out_local2global_size = local2global_size;
  *out_local2global = local2global;
  *out_local_n2e_ptr = local_n2e_ptr;
  *out_local_n2e_idx = local_n2e_idx;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int localize_element_indices(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    idx_t *const *const SMESH_RESTRICT elems, const ptrdiff_t local2global_size,
    const count_t *const SMESH_RESTRICT local_n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const idx_t *const SMESH_RESTRICT local2global,
    idx_t **const SMESH_RESTRICT local_elements) {
  const ptrdiff_t element_start =
      rank_start(n_global_elements, comm_size, comm_rank);

  for (int d = 0; d < nnodesxelem; ++d) {
    for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
      local_elements[d][i] = invalid_idx<idx_t>();
    }
  }

  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    const idx_t node = local2global[i];

    for (ptrdiff_t e = e_begin; e < e_end; ++e) {
      const element_idx_t element_idx = local_n2e_idx[e];
      const int element_owner =
          rank_owner(n_global_elements, element_idx, comm_size);
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
  return SMESH_SUCCESS;
}

// rearrange_local_nodes(...) rearrange node ordering based on the order:
// 1) owned, 2) shared, 3) ghosts and modify the local2global index accordingly
// Ownership is determine based on the smallest rank associated to the
// indicident element
//
template <typename idx_t, typename count_t, typename element_idx_t>
int rearrange_local_nodes(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    const ptrdiff_t local2global_size,
    const count_t *const SMESH_RESTRICT local_n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    idx_t *const SMESH_RESTRICT local2global,
    idx_t **const SMESH_RESTRICT local_elements) {

  ptrdiff_t n_owned = 0;
  ptrdiff_t n_shared = 0;

  for (ptrdiff_t i = 0; i < local2global_size; i++) {
    int owner = comm_size;
    int other = -1;
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    for (ptrdiff_t e = e_begin; e < e_end; ++e) {
      const element_idx_t element_idx = local_n2e_idx[e];
      const int element_owner =
          rank_owner(n_global_elements, element_idx, comm_size);

      owner = std::min(owner, element_owner);
      other = std::max(other, element_owner);
    }

    if (owner == comm_rank) {
      n_owned++;

      if (other != comm_rank) {
        n_shared++;
      }
    }
  }

  // const ptrdiff_t n_ghosts = local2global_size - n_owned;
  const ptrdiff_t n_owned_not_shared = n_owned - n_shared;

  idx_t *index_map = (idx_t *)malloc(local2global_size * sizeof(idx_t));

  ptrdiff_t count_ghosts = 0;
  ptrdiff_t count_owned_not_shared = 0;
  ptrdiff_t count_shared = 0;

  for (ptrdiff_t i = 0; i < local2global_size; i++) {
    int owner = comm_size;
    int other = -1;
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    for (ptrdiff_t e = e_begin; e < e_end; ++e) {
      const element_idx_t element_idx = local_n2e_idx[e];
      const int element_owner =
          rank_owner(n_global_elements, element_idx, comm_size);

      owner = std::min(owner, element_owner);
      other = std::max(other, element_owner);
    }

    if (owner == comm_rank) {
      if (other == comm_rank) {
        // Owned by this rank and not shared with others.
        index_map[i] = static_cast<idx_t>(count_owned_not_shared++);
      } else {
        // Owned by this rank but shared with others.
        index_map[i] = static_cast<idx_t>(n_owned_not_shared + count_shared++);
      }
    } else {
      // Not owned by this rank => ghost (includes nodes with no incident
      // elems).
      index_map[i] =
          static_cast<idx_t>(n_owned_not_shared + n_shared + count_ghosts++);
    }
  }

#ifndef NDEBUG
  SMESH_ASSERT(count_owned_not_shared == n_owned_not_shared);
  SMESH_ASSERT(count_shared == n_shared);
  SMESH_ASSERT(count_ghosts == (local2global_size - n_owned));
#endif

  const ptrdiff_t buff_max = std::max(local2global_size, n_local_elements);
  idx_t *buff = (idx_t *)malloc(buff_max * sizeof(idx_t));

  for (int d = 0; d < nnodesxelem; ++d) {
    memcpy(buff, local_elements[d], n_local_elements * sizeof(idx_t));
    for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
      local_elements[d][i] = index_map[buff[i]];
    }
  }

  memcpy(buff, local2global, local2global_size * sizeof(idx_t));
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    local2global[index_map[i]] = buff[i];
  }

  free(buff);
  free(index_map);
  return SMESH_SUCCESS;
}
} // namespace smesh
