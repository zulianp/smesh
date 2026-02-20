#include "smesh_decompose.hpp"

#include "smesh_base.hpp"
#include "smesh_communicator.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_sort.hpp"
#include "smesh_tracer.hpp"
#include <algorithm>
#include <mpi.h>
#include <stddef.h>

namespace smesh {

template <typename idx_t>
static bool is_sorted(const idx_t *const SMESH_RESTRICT arr,
                      const ptrdiff_t n) {

  bool ret = true;
  for (ptrdiff_t i = 0; i < n - 1; ++i) {
    if (arr[i] > arr[i + 1]) {
      printf("arr[%ld] = %ld > arr[%ld] = %ld\n", (long)i, (long)arr[i],
             (long)(i + 1), (long)arr[i + 1]);
      ret = false;
    }
  }
  return ret;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(MPI_Comm comm, const ptrdiff_t n_local_elements,
               const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
               const ptrdiff_t n_global_nodes, const int nnodesxelem,
               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
               count_t **out_n2eptr, element_idx_t **out_n2e_idx) {
  SMESH_TRACE_SCOPE("create_n2e");
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int *send_displs = (int *)calloc(size + 1, sizeof(int));
  int *send_count = (int *)calloc(size, sizeof(int));
  int *recv_displs = (int *)calloc((size + 1), sizeof(int));
  int *recv_count = (int *)calloc(size, sizeof(int));

  const ptrdiff_t nodes_start = rank_start(n_global_nodes, size, rank);
  const ptrdiff_t elements_start = rank_start(n_global_elements, size, rank);

  SMESH_ASSERT(n_local_nodes == rank_split(n_global_nodes, size, rank));

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

  // TODO: prepare and use compressed version of e2n with e2n_ptr and e2n_idx
  // (this saves repeating the element index in send_elements and send_nodes)

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
      SMESH_ASSERT(recv_nodes[i] >= 0 && recv_nodes[i] < n_local_nodes);
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
  SMESH_TRACE_SCOPE("localize_element_indices");
  const ptrdiff_t elements_start =
      rank_start(n_global_elements, comm_size, comm_rank);

  SMESH_ASSERT(is_sorted(local2global, local2global_size));

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
          if (node == elems[d][element_idx - elements_start]) {
            local_elements[d][element_idx - elements_start] =
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
// Attention: it invalidates local_n2e and local_n2e_idx
template <typename idx_t, typename count_t, typename element_idx_t>
int rearrange_local_nodes(const int comm_size, const int comm_rank,
                          const ptrdiff_t n_global_elements,
                          const ptrdiff_t n_local_elements,
                          const int nnodesxelem,
                          const ptrdiff_t local2global_size,
                          count_t *const SMESH_RESTRICT local_n2e_ptr,
                          element_idx_t *const SMESH_RESTRICT local_n2e_idx,
                          idx_t *const SMESH_RESTRICT local2global,
                          idx_t **const SMESH_RESTRICT local_elements,
                          ptrdiff_t *const SMESH_RESTRICT out_n_owned,
                          ptrdiff_t *const SMESH_RESTRICT out_n_shared,
                          ptrdiff_t *const SMESH_RESTRICT out_n_ghosts) {
  SMESH_TRACE_SCOPE("rearrange_local_nodes");

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

  const ptrdiff_t n2e_nnz =
      static_cast<ptrdiff_t>(local_n2e_ptr[local2global_size]);

  // One reusable scratch buffer (memory-parsimonious):
  // - used as idx_t[] for copying local_elements[d] and local2global
  // - used as count_t[] + element_idx_t[] simultaneously for n2e reordering
  const size_t primary_bytes =
      static_cast<size_t>(std::max(local2global_size, n_local_elements)) *
      sizeof(idx_t);

  const size_t ptr_bytes =
      static_cast<size_t>(local2global_size + 1) * sizeof(count_t);

  auto align_up = [](size_t off, size_t alignment) -> size_t {
    return (off + alignment - 1) & ~(alignment - 1);
  };

  const size_t idx_off = align_up(ptr_bytes, alignof(element_idx_t));
  const size_t idx_bytes = static_cast<size_t>(n2e_nnz) * sizeof(element_idx_t);
  const size_t n2e_bytes = idx_off + idx_bytes;

  const size_t scratch_bytes = std::max(primary_bytes, n2e_bytes);
  char *scratch = (char *)malloc(scratch_bytes);
  idx_t *buff = (idx_t *)scratch;

  for (int d = 0; d < nnodesxelem; ++d) {
    memcpy(buff, local_elements[d], n_local_elements * sizeof(idx_t));
    for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
      local_elements[d][i] = index_map[buff[i]];
    }
  }

  // Rearrage n2e
  memcpy(buff, local2global, local2global_size * sizeof(idx_t));
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    local2global[index_map[i]] = buff[i];
  }

  count_t *temp_local_n2e_ptr = (count_t *)scratch;
  temp_local_n2e_ptr[0] = 0;
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    temp_local_n2e_ptr[index_map[i] + 1] =
        local_n2e_ptr[i + 1] - local_n2e_ptr[i];
  }

  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    temp_local_n2e_ptr[i + 1] += temp_local_n2e_ptr[i];
  }

  element_idx_t *idx_buff = (element_idx_t *)(scratch + idx_off);
  memcpy(idx_buff, local_n2e_idx, n2e_nnz * sizeof(element_idx_t));
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    ptrdiff_t from_begin = local_n2e_ptr[i];
    idx_t to = index_map[i];
    ptrdiff_t to_begin = temp_local_n2e_ptr[to];
    ptrdiff_t to_end = temp_local_n2e_ptr[to + 1];

    SMESH_ASSERT(to_end - to_begin == local_n2e_ptr[i + 1] - from_begin);

    memcpy(local_n2e_idx + to_begin, idx_buff + from_begin,
           (to_end - to_begin) * sizeof(element_idx_t));
  }

  memcpy(local_n2e_ptr, temp_local_n2e_ptr,
         (local2global_size + 1) * sizeof(count_t));

  free(scratch);
  free(index_map);

  *out_n_owned = n_owned;
  *out_n_shared = n_shared;
  *out_n_ghosts = local2global_size - n_owned;

  SMESH_ASSERT(is_sorted(local2global, n_owned - n_shared));
  SMESH_ASSERT(is_sorted(local2global + n_owned - n_shared, n_shared));
  SMESH_ASSERT(is_sorted(local2global + n_owned, local2global_size - n_owned));

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int rearrange_local_elements(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    const ptrdiff_t local2global_size,
    count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    idx_t **const SMESH_RESTRICT local_elements, const ptrdiff_t n_owned_nodes,
    ptrdiff_t *const SMESH_RESTRICT n_owned_not_shared,
    element_idx_t *const SMESH_RESTRICT element_local_to_global) {
  SMESH_TRACE_SCOPE("rearrange_local_elements");
  const ptrdiff_t elements_start =
      rank_start(n_global_elements, comm_size, comm_rank);
  idx_t *old_to_new_map = (idx_t *)malloc(n_local_elements * sizeof(idx_t));
  ptrdiff_t shared_count = 0;

  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    bool is_shared = false;
    for (int d = 0; d < nnodesxelem; ++d) {
      const idx_t node = local_elements[d][i];
      if (node >= n_owned_nodes) {
        is_shared = true;
        break;
      }
    }

    shared_count += is_shared ? 1 : 0;
  }

  ptrdiff_t shared_offset = n_local_elements - shared_count;
  ptrdiff_t local_offset = 0;
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    bool is_shared = false;
    for (int d = 0; d < nnodesxelem; ++d) {
      const idx_t node = local_elements[d][i];
      if (node >= n_owned_nodes) {
        is_shared = true;
        break;
      }
    }

    old_to_new_map[i] = is_shared ? shared_offset++ : local_offset++;
  }

  idx_t *buff = (idx_t *)malloc(n_local_elements * sizeof(idx_t));
  for (int d = 0; d < nnodesxelem; ++d) {
    memcpy(buff, local_elements[d], n_local_elements * sizeof(idx_t));
    for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
      local_elements[d][old_to_new_map[i]] = buff[i];
    }
  }

  for (ptrdiff_t i = 0; i < local_n2e_ptr[local2global_size]; ++i) {
    if (rank_owner(n_global_elements, local_n2e_idx[i], comm_size) !=
        comm_rank) {
      // Remote ones are unchanged as we do not yet know the new global indices
      continue;
    }
    local_n2e_idx[i] =
        elements_start + old_to_new_map[local_n2e_idx[i] - elements_start];
  }

  *n_owned_not_shared = n_local_elements - shared_count;
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    element_local_to_global[old_to_new_map[i]] = i + elements_start;
  }

  free(old_to_new_map);
  free(buff);
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int expand_aura_elements_inconsistent(
    MPI_Comm comm, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements,
    //  const ptrdiff_t elements_n_shared,
    const int nnodesxelem, count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const idx_t *const SMESH_RESTRICT local2global,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT local_elements,
    const element_idx_t *const SMESH_RESTRICT element_local_to_global,
    const ptrdiff_t node_n_owned, const ptrdiff_t nodes_n_ghosts,
    idx_t **const SMESH_RESTRICT out_aura_elements,
    idx_t **const SMESH_RESTRICT out_aura_element_nodes,
    ptrdiff_t *const SMESH_RESTRICT out_n_aura) {
  SMESH_TRACE_SCOPE("expand_aura_elements_inconsistent");
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  const ptrdiff_t elements_start =
      rank_start(n_global_elements, comm_size, comm_rank);

  const ptrdiff_t n_local_nodes = node_n_owned + nodes_n_ghosts;
  const ptrdiff_t n2e_nnz =
      static_cast<ptrdiff_t>(local_n2e_ptr[n_local_nodes]);

  element_idx_t *remote_elements = (element_idx_t *)malloc(
      static_cast<size_t>(n2e_nnz) * sizeof(element_idx_t));
  ptrdiff_t n_remote = 0;
  for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    for (ptrdiff_t e = e_begin; e < e_end; ++e) {
      const element_idx_t element_idx = local_n2e_idx[e];
      const int element_owner =
          rank_owner(n_global_elements, element_idx, comm_size);
      if (element_owner != comm_rank) {
        remote_elements[n_remote++] = element_idx;
      }
    }
  }

  const ptrdiff_t n_unique_remote =
      static_cast<ptrdiff_t>(sort_and_unique(remote_elements, n_remote));

  int *send_count = (int *)calloc(comm_size, sizeof(int));
  int *send_displs = (int *)calloc(comm_size + 1, sizeof(int));
  for (ptrdiff_t i = 0; i < n_unique_remote; ++i) {
    const int owner =
        rank_owner(n_global_elements, remote_elements[i], comm_size);
    send_displs[owner + 1]++;
  }

  for (int r = 0; r < comm_size; ++r) {
    send_displs[r + 1] += send_displs[r];
  }

  element_idx_t *send_elements = (element_idx_t *)malloc(
      static_cast<size_t>(send_displs[comm_size]) * sizeof(element_idx_t));

  for (ptrdiff_t i = 0; i < n_unique_remote; ++i) {
    const element_idx_t element_idx = remote_elements[i];
    const int owner = rank_owner(n_global_elements, element_idx, comm_size);
    send_elements[send_displs[owner] + send_count[owner]++] = element_idx;
  }

  free(remote_elements);

  int *recv_count = (int *)calloc(comm_size, sizeof(int));
  int *recv_displs = (int *)malloc((comm_size + 1) * sizeof(int));
  SMESH_MPI_CATCH(
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm));

  recv_displs[0] = 0;
  for (int r = 0; r < comm_size; ++r) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  const ptrdiff_t recv_size = recv_displs[comm_size];
  element_idx_t *recv_elements = (element_idx_t *)malloc(
      static_cast<size_t>(recv_size) * sizeof(element_idx_t));

  SMESH_MPI_CATCH(MPI_Alltoallv(send_elements, send_count, send_displs,
                                smesh::mpi_type<element_idx_t>(), recv_elements,
                                recv_count, recv_displs,
                                smesh::mpi_type<element_idx_t>(), comm));

  element_idx_t *old_to_new = (element_idx_t *)malloc(
      static_cast<size_t>(n_local_elements) * sizeof(element_idx_t));
  for (ptrdiff_t new_idx = 0; new_idx < n_local_elements; ++new_idx) {
    const element_idx_t old_global = element_local_to_global[new_idx];
    const ptrdiff_t old_off =
        static_cast<ptrdiff_t>(old_global - elements_start);
    old_to_new[old_off] = static_cast<element_idx_t>(new_idx);
  }

  for (int d = 0; d < nnodesxelem; ++d) {
    idx_t *send_nodes =
        (idx_t *)malloc(static_cast<size_t>(recv_size) * sizeof(idx_t));
    for (ptrdiff_t i = 0; i < recv_size; ++i) {
      const element_idx_t element_old = recv_elements[i];
      const ptrdiff_t old_off =
          static_cast<ptrdiff_t>(element_old - elements_start);
      const element_idx_t local_e = old_to_new[old_off];
      send_nodes[i] = local2global[local_elements[d][local_e]];
    }

    idx_t *recv_nodes = (idx_t *)malloc(
        static_cast<size_t>(send_displs[comm_size]) * sizeof(idx_t));
    SMESH_MPI_CATCH(MPI_Alltoallv(
        send_nodes, recv_count, recv_displs, smesh::mpi_type<idx_t>(),
        recv_nodes, send_count, send_displs, smesh::mpi_type<idx_t>(), comm));

    out_aura_element_nodes[d] = recv_nodes;
    free(send_nodes);
  }

  free(old_to_new);
  free(recv_elements);
  free(recv_count);
  free(recv_displs);

  *out_aura_elements = reinterpret_cast<idx_t *>(send_elements);
  *out_n_aura = static_cast<ptrdiff_t>(send_displs[comm_size]);

  free(send_count);
  free(send_displs);

  return SMESH_SUCCESS;
}

template <typename idx_t>
int prepare_node_renumbering(MPI_Comm comm, const ptrdiff_t n_global_nodes,
                             const ptrdiff_t owned_nodes_start,
                             const ptrdiff_t n_owned_nodes,
                             const idx_t *const SMESH_RESTRICT local2global,
                             idx_t *const SMESH_RESTRICT global2owned) {
  SMESH_TRACE_SCOPE("prepare_node_renumbering");
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  const ptrdiff_t nodes_start =
      rank_start(n_global_nodes, comm_size, comm_rank);

  int *send_nodes_count = (int *)calloc(comm_size, sizeof(int));
  int *send_nodes_displs = (int *)calloc(comm_size + 1, sizeof(int));
  int *recv_nodes_count = (int *)calloc(comm_size, sizeof(int));
  int *recv_nodes_displs = (int *)calloc(comm_size + 1, sizeof(int));

  for (ptrdiff_t i = 0; i < n_owned_nodes; ++i) {
    const int owner = rank_owner(n_global_nodes, local2global[i], comm_size);
    send_nodes_displs[owner + 1]++;
  }

  for (ptrdiff_t i = 0; i < comm_size; ++i) {
    send_nodes_displs[i + 1] += send_nodes_displs[i];
  }

  idx_t *send_nodes =
      (idx_t *)malloc(send_nodes_displs[comm_size] * sizeof(idx_t));
  idx_t *send_nodes_mapping =
      (idx_t *)malloc(send_nodes_displs[comm_size] * sizeof(idx_t));
  for (ptrdiff_t i = 0; i < n_owned_nodes; ++i) {
    const int owner = rank_owner(n_global_nodes, local2global[i], comm_size);
    send_nodes_mapping[send_nodes_displs[owner] + send_nodes_count[owner]] =
        local2global[i];
    send_nodes[send_nodes_displs[owner] + send_nodes_count[owner]] =
        owned_nodes_start + i;
    send_nodes_count[owner]++;
  }

  MPI_Alltoall(send_nodes_count, 1, MPI_INT, recv_nodes_count, 1, MPI_INT,
               comm);
  recv_nodes_displs[0] = 0;
  for (ptrdiff_t i = 0; i < comm_size; ++i) {
    recv_nodes_displs[i + 1] += recv_nodes_displs[i] + recv_nodes_count[i];
  }

  idx_t *recv_nodes =
      (idx_t *)malloc(recv_nodes_displs[comm_size] * sizeof(idx_t));
  MPI_Alltoallv(send_nodes, send_nodes_count, send_nodes_displs,
                smesh::mpi_type<idx_t>(), recv_nodes, recv_nodes_count,
                recv_nodes_displs, smesh::mpi_type<idx_t>(), comm);

  idx_t *recv_nodes_mapping =
      (idx_t *)malloc(recv_nodes_displs[comm_size] * sizeof(idx_t));
  MPI_Alltoallv(send_nodes_mapping, send_nodes_count, send_nodes_displs,
                smesh::mpi_type<idx_t>(), recv_nodes_mapping, recv_nodes_count,
                recv_nodes_displs, smesh::mpi_type<idx_t>(), comm);

  for (ptrdiff_t i = 0; i < recv_nodes_displs[comm_size]; ++i) {
    global2owned[recv_nodes_mapping[i] - nodes_start] = recv_nodes[i];
  }

  free(send_nodes_count);
  free(send_nodes_displs);
  free(recv_nodes_count);
  free(recv_nodes_displs);
  free(send_nodes);
  free(send_nodes_mapping);
  free(recv_nodes);
  free(recv_nodes_mapping);
  return SMESH_SUCCESS;
}

int node_ownership_ranges(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                          ptrdiff_t *const SMESH_RESTRICT owned_nodes_ranges) {
  SMESH_TRACE_SCOPE("node_ownership_ranges");
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  owned_nodes_ranges[0] = 0;
  MPI_Allgather(&n_owned_nodes, 1, mpi_type<ptrdiff_t>(),
                &owned_nodes_ranges[1], 1, mpi_type<ptrdiff_t>(), comm);
  for (ptrdiff_t i = 0; i < comm_size; ++i) {
    owned_nodes_ranges[i + 1] += owned_nodes_ranges[i];
  }

  return SMESH_SUCCESS;
}

// TODO: Stich together the aura elements and the local elements
// Unify ghost nodes and create local aura node index
// - we can identify ghost nodes with n2e graph
// - from e2n graph we can identify the unsorted aura nodes (old global indices)

template <typename idx_t>
int stitch_aura_elements(
    MPI_Comm comm, 
    const ptrdiff_t n_owned_nodes, 
    const ptrdiff_t n_shared_nodes,
    const ptrdiff_t n_ghost_nodes,
    const idx_t *const SMESH_RESTRICT local2global, const int nnodesxelem,
    const ptrdiff_t n_aura_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT e2n_aura,
    const ptrdiff_t n_local_elements, idx_t **const SMESH_RESTRICT e2n_local,
    idx_t **const SMESH_RESTRICT n2n_local2global_out,
    ptrdiff_t *const SMESH_RESTRICT out_n_aura_nodes) {
  SMESH_TRACE_SCOPE("stitch_aura_elements");
  (void)comm;
  const ptrdiff_t n_local_nodes = n_owned_nodes + n_ghost_nodes;

#ifndef NDEBUG
  // Communicator::wrap(comm)->print_callback([&](std::ostream &os) {
  //   os << "local2global: " << n_owned_nodes << "\n";
  //   for (ptrdiff_t i = 0; i < n_owned_nodes; ++i) {
  //     os << local2global[i] << " ";
  //   }
  //   os << "\n";
  // });

  SMESH_ASSERT(is_sorted(local2global, n_owned_nodes - n_shared_nodes));
  SMESH_ASSERT(is_sorted(local2global + n_owned_nodes - n_shared_nodes, n_shared_nodes));
  SMESH_ASSERT(is_sorted(local2global + n_owned_nodes, n_ghost_nodes));
#endif

  const idx_t *const owned_begin = local2global;
  const idx_t *const owned_end = local2global + n_owned_nodes - n_shared_nodes;
  const idx_t *const shared_begin = local2global + n_owned_nodes - n_shared_nodes;
  const idx_t *const shared_end = local2global + n_owned_nodes;
  const idx_t *const ghost_begin = local2global + n_owned_nodes;
  const idx_t *const ghost_end = local2global + n_local_nodes;

  const ptrdiff_t aura_flat_n = n_aura_elements * nnodesxelem;
  idx_t *aura_nodes = (idx_t *)malloc(
      static_cast<size_t>(std::max<ptrdiff_t>(1, aura_flat_n)) * sizeof(idx_t));
  for (ptrdiff_t i = 0; i < n_aura_elements; ++i) {
    for (int d = 0; d < nnodesxelem; ++d) {
      aura_nodes[i * nnodesxelem + d] = e2n_aura[d][i];
    }
  }

  const ptrdiff_t n_unique_aura_nodes =
      static_cast<ptrdiff_t>(sort_and_unique(aura_nodes, aura_flat_n));

  idx_t *aura_new = (idx_t *)malloc(
      static_cast<size_t>(std::max<ptrdiff_t>(1, n_unique_aura_nodes)) *
      sizeof(idx_t));

  ptrdiff_t n_aura_nodes = 0;

  for (ptrdiff_t i = 0; i < n_unique_aura_nodes; ++i) {
    const idx_t g = aura_nodes[i];
    auto it = std::lower_bound(owned_begin, owned_end, g);
    const bool is_owned = (it != owned_end && *it == g);
    if (is_owned) {
      continue;
    }

    it = std::lower_bound(shared_begin, shared_end, g);
    const bool is_shared = (it != shared_end && *it == g);
    if (is_shared) {
      continue;
    }

    it = std::lower_bound(ghost_begin, ghost_end, g);
    const bool is_ghost = (it != ghost_end && *it == g);
    if (!is_ghost) {
      aura_new[n_aura_nodes++] = g;
    }
  }

  idx_t *n2n_local2global = (idx_t *)malloc(
      static_cast<size_t>(n_local_nodes + n_aura_nodes) * sizeof(idx_t));
  memcpy(n2n_local2global, local2global, n_local_nodes * sizeof(idx_t));
  memcpy(&n2n_local2global[n_local_nodes], aura_new,
         n_aura_nodes * sizeof(idx_t));

  auto find_local = [&](const idx_t g) -> idx_t {
    auto it = std::lower_bound(owned_begin, owned_end, g);
    if (it != owned_end && *it == g) {
      return static_cast<idx_t>(it - local2global);
    }

    it = std::lower_bound(shared_begin, shared_end, g);
    if (it != shared_end && *it == g) {
      return static_cast<idx_t>(it - local2global);
    }

    it = std::lower_bound(ghost_begin, ghost_end, g);
    if (it != ghost_end && *it == g) {
      return static_cast<idx_t>(it - local2global);
    }

    auto it2 = std::lower_bound(aura_new, aura_new + n_aura_nodes, g);
    return static_cast<idx_t>(n_local_nodes + std::distance(aura_new, it2));
  };

  for (int d = 0; d < nnodesxelem; ++d) {
    for (ptrdiff_t i = 0; i < n_aura_elements; ++i) {
      e2n_aura[d][i] = find_local(e2n_aura[d][i]);
    }

    e2n_local[d] = (idx_t *)realloc(
        e2n_local[d], static_cast<size_t>(n_local_elements + n_aura_elements) *
                          sizeof(idx_t));
    memcpy(e2n_local[d] + n_local_elements, e2n_aura[d],
           static_cast<size_t>(n_aura_elements) * sizeof(idx_t));
  }

  free(aura_nodes);
  free(aura_new);

  *n2n_local2global_out = n2n_local2global;
  *out_n_aura_nodes = n_aura_nodes;
  return SMESH_SUCCESS;
}

template <typename idx_t>
int collect_ghost_and_aura_import_indices(
    MPI_Comm comm, const ptrdiff_t n_owned_nodes, const ptrdiff_t n_ghost_nodes,
    const ptrdiff_t n_aura_nodes, const ptrdiff_t n_global_nodes,
    const idx_t *const SMESH_RESTRICT local2global,
    const idx_t *const SMESH_RESTRICT global2owned,
    const ptrdiff_t *const SMESH_RESTRICT owned_node_ranges,
    idx_t *const SMESH_RESTRICT ghost_and_aura_to_owned) {
  SMESH_TRACE_SCOPE("collect_ghost_and_aura_import_indices");
  (void)owned_node_ranges;
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  const ptrdiff_t n_import = n_ghost_nodes + n_aura_nodes;
  const ptrdiff_t nodes_start =
      rank_start(n_global_nodes, comm_size, comm_rank);

  int *send_count = (int *)calloc(comm_size, sizeof(int));
  int *send_displs = (int *)calloc(comm_size + 1, sizeof(int));
  int *recv_count = (int *)calloc(comm_size, sizeof(int));
  int *recv_displs = (int *)calloc(comm_size + 1, sizeof(int));

  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const idx_t g = local2global[n_owned_nodes + i];
    const int owner = rank_owner(n_global_nodes, g, comm_size);
    send_displs[owner + 1]++;
  }

  for (int r = 0; r < comm_size; ++r) {
    send_displs[r + 1] += send_displs[r];
  }

  idx_t *send_nodes =
      (idx_t *)malloc(static_cast<size_t>(n_import) * sizeof(idx_t));
  idx_t *send_pos =
      (idx_t *)malloc(static_cast<size_t>(n_import) * sizeof(idx_t));
  int *cursor = (int *)calloc(comm_size, sizeof(int));

  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const idx_t g = local2global[n_owned_nodes + i];
    const int owner = rank_owner(n_global_nodes, g, comm_size);
    const int slot = send_displs[owner] + cursor[owner]++;
    send_nodes[slot] = g;
    send_pos[slot] = static_cast<idx_t>(i);
    send_count[owner]++;
  }

  SMESH_MPI_CATCH(
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm));

  recv_displs[0] = 0;
  for (int r = 0; r < comm_size; ++r) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }

  const ptrdiff_t recv_total = recv_displs[comm_size];
  idx_t *recv_nodes =
      (idx_t *)malloc(static_cast<size_t>(recv_total) * sizeof(idx_t));
  idx_t *recv_pos =
      (idx_t *)malloc(static_cast<size_t>(recv_total) * sizeof(idx_t));

  SMESH_MPI_CATCH(MPI_Alltoallv(
      send_nodes, send_count, send_displs, smesh::mpi_type<idx_t>(), recv_nodes,
      recv_count, recv_displs, smesh::mpi_type<idx_t>(), comm));
  SMESH_MPI_CATCH(MPI_Alltoallv(send_pos, send_count, send_displs,
                                smesh::mpi_type<idx_t>(), recv_pos, recv_count,
                                recv_displs, smesh::mpi_type<idx_t>(), comm));

  for (ptrdiff_t i = 0; i < recv_total; ++i) {
    recv_nodes[i] =
        global2owned[static_cast<ptrdiff_t>(recv_nodes[i] - nodes_start)];
  }

  SMESH_MPI_CATCH(MPI_Alltoallv(
      recv_nodes, recv_count, recv_displs, smesh::mpi_type<idx_t>(), send_nodes,
      send_count, send_displs, smesh::mpi_type<idx_t>(), comm));
  SMESH_MPI_CATCH(MPI_Alltoallv(recv_pos, recv_count, recv_displs,
                                smesh::mpi_type<idx_t>(), send_pos, send_count,
                                send_displs, smesh::mpi_type<idx_t>(), comm));

  for (ptrdiff_t i = 0; i < n_import; ++i) {
    ghost_and_aura_to_owned[send_pos[i]] = send_nodes[i];
  }

  free(send_count);
  free(send_displs);
  free(recv_count);
  free(recv_displs);
  free(send_nodes);
  free(send_pos);
  free(recv_nodes);
  free(recv_pos);
  free(cursor);
  return SMESH_SUCCESS;
}

// TODOs
// What we have: every rank as the aura elements with old global node indices
// Aura elements have both ghost and aura nodes (no owned), hence they can be
// used to create Rank to rank connectivity graphs with renumberd global indices
// and construct import/export lists for nodal quantities as well as for
// elemental quantities The steps are: 1) Renumber global nodes (which have more
// complex ownership structure, see node_ownership_ranges) 2) Renumber the nodes
// in the aura elements to the new global indices 3) Localize the aura elements
// nodes (ghosts already have local indices) 4) Append aura elements to the
// local_elements array 5) Create the import/export lists for nodal quantities
// (support ghost/aura only or combined) 6) Create the import/export lists for
// elemental quantities (support ghost/aura only or combined)

} // namespace smesh
