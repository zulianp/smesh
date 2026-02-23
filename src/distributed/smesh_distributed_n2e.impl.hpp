#ifndef SMESH_DISTRIBUTED_N2E_IMPL_HPP
#define SMESH_DISTRIBUTED_N2E_IMPL_HPP

#include "smesh_distributed_base.hpp"
#include "smesh_tracer.hpp"

#include <mpi.h>
#include <stddef.h>

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e_small(
    MPI_Comm comm, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
    const ptrdiff_t n_global_nodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    count_t **out_n2eptr, element_idx_t **out_n2e_idx) {
  SMESH_TRACE_SCOPE("create_n2e_small");
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

template <typename T>
int all_to_allv_64(T *send_elements, i64 *large_send_count,
                   i64 *large_send_displs, T *recv_elements,
                   i64 *large_recv_count, i64 *large_recv_displs, MPI_Comm comm,
                   i64 max_chunk_size) {
  SMESH_TRACE_SCOPE("all_to_allv_64");

  int size;
  MPI_Comm_size(comm, &size);

  i64 max_send_count = 0;
  i64 buffer_size = 0;

  for (int r = 0; r < size; r++) {
    max_send_count = std::max(max_send_count, large_send_displs[r + 1]);
    buffer_size += std::min(max_chunk_size, large_send_displs[r + 1]);
  }

  i64 overall_max_send_count = 0;
  MPI_Allreduce(&max_send_count, &overall_max_send_count, 1, mpi_type<i64>(),
                MPI_MAX, comm);

  const i64 n_rounds = div_round_up(overall_max_send_count, max_chunk_size);

  int *send_displs = (int *)calloc(size + 1, sizeof(int));
  int *send_count = (int *)calloc(size, sizeof(int));
  int *recv_displs = (int *)calloc((size + 1), sizeof(int));
  int *recv_count = (int *)calloc(size, sizeof(int));

  if (n_rounds <= 1) {
    // Copy to 32 bits versions
    for (int r = 0; r < size; r++) {
      send_count[r] = large_send_count[r];
      recv_count[r] = large_recv_count[r];

      send_displs[r + 1] = large_send_displs[r + 1];
      recv_displs[r + 1] = large_recv_displs[r + 1];
    }

    MPI_Alltoallv(send_elements, send_count, send_displs, mpi_type<T>(),
                  recv_elements, recv_count, recv_displs, mpi_type<T>(), comm);

  } else {

    i64 *send_book_keeping = (i64 *)calloc(size, sizeof(i64));
    i64 *recv_book_keeping = (i64 *)calloc(size, sizeof(i64));

    T *send_buffer = (T *)malloc(buffer_size * sizeof(T));
    T *recv_buffer = (T *)malloc(buffer_size * sizeof(T));

    // Chunking for MPI communication
    for (i64 round = 0; round < n_rounds; round++) {
      // Determine the chunk size for each rank
      for (int r = 0; r < size; r++) {
        i64 send_begin = large_send_displs[r] + send_book_keeping[r];
        i64 send_end =
            std::min(max_chunk_size, large_send_count[r] - send_begin);
        send_book_keeping[r] += send_end - send_begin;
        send_count[r] = send_end - send_begin;
        send_displs[r + 1] = send_displs[r] + send_count[r];

        if (send_count[r])
          memcpy(&send_buffer[send_displs[r]], &send_elements[send_begin],
                 send_count[r] * sizeof(T));

        i64 recv_begin = large_recv_displs[r] + recv_book_keeping[r];
        i64 recv_end =
            std::min(max_chunk_size, large_recv_count[r] - recv_begin);
        recv_count[r] = recv_end - recv_begin;
        recv_displs[r + 1] = recv_displs[r] + recv_count[r];
      }

      SMESH_MPI_CATCH(MPI_Alltoallv(send_elements, send_count, send_displs,
                                    mpi_type<T>(), recv_elements, recv_count,
                                    recv_displs, mpi_type<T>(), comm));

      for (int r = 0; r < size; r++) {
        i64 recv_begin = large_recv_displs[r] + recv_book_keeping[r];
        i64 recv_end =
            std::min(max_chunk_size, large_recv_count[r] - recv_begin);

        if (recv_count[r])
          memcpy(&recv_elements[recv_begin], &recv_buffer[recv_displs[r]],
                 recv_count[r] * sizeof(T));
        recv_book_keeping[r] += recv_end - recv_begin;
      }
    }

    free(send_buffer);
    free(recv_buffer);
    free(send_book_keeping);
    free(recv_book_keeping);
  }

  free(send_displs);
  free(send_count);
  free(recv_displs);
  free(recv_count);
  return SMESH_SUCCESS;
}

// Creates n2e for large meshes using chucking for MPI communication
template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e_large(
    MPI_Comm comm, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
    const ptrdiff_t n_global_nodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    count_t **out_n2eptr, element_idx_t **out_n2e_idx) {
  SMESH_TRACE_SCOPE("create_n2e_large");
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  i64 *send_displs = (i64 *)calloc(size + 1, sizeof(i64));
  i64 *send_count = (i64 *)calloc(size, sizeof(i64));
  i64 *recv_displs = (i64 *)calloc((size + 1), sizeof(i64));
  i64 *recv_count = (i64 *)calloc(size, sizeof(i64));

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
  SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count,
                               1, mpi_type<i64>(), comm));

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

  //   const i64 max_chunk_size = std::numeric_limits<i32>::max() / comm_size;
  const i64 max_chunk_size = 1024;

  SMESH_MPI_CATCH(all_to_allv_64(send_elements, send_count, send_displs,
                                 recv_elements, recv_count, recv_displs, comm,
                                 max_chunk_size));

  SMESH_MPI_CATCH(all_to_allv_64(send_nodes, send_count, send_displs,
                                 recv_nodes, recv_count, recv_displs, comm,
                                 max_chunk_size));

  count_t *n2e_ptr = (count_t *)calloc((n_local_nodes + 1), sizeof(count_t));

  for (ptrdiff_t r = 0; r < size; r++) {
    i64 begin = recv_displs[r];
    i64 end = recv_displs[r + 1];
    for (i64 i = begin; i < end; i++) {
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

} // namespace smesh

#endif // SMESH_DISTRIBUTED_N2E_IMPL_HPP