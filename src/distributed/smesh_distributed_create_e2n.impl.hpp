#ifndef SMESH_DISTRIBUTED_CREATE_E2N_IMPL_HPP
#define SMESH_DISTRIBUTED_CREATE_E2N_IMPL_HPP

#include "smesh_alltoallv.impl.hpp"
#include "smesh_alloc.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_aura.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_graph.hpp"
#include "smesh_sort.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <mpi.h>

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
static int create_n2e_from_e2n(
    MPI_Comm comm, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
    const ptrdiff_t n_global_nodes,
    const ptrdiff_t *const SMESH_RESTRICT e2n_ptr,
    const idx_t *const SMESH_RESTRICT e2n_idx, count_t **out_n2eptr,
    element_idx_t **out_n2e_idx) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  i64 *send_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
  i64 *send_count = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
  i64 *recv_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
  i64 *recv_count = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));

  const ptrdiff_t nodes_start = rank_start(n_global_nodes, size, rank);
  const ptrdiff_t elements_start = rank_start(n_global_elements, size, rank);
  SMESH_ASSERT(n_local_nodes == rank_split(n_global_nodes, size, rank));

  const ptrdiff_t nnz = n_local_elements > 0 ? e2n_ptr[n_local_elements] : 0;
  for (ptrdiff_t k = 0; k < nnz; ++k) {
    const int p = rank_owner(n_global_nodes, e2n_idx[k], size);
    send_displs[p + 1]++;
  }

  SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count,
                               1, mpi_type<i64>(), comm));

  recv_displs[0] = 0;
  for (int r = 0; r < size; ++r) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }
  send_displs[0] = 0;
  for (int r = 0; r < size; ++r) {
    send_displs[r + 1] += send_displs[r];
  }

  const ptrdiff_t send_size = send_displs[size];
  const ptrdiff_t recv_size = recv_displs[size];

  element_idx_t *send_elements =
      (element_idx_t *)SMESH_ALLOC((size_t)send_size * sizeof(element_idx_t));
  idx_t *send_nodes = (idx_t *)SMESH_ALLOC((size_t)send_size * sizeof(idx_t));
  element_idx_t *recv_elements =
      (element_idx_t *)SMESH_ALLOC((size_t)recv_size * sizeof(element_idx_t));
  idx_t *recv_nodes = (idx_t *)SMESH_ALLOC((size_t)recv_size * sizeof(idx_t));

  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    const ptrdiff_t k0 = e2n_ptr[e];
    const ptrdiff_t k1 = e2n_ptr[e + 1];
    const element_idx_t eid =
        static_cast<element_idx_t>(elements_start + e);
    for (ptrdiff_t k = k0; k < k1; ++k) {
      const idx_t node = e2n_idx[k];
      const int p = rank_owner(n_global_nodes, node, size);
      send_elements[send_displs[p] + send_count[p]] = eid;
      send_nodes[send_displs[p] + send_count[p]] = node;
      send_count[p]++;
    }
  }

  const i64 max_chunk_size = std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64(send_elements, send_count, send_displs,
                                 recv_elements, recv_count, recv_displs, comm,
                                 max_chunk_size));
  SMESH_MPI_CATCH(all_to_allv_64(send_nodes, send_count, send_displs, recv_nodes,
                                 recv_count, recv_displs, comm, max_chunk_size));

  count_t *n2e_ptr =
      (count_t *)SMESH_CALLOC((size_t)n_local_nodes + 1, sizeof(count_t));
  for (ptrdiff_t r = 0; r < size; ++r) {
    for (i64 i = recv_displs[r]; i < recv_displs[r + 1]; ++i) {
      recv_nodes[i] = static_cast<idx_t>(
          static_cast<ptrdiff_t>(recv_nodes[i]) - nodes_start);
      n2e_ptr[recv_nodes[i] + 1]++;
    }
  }
  for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
    n2e_ptr[i + 1] += n2e_ptr[i];
  }

  element_idx_t *n2e_idx = (element_idx_t *)SMESH_ALLOC(
      (size_t)n2e_ptr[n_local_nodes] * sizeof(element_idx_t));
  i64 *book_keeping = (i64 *)SMESH_CALLOC((size_t)n_local_nodes, sizeof(i64));
  for (ptrdiff_t r = 0; r < size; ++r) {
    for (i64 i = recv_displs[r]; i < recv_displs[r + 1]; ++i) {
      n2e_idx[n2e_ptr[recv_nodes[i]] + book_keeping[recv_nodes[i]]++] =
          recv_elements[i];
    }
  }

  *out_n2eptr = n2e_ptr;
  *out_n2e_idx = n2e_idx;
  SMESH_FREE(send_displs);
  SMESH_FREE(send_count);
  SMESH_FREE(recv_displs);
  SMESH_FREE(recv_count);
  SMESH_FREE(send_elements);
  SMESH_FREE(send_nodes);
  SMESH_FREE(recv_elements);
  SMESH_FREE(recv_nodes);
  SMESH_FREE(book_keeping);
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t,
          typename large_idx_t>
static int localize_e2n(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements,
    const ptrdiff_t *const SMESH_RESTRICT e2n_ptr,
    const idx_t *const SMESH_RESTRICT e2n_idx_global,
    const ptrdiff_t local2global_size,
    const count_t *const SMESH_RESTRICT local_n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const large_idx_t *const SMESH_RESTRICT local2global,
    idx_t *const SMESH_RESTRICT e2n_idx_local) {
  const ptrdiff_t elements_start =
      rank_start(n_global_elements, comm_size, comm_rank);
  const ptrdiff_t nnz =
      n_local_elements > 0 ? e2n_ptr[n_local_elements] : 0;
  for (ptrdiff_t k = 0; k < nnz; ++k) {
    e2n_idx_local[k] = invalid_idx<idx_t>();
  }

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    const large_idx_t node = local2global[i];
    for (count_t e = e_begin; e < e_end; ++e) {
      const element_idx_t element_idx = local_n2e_idx[e];
      if (rank_owner(n_global_elements, element_idx, comm_size) != comm_rank) {
        continue;
      }
      const ptrdiff_t slot =
          static_cast<ptrdiff_t>(element_idx) - elements_start;
      const ptrdiff_t k0 = e2n_ptr[slot];
      const ptrdiff_t k1 = e2n_ptr[slot + 1];
      for (ptrdiff_t k = k0; k < k1; ++k) {
        if (node == static_cast<large_idx_t>(e2n_idx_global[k])) {
          e2n_idx_local[k] = static_cast<idx_t>(i);
          break;
        }
      }
    }
  }
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t,
          typename large_idx_t>
static int rearrange_local_nodes_e2n(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const ptrdiff_t local2global_size,
    count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    large_idx_t *const SMESH_RESTRICT local2global,
    const ptrdiff_t *const SMESH_RESTRICT e2n_ptr,
    idx_t *const SMESH_RESTRICT e2n_idx, ptrdiff_t *const SMESH_RESTRICT out_n_owned,
    ptrdiff_t *const SMESH_RESTRICT out_n_shared,
    ptrdiff_t *const SMESH_RESTRICT out_n_ghosts) {
  ptrdiff_t n_owned = 0;
  ptrdiff_t n_shared = 0;

#pragma omp parallel for reduction(+ : n_owned, n_shared)
  for (ptrdiff_t i = 0; i < local2global_size; i++) {
    int owner = comm_size;
    int other = -1;
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    if (e_begin == e_end) {
      owner = comm_rank;
      other = comm_rank;
    }
    for (count_t e = e_begin; e < e_end; ++e) {
      const int element_owner =
          rank_owner(n_global_elements, local_n2e_idx[e], comm_size);
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

  const ptrdiff_t n_owned_not_shared = n_owned - n_shared;
  idx_t *index_map = (idx_t *)SMESH_ALLOC((size_t)local2global_size * sizeof(idx_t));
  ptrdiff_t count_ghosts = 0;
  ptrdiff_t count_owned_not_shared = 0;
  ptrdiff_t count_shared = 0;

  for (ptrdiff_t i = 0; i < local2global_size; i++) {
    int owner = comm_size;
    int other = -1;
    const count_t e_begin = local_n2e_ptr[i];
    const count_t e_end = local_n2e_ptr[i + 1];
    if (e_begin == e_end) {
      owner = comm_rank;
      other = comm_rank;
    }
    for (count_t e = e_begin; e < e_end; ++e) {
      const int element_owner =
          rank_owner(n_global_elements, local_n2e_idx[e], comm_size);
      owner = std::min(owner, element_owner);
      other = std::max(other, element_owner);
    }
    if (owner == comm_rank) {
      if (other == comm_rank) {
        index_map[i] = static_cast<idx_t>(count_owned_not_shared++);
      } else {
        index_map[i] =
            static_cast<idx_t>(n_owned_not_shared + count_shared++);
      }
    } else {
      index_map[i] = static_cast<idx_t>(n_owned_not_shared + n_shared +
                                        count_ghosts++);
    }
  }

  const ptrdiff_t nnz =
      n_local_elements > 0 ? e2n_ptr[n_local_elements] : 0;
  for (ptrdiff_t k = 0; k < nnz; ++k) {
    e2n_idx[k] = index_map[e2n_idx[k]];
  }

  large_idx_t *l2g_buff = (large_idx_t *)SMESH_ALLOC(
      (size_t)local2global_size * sizeof(large_idx_t));
  memcpy(l2g_buff, local2global, (size_t)local2global_size * sizeof(large_idx_t));
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    local2global[index_map[i]] = l2g_buff[i];
  }
  SMESH_FREE(l2g_buff);

  const ptrdiff_t n2e_nnz =
      static_cast<ptrdiff_t>(local_n2e_ptr[local2global_size]);
  count_t *temp_ptr =
      (count_t *)SMESH_ALLOC((size_t)(local2global_size + 1) * sizeof(count_t));
  element_idx_t *temp_idx =
      (element_idx_t *)SMESH_ALLOC((size_t)n2e_nnz * sizeof(element_idx_t));
  temp_ptr[0] = 0;
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    const ptrdiff_t dst = static_cast<ptrdiff_t>(index_map[i]);
    temp_ptr[dst + 1] = local_n2e_ptr[i + 1] - local_n2e_ptr[i];
  }
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    temp_ptr[i + 1] += temp_ptr[i];
  }
  for (ptrdiff_t i = 0; i < local2global_size; ++i) {
    const ptrdiff_t dst = static_cast<ptrdiff_t>(index_map[i]);
    const count_t n = local_n2e_ptr[i + 1] - local_n2e_ptr[i];
    memcpy(temp_idx + temp_ptr[dst], local_n2e_idx + local_n2e_ptr[i],
           (size_t)n * sizeof(element_idx_t));
  }
  memcpy(local_n2e_ptr, temp_ptr,
         (size_t)(local2global_size + 1) * sizeof(count_t));
  memcpy(local_n2e_idx, temp_idx, (size_t)n2e_nnz * sizeof(element_idx_t));
  SMESH_FREE(temp_ptr);
  SMESH_FREE(temp_idx);
  SMESH_FREE(index_map);

  *out_n_owned = n_owned;
  *out_n_shared = n_shared;
  *out_n_ghosts = local2global_size - n_owned;
  (void)e2n_ptr;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t,
          typename large_idx_t>
static int rearrange_local_elements_e2n(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const ptrdiff_t local2global_size,
    count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    ptrdiff_t *const SMESH_RESTRICT e2n_ptr, idx_t *const SMESH_RESTRICT e2n_idx,
    const ptrdiff_t n_owned_nodes,
    ptrdiff_t *const SMESH_RESTRICT n_owned_not_shared,
    large_idx_t *const SMESH_RESTRICT element_local_to_global) {
  const ptrdiff_t elements_start =
      rank_start(n_global_elements, comm_size, comm_rank);
  idx_t *old_to_new = (idx_t *)SMESH_ALLOC((size_t)n_local_elements * sizeof(idx_t));
  ptrdiff_t shared_count = 0;
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    bool is_shared = false;
    for (ptrdiff_t k = e2n_ptr[i]; k < e2n_ptr[i + 1]; ++k) {
      if (e2n_idx[k] >= n_owned_nodes) {
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
    for (ptrdiff_t k = e2n_ptr[i]; k < e2n_ptr[i + 1]; ++k) {
      if (e2n_idx[k] >= n_owned_nodes) {
        is_shared = true;
        break;
      }
    }
    old_to_new[i] = is_shared ? static_cast<idx_t>(shared_offset++)
                              : static_cast<idx_t>(local_offset++);
  }

  const ptrdiff_t nnz = n_local_elements > 0 ? e2n_ptr[n_local_elements] : 0;
  ptrdiff_t *new_ptr =
      (ptrdiff_t *)SMESH_ALLOC((size_t)(n_local_elements + 1) * sizeof(ptrdiff_t));
  idx_t *new_idx = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(nnz, 1) *
                                        sizeof(idx_t));
  new_ptr[0] = 0;
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    new_ptr[old_to_new[i] + 1] = e2n_ptr[i + 1] - e2n_ptr[i];
  }
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    new_ptr[i + 1] += new_ptr[i];
  }
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    const ptrdiff_t dst = static_cast<ptrdiff_t>(old_to_new[i]);
    const ptrdiff_t n = e2n_ptr[i + 1] - e2n_ptr[i];
    memcpy(new_idx + new_ptr[dst], e2n_idx + e2n_ptr[i],
           (size_t)n * sizeof(idx_t));
  }
  memcpy(e2n_ptr, new_ptr, (size_t)(n_local_elements + 1) * sizeof(ptrdiff_t));
  if (nnz > 0) {
    memcpy(e2n_idx, new_idx, (size_t)nnz * sizeof(idx_t));
  }
  SMESH_FREE(new_ptr);
  SMESH_FREE(new_idx);

  for (ptrdiff_t i = 0; i < local_n2e_ptr[local2global_size]; ++i) {
    if (rank_owner(n_global_elements, local_n2e_idx[i], comm_size) !=
        comm_rank) {
      continue;
    }
    local_n2e_idx[i] = static_cast<element_idx_t>(
        elements_start +
        old_to_new[static_cast<ptrdiff_t>(local_n2e_idx[i]) - elements_start]);
  }

  *n_owned_not_shared = n_local_elements - shared_count;
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    element_local_to_global[old_to_new[i]] =
        static_cast<large_idx_t>(i + elements_start);
  }
  SMESH_FREE(old_to_new);
  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t,
          typename large_idx_t>
static int expand_aura_elements_e2n(
    MPI_Comm comm, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements,
    const ptrdiff_t *const SMESH_RESTRICT e2n_ptr,
    const idx_t *const SMESH_RESTRICT e2n_idx,
    const count_t *const SMESH_RESTRICT local_n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const large_idx_t *const SMESH_RESTRICT local2global,
    const large_idx_t *const SMESH_RESTRICT element_local_to_global,
    const ptrdiff_t node_n_owned, large_idx_t **const SMESH_RESTRICT out_aura_mapping,
    ptrdiff_t **const SMESH_RESTRICT out_aura_e2n_ptr,
    idx_t **const SMESH_RESTRICT out_aura_e2n_idx,
    ptrdiff_t *const SMESH_RESTRICT out_n_aura) {
  int comm_rank = 0;
  int comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  const ptrdiff_t elements_start =
      rank_start(n_global_elements, comm_size, comm_rank);

  const ptrdiff_t n2e_nnz =
      static_cast<ptrdiff_t>(local_n2e_ptr[node_n_owned]);
  element_idx_t *remote_elements = (element_idx_t *)SMESH_ALLOC(
      (size_t)std::max<ptrdiff_t>(n2e_nnz, 1) * sizeof(element_idx_t));
  ptrdiff_t n_remote = 0;
  for (ptrdiff_t i = 0; i < node_n_owned; ++i) {
    for (count_t e = local_n2e_ptr[i]; e < local_n2e_ptr[i + 1]; ++e) {
      const element_idx_t element_idx = local_n2e_idx[e];
      if (rank_owner(n_global_elements, element_idx, comm_size) != comm_rank) {
        remote_elements[n_remote++] = element_idx;
      }
    }
  }
  const ptrdiff_t n_unique_remote =
      static_cast<ptrdiff_t>(sort_and_unique(remote_elements, (size_t)n_remote));

  i64 *send_count = (i64 *)SMESH_CALLOC((size_t)comm_size, sizeof(i64));
  i64 *send_displs = (i64 *)SMESH_CALLOC((size_t)comm_size + 1, sizeof(i64));
  for (ptrdiff_t i = 0; i < n_unique_remote; ++i) {
    send_displs[rank_owner(n_global_elements, remote_elements[i], comm_size) +
                1]++;
  }
  for (int r = 0; r < comm_size; ++r) {
    send_displs[r + 1] += send_displs[r];
  }
  element_idx_t *send_elements = (element_idx_t *)SMESH_ALLOC(
      (size_t)std::max<i64>(send_displs[comm_size], 1) * sizeof(element_idx_t));
  for (ptrdiff_t i = 0; i < n_unique_remote; ++i) {
    const int owner =
        rank_owner(n_global_elements, remote_elements[i], comm_size);
    send_elements[send_displs[owner] + send_count[owner]++] =
        remote_elements[i];
  }
  SMESH_FREE(remote_elements);

  i64 *recv_count = (i64 *)SMESH_CALLOC((size_t)comm_size, sizeof(i64));
  i64 *recv_displs = (i64 *)SMESH_ALLOC(((size_t)comm_size + 1) * sizeof(i64));
  SMESH_MPI_CATCH(MPI_Alltoall(send_count, 1, mpi_type<i64>(), recv_count, 1,
                               mpi_type<i64>(), comm));
  recv_displs[0] = 0;
  for (int r = 0; r < comm_size; ++r) {
    recv_displs[r + 1] = recv_displs[r] + recv_count[r];
  }
  const i64 recv_size = recv_displs[comm_size];
  element_idx_t *recv_elements = (element_idx_t *)SMESH_ALLOC(
      (size_t)std::max<i64>(recv_size, 1) * sizeof(element_idx_t));
  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / comm_size;
  SMESH_MPI_CATCH(all_to_allv_64(send_elements, send_count, send_displs,
                                 recv_elements, recv_count, recv_displs, comm,
                                 max_chunk_size));

  element_idx_t *old_to_new = (element_idx_t *)SMESH_ALLOC(
      (size_t)n_local_elements * sizeof(element_idx_t));
  for (ptrdiff_t new_idx = 0; new_idx < n_local_elements; ++new_idx) {
    const ptrdiff_t old_off =
        static_cast<ptrdiff_t>(element_local_to_global[new_idx] -
                               elements_start);
    old_to_new[old_off] = static_cast<element_idx_t>(new_idx);
  }

  i32 *send_nxe = (i32 *)SMESH_ALLOC((size_t)std::max<i64>(recv_size, 1) *
                                     sizeof(i32));
  i64 *node_send_count = (i64 *)SMESH_CALLOC((size_t)comm_size, sizeof(i64));
  for (i64 i = 0; i < recv_size; ++i) {
    const ptrdiff_t old_off =
        static_cast<ptrdiff_t>(recv_elements[i]) - elements_start;
    const ptrdiff_t local_e = static_cast<ptrdiff_t>(old_to_new[old_off]);
    const i32 nxe = static_cast<i32>(e2n_ptr[local_e + 1] - e2n_ptr[local_e]);
    send_nxe[i] = nxe;
  }
  i64 cursor = 0;
  for (int r = 0; r < comm_size; ++r) {
    for (i64 i = 0; i < recv_count[r]; ++i) {
      node_send_count[r] += send_nxe[cursor++];
    }
  }
  i64 *node_send_displs =
      (i64 *)SMESH_ALLOC(((size_t)comm_size + 1) * sizeof(i64));
  node_send_displs[0] = 0;
  for (int r = 0; r < comm_size; ++r) {
    node_send_displs[r + 1] = node_send_displs[r] + node_send_count[r];
  }
  idx_t *send_nodes = (idx_t *)SMESH_ALLOC(
      (size_t)std::max<i64>(node_send_displs[comm_size], 1) * sizeof(idx_t));
  ptrdiff_t w = 0;
  for (i64 i = 0; i < recv_size; ++i) {
    const ptrdiff_t old_off =
        static_cast<ptrdiff_t>(recv_elements[i]) - elements_start;
    const ptrdiff_t local_e = static_cast<ptrdiff_t>(old_to_new[old_off]);
    for (ptrdiff_t k = e2n_ptr[local_e]; k < e2n_ptr[local_e + 1]; ++k) {
      send_nodes[w++] = static_cast<idx_t>(local2global[e2n_idx[k]]);
    }
  }
  SMESH_FREE(old_to_new);
  SMESH_FREE(recv_elements);

  i32 *recv_nxe = (i32 *)SMESH_ALLOC(
      (size_t)std::max<i64>(send_displs[comm_size], 1) * sizeof(i32));
  SMESH_MPI_CATCH(all_to_allv_64(send_nxe, recv_count, recv_displs, recv_nxe,
                                 send_count, send_displs, comm, max_chunk_size));
  SMESH_FREE(send_nxe);

  i64 *node_recv_count = (i64 *)SMESH_CALLOC((size_t)comm_size, sizeof(i64));
  cursor = 0;
  for (int r = 0; r < comm_size; ++r) {
    for (i64 i = 0; i < send_count[r]; ++i) {
      node_recv_count[r] += recv_nxe[cursor++];
    }
  }
  i64 *node_recv_displs =
      (i64 *)SMESH_ALLOC(((size_t)comm_size + 1) * sizeof(i64));
  node_recv_displs[0] = 0;
  for (int r = 0; r < comm_size; ++r) {
    node_recv_displs[r + 1] = node_recv_displs[r] + node_recv_count[r];
  }
  idx_t *recv_nodes = (idx_t *)SMESH_ALLOC(
      (size_t)std::max<i64>(node_recv_displs[comm_size], 1) * sizeof(idx_t));
  SMESH_MPI_CATCH(all_to_allv_64(send_nodes, node_send_count, node_send_displs,
                                 recv_nodes, node_recv_count, node_recv_displs,
                                 comm, max_chunk_size));
  SMESH_FREE(send_nodes);
  SMESH_FREE(node_send_count);
  SMESH_FREE(node_send_displs);
  SMESH_FREE(node_recv_count);
  SMESH_FREE(node_recv_displs);
  SMESH_FREE(recv_count);
  SMESH_FREE(recv_displs);

  const ptrdiff_t n_aura = static_cast<ptrdiff_t>(send_displs[comm_size]);
  ptrdiff_t *aura_ptr =
      (ptrdiff_t *)SMESH_ALLOC((size_t)(n_aura + 1) * sizeof(ptrdiff_t));
  aura_ptr[0] = 0;
  for (ptrdiff_t i = 0; i < n_aura; ++i) {
    aura_ptr[i + 1] = aura_ptr[i] + recv_nxe[i];
  }
  SMESH_FREE(recv_nxe);
  const ptrdiff_t aura_nnz = aura_ptr[n_aura];
  idx_t *aura_idx = recv_nodes;
  if (aura_nnz == 0) {
    SMESH_FREE(recv_nodes);
    aura_idx = (idx_t *)SMESH_ALLOC(sizeof(idx_t));
  }

  large_idx_t *aura_mapping = (large_idx_t *)SMESH_ALLOC(
      (size_t)std::max<ptrdiff_t>(n_aura, 1) * sizeof(large_idx_t));
  for (ptrdiff_t i = 0; i < n_aura; ++i) {
    aura_mapping[i] = static_cast<large_idx_t>(send_elements[i]);
  }
  SMESH_FREE(send_elements);
  SMESH_FREE(send_count);
  SMESH_FREE(send_displs);

  *out_aura_mapping = aura_mapping;
  *out_aura_e2n_ptr = aura_ptr;
  *out_aura_e2n_idx = aura_idx;
  *out_n_aura = n_aura;
  (void)e2n_idx;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename large_idx_t>
static int stich_aura_elements_e2n(
    const ptrdiff_t n_owned_nodes, const ptrdiff_t n_shared_nodes,
    const ptrdiff_t n_ghost_nodes,
    const large_idx_t *const SMESH_RESTRICT local2global,
    const ptrdiff_t n_local_elements, ptrdiff_t **const SMESH_RESTRICT e2n_ptr_io,
    idx_t **const SMESH_RESTRICT e2n_idx_io, const ptrdiff_t n_aura_elements,
    const ptrdiff_t *const SMESH_RESTRICT aura_e2n_ptr,
    idx_t *const SMESH_RESTRICT aura_e2n_idx,
    large_idx_t **const SMESH_RESTRICT n2n_local2global_out,
    ptrdiff_t *const SMESH_RESTRICT out_n_aura_nodes) {
  const ptrdiff_t n_local_nodes = n_owned_nodes + n_ghost_nodes;
  const large_idx_t *const owned_begin = local2global;
  const large_idx_t *const owned_end =
      local2global + n_owned_nodes - n_shared_nodes;
  const large_idx_t *const shared_begin =
      local2global + n_owned_nodes - n_shared_nodes;
  const large_idx_t *const shared_end = local2global + n_owned_nodes;
  const large_idx_t *const ghost_begin = local2global + n_owned_nodes;
  const large_idx_t *const ghost_end = local2global + n_local_nodes;

  const ptrdiff_t aura_nnz =
      n_aura_elements > 0 ? aura_e2n_ptr[n_aura_elements] : 0;
  large_idx_t *aura_nodes = (large_idx_t *)SMESH_ALLOC(
      (size_t)std::max<ptrdiff_t>(aura_nnz, 1) * sizeof(large_idx_t));
  for (ptrdiff_t k = 0; k < aura_nnz; ++k) {
    aura_nodes[k] = static_cast<large_idx_t>(aura_e2n_idx[k]);
  }
  const ptrdiff_t n_unique =
      static_cast<ptrdiff_t>(sort_and_unique(aura_nodes, (size_t)aura_nnz));

  large_idx_t *aura_new = (large_idx_t *)SMESH_ALLOC(
      (size_t)std::max<ptrdiff_t>(n_unique, 1) * sizeof(large_idx_t));
  ptrdiff_t n_aura_nodes = 0;
  for (ptrdiff_t i = 0; i < n_unique; ++i) {
    const large_idx_t g = aura_nodes[i];
    auto it = std::lower_bound(owned_begin, owned_end, g);
    if (it != owned_end && *it == g) {
      continue;
    }
    it = std::lower_bound(shared_begin, shared_end, g);
    if (it != shared_end && *it == g) {
      continue;
    }
    it = std::lower_bound(ghost_begin, ghost_end, g);
    if (it != ghost_end && *it == g) {
      continue;
    }
    aura_new[n_aura_nodes++] = g;
  }

  large_idx_t *n2n_local2global = (large_idx_t *)SMESH_ALLOC(
      (size_t)(n_local_nodes + n_aura_nodes) * sizeof(large_idx_t));
  memcpy(n2n_local2global, local2global,
         (size_t)n_local_nodes * sizeof(large_idx_t));
  memcpy(n2n_local2global + n_local_nodes, aura_new,
         (size_t)n_aura_nodes * sizeof(large_idx_t));

  auto find_local = [&](const large_idx_t g) -> idx_t {
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

  for (ptrdiff_t k = 0; k < aura_nnz; ++k) {
    aura_e2n_idx[k] = find_local(static_cast<large_idx_t>(aura_e2n_idx[k]));
  }

  const ptrdiff_t local_nnz =
      n_local_elements > 0 ? (*e2n_ptr_io)[n_local_elements] : 0;
  const ptrdiff_t n_total = n_local_elements + n_aura_elements;
  ptrdiff_t *new_ptr =
      (ptrdiff_t *)SMESH_ALLOC((size_t)(n_total + 1) * sizeof(ptrdiff_t));
  memcpy(new_ptr, *e2n_ptr_io,
         (size_t)(n_local_elements + 1) * sizeof(ptrdiff_t));
  for (ptrdiff_t i = 0; i < n_aura_elements; ++i) {
    new_ptr[n_local_elements + i + 1] =
        new_ptr[n_local_elements] + aura_e2n_ptr[i + 1];
  }
  idx_t *new_idx = (idx_t *)SMESH_ALLOC(
      (size_t)std::max<ptrdiff_t>(local_nnz + aura_nnz, 1) * sizeof(idx_t));
  if (local_nnz > 0) {
    memcpy(new_idx, *e2n_idx_io, (size_t)local_nnz * sizeof(idx_t));
  }
  if (aura_nnz > 0) {
    memcpy(new_idx + local_nnz, aura_e2n_idx, (size_t)aura_nnz * sizeof(idx_t));
  }
  SMESH_FREE(*e2n_ptr_io);
  SMESH_FREE(*e2n_idx_io);
  *e2n_ptr_io = new_ptr;
  *e2n_idx_io = new_idx;

  SMESH_FREE(aura_nodes);
  SMESH_FREE(aura_new);
  *n2n_local2global_out = n2n_local2global;
  *out_n_aura_nodes = n_aura_nodes;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename mapping_t>
static int group_ghost_and_aura_e2n(
    const int comm_size, const ptrdiff_t n_owned, const ptrdiff_t n_ghosts,
    const ptrdiff_t n_aura_nodes, mapping_t *const SMESH_RESTRICT local2global,
    idx_t *const SMESH_RESTRICT ghost_and_aura_to_owned,
    int *const SMESH_RESTRICT owner, const ptrdiff_t n_total_elements,
    const ptrdiff_t *const SMESH_RESTRICT e2n_ptr,
    idx_t *const SMESH_RESTRICT e2n_idx) {
  const ptrdiff_t n_import = n_ghosts + n_aura_nodes;
  if (comm_size <= 1 || n_import <= 0) {
    return SMESH_SUCCESS;
  }

  ptrdiff_t *ghost_displs =
      (ptrdiff_t *)SMESH_CALLOC((size_t)comm_size + 1, sizeof(ptrdiff_t));
  ptrdiff_t *aura_displs =
      (ptrdiff_t *)SMESH_CALLOC((size_t)comm_size + 1, sizeof(ptrdiff_t));
  ptrdiff_t *ghost_cursor =
      (ptrdiff_t *)SMESH_CALLOC((size_t)comm_size, sizeof(ptrdiff_t));
  ptrdiff_t *aura_cursor =
      (ptrdiff_t *)SMESH_CALLOC((size_t)comm_size, sizeof(ptrdiff_t));

  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const int r = owner[n_owned + i];
    if (i < n_ghosts) {
      ghost_displs[(size_t)r + 1]++;
    } else {
      aura_displs[(size_t)r + 1]++;
    }
  }
  for (int r = 0; r < comm_size; ++r) {
    ghost_displs[(size_t)r + 1] += ghost_displs[(size_t)r];
    aura_displs[(size_t)r + 1] += aura_displs[(size_t)r];
  }
  for (int r = 0; r <= comm_size; ++r) {
    aura_displs[(size_t)r] += n_ghosts;
  }

  idx_t *old_to_new = (idx_t *)SMESH_ALLOC((size_t)n_import * sizeof(idx_t));
  idx_t *ghost_tmp = (idx_t *)SMESH_ALLOC((size_t)n_import * sizeof(idx_t));
  mapping_t *l2g_tmp =
      (mapping_t *)SMESH_ALLOC((size_t)n_import * sizeof(mapping_t));
  int *owner_tmp = (int *)SMESH_ALLOC((size_t)n_import * sizeof(int));

  for (ptrdiff_t i = 0; i < n_import; ++i) {
    const int r = owner[n_owned + i];
    const ptrdiff_t pos =
        (i < n_ghosts) ? (ghost_displs[(size_t)r] + ghost_cursor[(size_t)r]++)
                       : (aura_displs[(size_t)r] + aura_cursor[(size_t)r]++);
    old_to_new[i] = static_cast<idx_t>(pos);
    ghost_tmp[pos] = ghost_and_aura_to_owned[i];
    l2g_tmp[pos] = local2global[n_owned + i];
    owner_tmp[pos] = r;
  }
  for (ptrdiff_t pos = 0; pos < n_import; ++pos) {
    ghost_and_aura_to_owned[pos] = ghost_tmp[pos];
    local2global[n_owned + pos] = l2g_tmp[pos];
    owner[n_owned + pos] = owner_tmp[pos];
  }

  const ptrdiff_t nnz = n_total_elements > 0 ? e2n_ptr[n_total_elements] : 0;
  const idx_t owned_base = static_cast<idx_t>(n_owned);
  for (ptrdiff_t k = 0; k < nnz; ++k) {
    if (e2n_idx[k] >= owned_base) {
      e2n_idx[k] = owned_base + old_to_new[e2n_idx[k] - owned_base];
    }
  }

  SMESH_FREE(owner_tmp);
  SMESH_FREE(l2g_tmp);
  SMESH_FREE(ghost_tmp);
  SMESH_FREE(old_to_new);
  SMESH_FREE(ghost_displs);
  SMESH_FREE(aura_displs);
  SMESH_FREE(ghost_cursor);
  SMESH_FREE(aura_cursor);
  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t, typename large_idx_t>
static int mesh_create_parallel_e2n(
    const MPI_Comm comm, const int comm_size, const int comm_rank,
    ptrdiff_t *e2n_ptr, idx_t *e2n_idx, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements, const int spatial_dim, geom_t **points,
    const ptrdiff_t n_local2global, const ptrdiff_t n_global_nodes,
    ptrdiff_t *n_global_elements_out, ptrdiff_t *n_owned_elements_out,
    ptrdiff_t *n_shared_elements_out, ptrdiff_t *n_ghost_elements_out,
    large_idx_t **element_mapping_out, large_idx_t **aura_element_mapping_out,
    ptrdiff_t **e2n_ptr_out, idx_t **e2n_idx_out, int *spatial_dim_out,
    ptrdiff_t *n_global_nodes_out, ptrdiff_t *n_owned_nodes_out,
    ptrdiff_t *n_shared_nodes_out, ptrdiff_t *n_ghost_nodes_out,
    ptrdiff_t *n_aura_nodes_out, large_idx_t **node_mapping_out,
    geom_t ***points_out, int **node_owner_out, ptrdiff_t **node_offsets_out,
    idx_t **ghosts_out) {
  large_idx_t *n2eptr = nullptr;
  idx_t *n2e_idx = nullptr;
  if (create_n2e_from_e2n<idx_t, large_idx_t, idx_t>(
          comm, n_local_elements, n_global_elements, n_local2global,
          n_global_nodes, e2n_ptr, e2n_idx, &n2eptr, &n2e_idx) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }
  sort_n2e<large_idx_t, idx_t>(n_local2global, n2eptr, n2e_idx);

  ptrdiff_t local2global_size = 0;
  large_idx_t *local2global = nullptr;
  large_idx_t *local_n2e_ptr = nullptr;
  idx_t *local_n2e_idx = nullptr;
  redistribute_n2e(comm, comm_size, comm_rank, n_local2global, n_global_nodes,
                   n_global_elements, n2eptr, n2e_idx, &local2global_size,
                   &local2global, &local_n2e_ptr, &local_n2e_idx);
  SMESH_FREE(n2eptr);
  SMESH_FREE(n2e_idx);

  const ptrdiff_t nnz =
      n_local_elements > 0 ? e2n_ptr[n_local_elements] : 0;
  idx_t *local_e2n_idx =
      (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(nnz, 1) * sizeof(idx_t));
  localize_e2n<idx_t, large_idx_t, idx_t, large_idx_t>(
      comm_size, comm_rank, n_global_elements, n_local_elements, e2n_ptr,
      e2n_idx, local2global_size, local_n2e_ptr, local_n2e_idx, local2global,
      local_e2n_idx);
  SMESH_FREE(e2n_idx);
  e2n_idx = local_e2n_idx;

  ptrdiff_t n_owned = 0;
  ptrdiff_t n_shared = 0;
  ptrdiff_t n_ghosts = 0;
  rearrange_local_nodes_e2n<idx_t, large_idx_t, idx_t, large_idx_t>(
      comm_size, comm_rank, n_global_elements, n_local_elements,
      local2global_size, local_n2e_ptr, local_n2e_idx, local2global, e2n_ptr,
      e2n_idx, &n_owned, &n_shared, &n_ghosts);

  large_idx_t *element_mapping = (large_idx_t *)SMESH_ALLOC(
      (size_t)n_local_elements * sizeof(large_idx_t));
  ptrdiff_t n_owned_not_shared = 0;
  rearrange_local_elements_e2n<idx_t, large_idx_t, idx_t, large_idx_t>(
      comm_size, comm_rank, n_global_elements, n_local_elements,
      local2global_size, local_n2e_ptr, local_n2e_idx, e2n_ptr, e2n_idx, n_owned,
      &n_owned_not_shared, element_mapping);

  large_idx_t *aura_element_mapping = nullptr;
  ptrdiff_t *aura_e2n_ptr = nullptr;
  idx_t *aura_e2n_idx = nullptr;
  ptrdiff_t n_aura_elements = 0;
  expand_aura_elements_e2n<idx_t, large_idx_t, idx_t, large_idx_t>(
      comm, n_global_elements, n_local_elements, e2n_ptr, e2n_idx, local_n2e_ptr,
      local_n2e_idx, local2global, element_mapping, n_owned,
      &aura_element_mapping, &aura_e2n_ptr, &aura_e2n_idx, &n_aura_elements);
  SMESH_FREE(local_n2e_ptr);
  SMESH_FREE(local_n2e_idx);

  long long owned_nodes_start_ll = 0;
  long long n_owned_ll = (long long)n_owned;
  SMESH_MPI_CATCH(MPI_Exscan(&n_owned_ll, &owned_nodes_start_ll, 1, MPI_LONG_LONG,
                             MPI_SUM, comm));
  if (!comm_rank) {
    owned_nodes_start_ll = 0;
  }
  const ptrdiff_t owned_nodes_start =
      static_cast<ptrdiff_t>(owned_nodes_start_ll);

  idx_t *global2owned = (idx_t *)SMESH_CALLOC(
      (size_t)rank_split(n_global_nodes, comm_size, comm_rank), sizeof(idx_t));
  prepare_node_renumbering(comm, n_global_nodes, owned_nodes_start, n_owned,
                           local2global, global2owned);

  ptrdiff_t *owned_node_ranges =
      (ptrdiff_t *)SMESH_ALLOC((size_t)(comm_size + 1) * sizeof(ptrdiff_t));
  node_ownership_ranges(comm, n_owned, owned_node_ranges);

  large_idx_t *local2global_with_aura = nullptr;
  ptrdiff_t n_aura_nodes = 0;
  stich_aura_elements_e2n<idx_t, large_idx_t>(
      n_owned, n_shared, n_ghosts, local2global, n_local_elements, &e2n_ptr,
      &e2n_idx, n_aura_elements, aura_e2n_ptr, aura_e2n_idx,
      &local2global_with_aura, &n_aura_nodes);
  SMESH_FREE(aura_e2n_ptr);
  SMESH_FREE(aura_e2n_idx);
  SMESH_FREE(local2global);
  local2global = local2global_with_aura;
  local2global_size = n_owned + n_ghosts + n_aura_nodes;

  idx_t *ghost_and_aura_to_owned = (idx_t *)SMESH_ALLOC(
      (size_t)std::max<ptrdiff_t>(n_ghosts + n_aura_nodes, 1) * sizeof(idx_t));
  collect_ghost_and_aura_import_indices(
      comm, n_owned, n_ghosts, n_aura_nodes, n_global_nodes, local2global,
      global2owned, owned_node_ranges, ghost_and_aura_to_owned);

  node_ownership_ranges(comm, n_owned, owned_node_ranges);
  int *owner = (int *)SMESH_ALLOC(
      (size_t)(n_owned + n_ghosts + n_aura_nodes) * sizeof(int));
  determine_ownership(comm_size, comm_rank, n_owned, n_ghosts, n_aura_nodes,
                      ghost_and_aura_to_owned, owned_node_ranges, owner);

  const ptrdiff_t n_total_elements = n_local_elements + n_aura_elements;
  group_ghost_and_aura_e2n<idx_t, large_idx_t>(
      comm_size, n_owned, n_ghosts, n_aura_nodes, local2global,
      ghost_and_aura_to_owned, owner, n_total_elements, e2n_ptr, e2n_idx);

  const ptrdiff_t n_local_nodes = n_owned + n_ghosts + n_aura_nodes;
  geom_t **local_points = (geom_t **)SMESH_ALLOC((size_t)spatial_dim * sizeof(geom_t *));
  for (int d = 0; d < spatial_dim; ++d) {
    local_points[d] =
        (geom_t *)SMESH_ALLOC((size_t)n_local_nodes * sizeof(geom_t));
    gather_mapped_field(comm, n_local_nodes, n_global_nodes, local2global,
                        smesh::mpi_type<geom_t>(), points[d], local_points[d]);
    SMESH_FREE(points[d]);
  }
  SMESH_FREE(points);

  ptrdiff_t n_shared_elements = 0;
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    for (ptrdiff_t k = e2n_ptr[i]; k < e2n_ptr[i + 1]; ++k) {
      if (owner[e2n_idx[k]] != comm_rank) {
        n_shared_elements++;
        break;
      }
    }
  }

  *n_global_elements_out = n_global_elements;
  *n_owned_elements_out = n_local_elements;
  *n_shared_elements_out = n_shared_elements;
  *n_ghost_elements_out = n_aura_elements;
  *element_mapping_out = element_mapping;
  *aura_element_mapping_out = aura_element_mapping;
  *e2n_ptr_out = e2n_ptr;
  *e2n_idx_out = e2n_idx;
  *spatial_dim_out = spatial_dim;
  *n_global_nodes_out = n_global_nodes;
  *n_owned_nodes_out = n_owned;
  *n_shared_nodes_out = n_shared;
  *n_ghost_nodes_out = n_ghosts;
  *n_aura_nodes_out = n_aura_nodes;
  *node_mapping_out = local2global;
  *points_out = local_points;
  *node_owner_out = owner;
  *node_offsets_out = owned_node_ranges;
  *ghosts_out = ghost_and_aura_to_owned;
  SMESH_FREE(global2owned);
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif

