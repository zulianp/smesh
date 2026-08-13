#ifndef SMESH_DISTRIBUTED_REORDER_IMPL_HPP
#define SMESH_DISTRIBUTED_REORDER_IMPL_HPP

#include "smesh_distributed_reorder.hpp"

#include "smesh_decompose.hpp"
#include "smesh_distributed_aura.hpp"
#include "smesh_distributed_read.impl.hpp"
#include "smesh_sfc.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <type_traits>
#include <vector>

extern "C" {
#include "mpi-sort.h"
}

namespace smesh {

template <typename idx_t, typename geom_t>
static int gather_centroids_from_soa(
    MPI_Comm comm, const int nxe, const ptrdiff_t n_local,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    geom_t *const SMESH_RESTRICT cx, geom_t *const SMESH_RESTRICT cy,
    geom_t *const SMESH_RESTRICT cz) {
  geom_t *const tmp = (geom_t *)SMESH_ALLOC(
      (size_t)(n_local > 0 ? n_local : 1) * sizeof(geom_t));
  if (!tmp) {
    return SMESH_FAILURE;
  }

  for (ptrdiff_t e = 0; e < n_local; ++e) {
    cx[e] = 0;
    cy[e] = 0;
    cz[e] = 0;
  }

  geom_t *const out = n_local > 0 ? tmp : nullptr;
  for (int d = 0; d < nxe; ++d) {
    idx_t *const mapping = n_local > 0 ? elements[d] : nullptr;

    if (gather_mapped_field(comm, n_local, n_global_nodes, mapping,
                            mpi_type<geom_t>(), points[0],
                            out) != SMESH_SUCCESS) {
      SMESH_FREE(tmp);
      return SMESH_FAILURE;
    }
    for (ptrdiff_t e = 0; e < n_local; ++e) {
      cx[e] += tmp[e];
    }

    if (gather_mapped_field(comm, n_local, n_global_nodes, mapping,
                            mpi_type<geom_t>(), points[1],
                            out) != SMESH_SUCCESS) {
      SMESH_FREE(tmp);
      return SMESH_FAILURE;
    }
    for (ptrdiff_t e = 0; e < n_local; ++e) {
      cy[e] += tmp[e];
    }

    if (gather_mapped_field(comm, n_local, n_global_nodes, mapping,
                            mpi_type<geom_t>(), points[2],
                            out) != SMESH_SUCCESS) {
      SMESH_FREE(tmp);
      return SMESH_FAILURE;
    }
    for (ptrdiff_t e = 0; e < n_local; ++e) {
      cz[e] += tmp[e];
    }
  }

  if (n_local > 0) {
    const geom_t inv = geom_t(1) / geom_t(nxe);
    for (ptrdiff_t e = 0; e < n_local; ++e) {
      cx[e] *= inv;
      cy[e] *= inv;
      cz[e] *= inv;
    }
  }

  SMESH_FREE(tmp);
  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t, typename Ordering>
int distributed_reorder_elements(
    MPI_Comm comm, const int nnodesxelem, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    large_idx_t *const SMESH_RESTRICT sorted_ids, Ordering ordering) {
  if (n_local_elements == 0) {
    return SMESH_SUCCESS;
  }
  if (!sorted_ids) {
    return SMESH_FAILURE;
  }

  int rank = 0;
  int size = 1;
  SMESH_MPI_CATCH(MPI_Comm_rank(comm, &rank));
  SMESH_MPI_CATCH(MPI_Comm_size(comm, &size));

  const ptrdiff_t n_sorted_elements = rank_split(n_global_elements, size, rank);
  if (n_sorted_elements != n_local_elements) {
    SMESH_ERROR("In-place distributed element sorting requires rank-split "
                "element ownership");
    return SMESH_FAILURE;
  }

  std::vector<geom_t> cx((size_t)n_local_elements);
  std::vector<geom_t> cy((size_t)n_local_elements);
  std::vector<geom_t> cz((size_t)n_local_elements);
  if (gather_centroids_from_soa(comm, nnodesxelem, n_local_elements, elements,
                                n_global_nodes, points, cx.data(), cy.data(),
                                cz.data()) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max()};
  geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest()};
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    local_min[0] = std::min(local_min[0], cx[(size_t)e]);
    local_min[1] = std::min(local_min[1], cy[(size_t)e]);
    local_min[2] = std::min(local_min[2], cz[(size_t)e]);
    local_max[0] = std::max(local_max[0], cx[(size_t)e]);
    local_max[1] = std::max(local_max[1], cy[(size_t)e]);
    local_max[2] = std::max(local_max[2], cz[(size_t)e]);
  }

  geom_t global_min[3];
  geom_t global_max[3];
  SMESH_MPI_CATCH(MPI_Allreduce(local_min, global_min, 3, mpi_type<geom_t>(),
                                MPI_MIN, comm));
  SMESH_MPI_CATCH(MPI_Allreduce(local_max, global_max, 3, mpi_type<geom_t>(),
                                MPI_MAX, comm));

  std::vector<u32> send_keys((size_t)n_local_elements);
  if (ordering(n_local_elements, cx.data(), cy.data(), cz.data(), global_min[0],
               global_max[0], global_min[1], global_max[1], global_min[2],
               global_max[2], send_keys.data()) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  const ptrdiff_t element_start = rank_start(n_global_elements, size, rank);
  std::vector<large_idx_t> send_ids((size_t)n_local_elements);
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    send_ids[(size_t)e] = static_cast<large_idx_t>(element_start + e);
  }

  std::vector<u32> sorted_keys((size_t)n_sorted_elements);
  if (MPI_Sort_bykey(send_keys.data(), send_ids.data(),
                     static_cast<int>(n_local_elements), mpi_type<u32>(),
                     mpi_type<large_idx_t>(), sorted_keys.data(), sorted_ids,
                     static_cast<int>(n_sorted_elements),
                     comm) != MPI_SUCCESS) {
    return SMESH_FAILURE;
  }

#ifndef NDEBUG
  for (ptrdiff_t e = 1; e < n_sorted_elements; ++e) {
    SMESH_ASSERT(sorted_keys[(size_t)(e - 1)] <= sorted_keys[(size_t)e]);
  }

  if (rank + 1 < size && n_sorted_elements > 0) {
    const u32 local_last = sorted_keys.back();
    u32 next_first = 0;
    SMESH_MPI_CATCH(MPI_Sendrecv(&local_last, 1, mpi_type<u32>(), rank + 1, 0,
                                 &next_first, 1, mpi_type<u32>(), rank + 1, 1,
                                 comm, MPI_STATUS_IGNORE));
    SMESH_ASSERT(local_last <= next_first);
  }

  if (rank > 0 && n_sorted_elements > 0) {
    const u32 local_first = sorted_keys.front();
    u32 prev_last = 0;
    SMESH_MPI_CATCH(MPI_Sendrecv(&local_first, 1, mpi_type<u32>(), rank - 1, 1,
                                 &prev_last, 1, mpi_type<u32>(), rank - 1, 0,
                                 comm, MPI_STATUS_IGNORE));
    SMESH_ASSERT(prev_last <= local_first);
  }
#endif

  std::vector<idx_t> sorted_elements((size_t)n_sorted_elements);
  for (int d = 0; d < nnodesxelem; ++d) {
    if (gather_mapped_field(comm, n_sorted_elements, n_global_elements,
                            sorted_ids, mpi_type<idx_t>(), elements[d],
                            sorted_elements.data()) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    std::memcpy(elements[d], sorted_elements.data(),
                (size_t)n_sorted_elements * sizeof(idx_t));
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t, typename Ordering>
int distributed_assign_elements_sfc_multiblock(
    MPI_Comm comm, const block_idx_t n_blocks,
    const int *const SMESH_RESTRICT nnodesxelem,
    const ptrdiff_t *const SMESH_RESTRICT n_local_elements,
    const ptrdiff_t *const SMESH_RESTRICT n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
        *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    ptrdiff_t *const SMESH_RESTRICT n_assigned_out,
    large_idx_t **const SMESH_RESTRICT sorted_concat_ids_out,
    idx_t ***const SMESH_RESTRICT assigned_elements_out, const bool use_sfc,
    Ordering ordering) {
  if (!n_assigned_out || !sorted_concat_ids_out || !assigned_elements_out ||
      n_blocks <= 0) {
    return SMESH_FAILURE;
  }

  int rank = 0;
  int size = 1;
  SMESH_MPI_CATCH(MPI_Comm_rank(comm, &rank));
  SMESH_MPI_CATCH(MPI_Comm_size(comm, &size));

  const int nxe0 = nnodesxelem[0];
  ptrdiff_t n_global_total = 0;
  ptrdiff_t n_send = 0;
  std::vector<ptrdiff_t> concat_offset((size_t)n_blocks + 1, 0);
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    if (nnodesxelem[b] != nxe0) {
      SMESH_ERROR("distributed_assign_elements_sfc_multiblock: heterogeneous "
                  "nnodesxelem is not supported in A1 (block %d has %d, "
                  "expected %d)\n",
                  (int)b, nnodesxelem[b], nxe0);
      return SMESH_FAILURE;
    }
    concat_offset[(size_t)b] = n_global_total;
    n_global_total += n_global_elements[b];
    n_send += n_local_elements[b];
  }
  concat_offset[(size_t)n_blocks] = n_global_total;

  if (n_global_total < size) {
    SMESH_ERROR("distributed_assign_elements_sfc_multiblock: need at least "
                "one element per rank (n_global=%ld, size=%d)\n",
                (long)n_global_total, size);
    return SMESH_FAILURE;
  }

  const ptrdiff_t n_recv = rank_split(n_global_total, size, rank);
  *n_assigned_out = n_recv;

  std::vector<u32> send_keys((size_t)std::max<ptrdiff_t>(n_send, 1));
  std::vector<large_idx_t> send_ids((size_t)std::max<ptrdiff_t>(n_send, 1));

  ptrdiff_t cursor = 0;
  if (use_sfc) {
    std::vector<geom_t> cx((size_t)std::max<ptrdiff_t>(n_send, 1));
    std::vector<geom_t> cy((size_t)std::max<ptrdiff_t>(n_send, 1));
    std::vector<geom_t> cz((size_t)std::max<ptrdiff_t>(n_send, 1));
    geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                           std::numeric_limits<geom_t>::max(),
                           std::numeric_limits<geom_t>::max()};
    geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                           std::numeric_limits<geom_t>::lowest(),
                           std::numeric_limits<geom_t>::lowest()};

    for (block_idx_t b = 0; b < n_blocks; ++b) {
      const ptrdiff_t n_local_b = n_local_elements[b];
      const ptrdiff_t n_global_b = n_global_elements[b];
      const int nxe = nnodesxelem[b];

      if (gather_centroids_from_soa(
              comm, nxe, n_local_b, elements[b], n_global_nodes, points,
              n_local_b > 0 ? cx.data() + cursor : nullptr,
              n_local_b > 0 ? cy.data() + cursor : nullptr,
              n_local_b > 0 ? cz.data() + cursor : nullptr) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }

      const ptrdiff_t element_start = rank_start(n_global_b, size, rank);
      for (ptrdiff_t e = 0; e < n_local_b; ++e) {
        const geom_t x = cx[(size_t)(cursor + e)];
        const geom_t y = cy[(size_t)(cursor + e)];
        const geom_t z = cz[(size_t)(cursor + e)];
        local_min[0] = std::min(local_min[0], x);
        local_min[1] = std::min(local_min[1], y);
        local_min[2] = std::min(local_min[2], z);
        local_max[0] = std::max(local_max[0], x);
        local_max[1] = std::max(local_max[1], y);
        local_max[2] = std::max(local_max[2], z);
        send_ids[(size_t)(cursor + e)] = static_cast<large_idx_t>(
            concat_offset[(size_t)b] + element_start + e);
      }
      cursor += n_local_b;
    }

    geom_t global_min[3];
    geom_t global_max[3];
    SMESH_MPI_CATCH(MPI_Allreduce(local_min, global_min, 3, mpi_type<geom_t>(),
                                  MPI_MIN, comm));
    SMESH_MPI_CATCH(MPI_Allreduce(local_max, global_max, 3, mpi_type<geom_t>(),
                                  MPI_MAX, comm));
    if (n_send > 0) {
      if (ordering(n_send, cx.data(), cy.data(), cz.data(), global_min[0],
                   global_max[0], global_min[1], global_max[1], global_min[2],
                   global_max[2], send_keys.data()) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
    }
  } else {
    for (block_idx_t b = 0; b < n_blocks; ++b) {
      const ptrdiff_t n_local_b = n_local_elements[b];
      if (n_local_b == 0) {
        continue;
      }
      const ptrdiff_t element_start =
          rank_start(n_global_elements[b], size, rank);
      for (ptrdiff_t e = 0; e < n_local_b; ++e) {
        send_ids[(size_t)cursor] = static_cast<large_idx_t>(
            concat_offset[(size_t)b] + element_start + e);
        ++cursor;
      }
    }
  }

  large_idx_t *sorted_ids =
      (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_recv, 1) *
                                 sizeof(large_idx_t));
  if (use_sfc) {
    std::vector<u32> sorted_keys((size_t)std::max<ptrdiff_t>(n_recv, 1));
    if (MPI_Sort_bykey(n_send > 0 ? send_keys.data() : nullptr,
                       n_send > 0 ? send_ids.data() : nullptr,
                       static_cast<int>(n_send), mpi_type<u32>(),
                       mpi_type<large_idx_t>(),
                       n_recv > 0 ? sorted_keys.data() : nullptr,
                       n_recv > 0 ? sorted_ids : nullptr,
                       static_cast<int>(n_recv), comm) != MPI_SUCCESS) {
      SMESH_FREE(sorted_ids);
      return SMESH_FAILURE;
    }
  } else {
    // File/concat order: sort by concat_id → rank_split(sum) ownership.
    std::vector<large_idx_t> send_keys_id(
        (size_t)std::max<ptrdiff_t>(n_send, 1));
    for (ptrdiff_t i = 0; i < n_send; ++i) {
      send_keys_id[(size_t)i] = send_ids[(size_t)i];
    }
    large_idx_t *sorted_keys =
        (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_recv, 1) *
                                   sizeof(large_idx_t));
    if (MPI_Sort_bykey(n_send > 0 ? send_keys_id.data() : nullptr,
                       n_send > 0 ? send_ids.data() : nullptr,
                       static_cast<int>(n_send), mpi_type<large_idx_t>(),
                       mpi_type<large_idx_t>(),
                       n_recv > 0 ? sorted_keys : nullptr,
                       n_recv > 0 ? sorted_ids : nullptr,
                       static_cast<int>(n_recv), comm) != MPI_SUCCESS) {
      SMESH_FREE(sorted_keys);
      SMESH_FREE(sorted_ids);
      return SMESH_FAILURE;
    }
    SMESH_FREE(sorted_keys);
  }

  idx_t **assigned =
      (idx_t **)SMESH_ALLOC((size_t)nxe0 * sizeof(idx_t *));
  for (int d = 0; d < nxe0; ++d) {
    assigned[d] =
        (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_recv, 1) *
                             sizeof(idx_t));
  }

  // Gather connectivity per original block into assigned SoA.
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    std::vector<large_idx_t> block_ids;
    std::vector<ptrdiff_t> dest_pos;
    block_ids.reserve((size_t)n_recv);
    dest_pos.reserve((size_t)n_recv);
    const ptrdiff_t off = concat_offset[(size_t)b];
    const ptrdiff_t n_global_b = n_global_elements[b];
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
      const large_idx_t cid = sorted_ids[i];
      if (cid < static_cast<large_idx_t>(off) ||
          cid >= static_cast<large_idx_t>(off + n_global_b)) {
        continue;
      }
      block_ids.push_back(cid - static_cast<large_idx_t>(off));
      dest_pos.push_back(i);
    }

    const ptrdiff_t n_b = static_cast<ptrdiff_t>(block_ids.size());
    if (n_b == 0) {
      // Still need a collective gather_mapped_field participation? gather
      // uses Alltoall internally so all ranks must call it the same number
      // of times. Call with n_local=0.
      std::vector<idx_t> empty_tmp;
      for (int d = 0; d < nxe0; ++d) {
        if (gather_mapped_field(comm, 0, n_global_b, (large_idx_t *)nullptr,
                                mpi_type<idx_t>(), elements[b][d],
                                (idx_t *)nullptr) != SMESH_SUCCESS) {
          for (int dd = 0; dd < nxe0; ++dd) {
            SMESH_FREE(assigned[dd]);
          }
          SMESH_FREE(assigned);
          SMESH_FREE(sorted_ids);
          return SMESH_FAILURE;
        }
      }
      continue;
    }

    std::vector<idx_t> tmp((size_t)n_b);
    for (int d = 0; d < nxe0; ++d) {
      if (gather_mapped_field(comm, n_b, n_global_b, block_ids.data(),
                              mpi_type<idx_t>(), elements[b][d],
                              tmp.data()) != SMESH_SUCCESS) {
        for (int dd = 0; dd < nxe0; ++dd) {
          SMESH_FREE(assigned[dd]);
        }
        SMESH_FREE(assigned);
        SMESH_FREE(sorted_ids);
        return SMESH_FAILURE;
      }
      for (ptrdiff_t j = 0; j < n_b; ++j) {
        assigned[d][dest_pos[(size_t)j]] = tmp[(size_t)j];
      }
    }
  }

  *sorted_concat_ids_out = sorted_ids;
  *assigned_elements_out = assigned;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t, typename global_idx_t,
          typename Ordering>
int mesh_from_folder_reordered_basic(
    const MPI_Comm comm, const Path &folder,
    // Elements
    int *nnodesxelem_out, ptrdiff_t *n_global_elements_out,
    ptrdiff_t *n_owned_elements_out, ptrdiff_t *n_shared_elements_out,
    ptrdiff_t *n_ghost_elements_out, global_idx_t **element_mapping_out,
    global_idx_t **aura_element_mapping_out, idx_t ***elements_out,
    // Nodes
    int *spatial_dim_out, ptrdiff_t *n_global_nodes_out,
    ptrdiff_t *n_owned_nodes_out, ptrdiff_t *n_shared_nodes_out,
    ptrdiff_t *n_ghost_nodes_out, ptrdiff_t *n_aura_nodes_out,
    global_idx_t **node_mapping_out, geom_t ***points_out,
    // Distributed connectivities
    int **node_owner_out, ptrdiff_t **node_offsets_out, idx_t **ghosts_out,
    Ordering ordering) {
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  int nnodesxelem;
  idx_t **elems = nullptr;
  ptrdiff_t n_local_elements;
  ptrdiff_t n_global_elements;
  if (mesh_block_from_folder<idx_t>(comm, folder, &nnodesxelem, &elems,
                                    &n_local_elements,
                                    &n_global_elements) != SMESH_SUCCESS) {
    SMESH_ERROR("Failed to read mesh blocks\n");
    return SMESH_FAILURE;
  }

  int spatial_dim;
  geom_t **points = nullptr;
  ptrdiff_t n_local2global;
  ptrdiff_t n_global_nodes;
  if (mesh_coordinates_from_folder(comm, folder, &spatial_dim, &points,
                                   &n_local2global,
                                   &n_global_nodes) != SMESH_SUCCESS) {
    for (int d = 0; d < nnodesxelem; ++d) {
      SMESH_FREE(elems[d]);
    }
    SMESH_FREE(elems);

    SMESH_ERROR("Failed to read coordinates\n");
    return SMESH_FAILURE;
  }

  smesh::large_idx_t *sorted_ids = (smesh::large_idx_t *)SMESH_ALLOC(
      n_local_elements * sizeof(smesh::large_idx_t));
  if (distributed_reorder_elements<idx_t, geom_t, Ordering>(
          comm, nnodesxelem, n_local_elements, n_global_elements, elems,
          n_global_nodes, points, sorted_ids, ordering) != SMESH_SUCCESS) {
    for (int d = 0; d < nnodesxelem; ++d) {
      SMESH_FREE(elems[d]);
    }
    SMESH_FREE(elems);
    for (int d = 0; d < spatial_dim; ++d) {
      SMESH_FREE(points[d]);
    }
    SMESH_FREE(points);
    SMESH_FREE(sorted_ids);
    return SMESH_FAILURE;
  }

  const int ret = mesh_create_parallel<idx_t, geom_t, global_idx_t>(
      comm, comm_size, comm_rank, nnodesxelem, elems, n_local_elements,
      n_global_elements, spatial_dim, points, n_local2global, n_global_nodes,
      nullptr, nnodesxelem_out, n_global_elements_out, n_owned_elements_out,
      n_shared_elements_out, n_ghost_elements_out, element_mapping_out,
      aura_element_mapping_out, elements_out, spatial_dim_out,
      n_global_nodes_out, n_owned_nodes_out, n_shared_nodes_out,
      n_ghost_nodes_out, n_aura_nodes_out, node_mapping_out, points_out,
      node_owner_out, node_offsets_out, ghosts_out);
  if (ret == SMESH_SUCCESS) {
    global_idx_t *const element_mapping = *element_mapping_out;
    for (ptrdiff_t i = 0; i < *n_owned_elements_out; ++i) {
      const ptrdiff_t local_id =
          static_cast<ptrdiff_t>(element_mapping[i]) -
          rank_start(n_global_elements, comm_size, comm_rank);
      element_mapping[i] = static_cast<global_idx_t>(sorted_ids[local_id]);
    }

    const ptrdiff_t n_aura_elements = *n_ghost_elements_out;
    global_idx_t *const aura_element_mapping = *aura_element_mapping_out;
    smesh::large_idx_t *aura_element_mapping_large =
        n_aura_elements > 0 ? (smesh::large_idx_t *)SMESH_ALLOC(
                                  n_aura_elements * sizeof(smesh::large_idx_t))
                            : nullptr;
    gather_mapped_field(comm, n_aura_elements, n_global_elements,
                        aura_element_mapping,
                        smesh::mpi_type<smesh::large_idx_t>(), sorted_ids,
                        aura_element_mapping_large);
    for (ptrdiff_t i = 0; i < n_aura_elements; ++i) {
      aura_element_mapping[i] =
          static_cast<global_idx_t>(aura_element_mapping_large[i]);
    }
    if (aura_element_mapping_large) {
      SMESH_FREE(aura_element_mapping_large);
    }
  }
  SMESH_FREE(sorted_ids);
  return ret;
}

template <typename idx_t, typename geom_t, typename global_idx_t,
          typename Ordering>
int mesh_from_folder_reordered(
    const MPI_Comm comm, const Path &folder,
    // Elements
    int *nnodesxelem_out, ptrdiff_t *n_global_elements_out,
    ptrdiff_t *n_owned_elements_out, ptrdiff_t *n_shared_elements_out,
    ptrdiff_t *n_ghost_elements_out, global_idx_t **element_mapping_out,
    global_idx_t **aura_element_mapping_out, idx_t ***elements_out,
    // Nodes
    int *spatial_dim_out, ptrdiff_t *n_global_nodes_out,
    ptrdiff_t *n_owned_nodes_out, ptrdiff_t *n_shared_nodes_out,
    ptrdiff_t *n_ghost_nodes_out, ptrdiff_t *n_aura_nodes_out,
    global_idx_t **node_mapping_out, geom_t ***points_out,
    // Distributed connectivities
    int **node_owner_out, ptrdiff_t **node_offsets_out, idx_t **ghosts_out,
    Ordering ordering) {
  if constexpr (std::is_same_v<global_idx_t, idx_t>) {
    return mesh_from_folder_reordered_basic<idx_t, geom_t, global_idx_t,
                                            Ordering>(
        comm, folder, nnodesxelem_out, n_global_elements_out,
        n_owned_elements_out, n_shared_elements_out, n_ghost_elements_out,
        element_mapping_out, aura_element_mapping_out, elements_out,
        spatial_dim_out, n_global_nodes_out, n_owned_nodes_out,
        n_shared_nodes_out, n_ghost_nodes_out, n_aura_nodes_out,
        node_mapping_out, points_out, node_owner_out, node_offsets_out,
        ghosts_out, ordering);
  } else {
    std::vector<Path> i_files =
        detect_files(folder / "i0.*", {"raw", "int16", "int32", "int64"});
    std::vector<Path> x_files =
        detect_files(folder / "x*.*", {"raw", "float16", "float32", "float64"});

    if (i_files.empty() || x_files.empty()) {
      SMESH_ERROR("No mesh files found in folder %s\n", folder.c_str());
      return SMESH_FAILURE;
    }

    const size_t num_e_idx = file_size(i_files[0]) /
                             num_bytes(to_integer_type(i_files[0].extension()));
    const size_t num_nodes =
        file_size(x_files[0]) / num_bytes(to_real_type(x_files[0].extension()));
    const size_t max_idx = std::max(num_e_idx, num_nodes);

    if (max_idx < std::numeric_limits<idx_t>::max()) {
      return mesh_from_folder_reordered_basic<idx_t, geom_t, global_idx_t,
                                              Ordering>(
          comm, folder, nnodesxelem_out, n_global_elements_out,
          n_owned_elements_out, n_shared_elements_out, n_ghost_elements_out,
          element_mapping_out, aura_element_mapping_out, elements_out,
          spatial_dim_out, n_global_nodes_out, n_owned_nodes_out,
          n_shared_nodes_out, n_ghost_nodes_out, n_aura_nodes_out,
          node_mapping_out, points_out, node_owner_out, node_offsets_out,
          ghosts_out, ordering);
    }

    int rank;
    MPI_Comm_rank(comm, &rank);
    if (!rank) {
      printf("Found: %ld Nodes, Using large index type\n", (long)num_nodes);
      fflush(stdout);
    }

    static_assert(sizeof(global_idx_t) > sizeof(idx_t),
                  "global_idx_t must be larger than idx_t!");

    global_idx_t **elements;
    global_idx_t *ghosts;

    if (mesh_from_folder_reordered_basic<global_idx_t, geom_t, global_idx_t,
                                         Ordering>(
            comm, folder, nnodesxelem_out, n_global_elements_out,
            n_owned_elements_out, n_shared_elements_out, n_ghost_elements_out,
            element_mapping_out, aura_element_mapping_out, &elements,
            spatial_dim_out, n_global_nodes_out, n_owned_nodes_out,
            n_shared_nodes_out, n_ghost_nodes_out, n_aura_nodes_out,
            node_mapping_out, points_out, node_owner_out, node_offsets_out,
            &ghosts, ordering) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    const ptrdiff_t nnodesxelem = *nnodesxelem_out;
    const ptrdiff_t n_local_elements = *n_owned_elements_out;
    const ptrdiff_t n_local_nodes =
        *n_owned_nodes_out + *n_ghost_nodes_out + *n_aura_nodes_out;

    if (static_cast<long double>(n_local_nodes) >
        static_cast<long double>(std::numeric_limits<idx_t>::max())) {
      SMESH_ERROR(
          "Distributed read requires %ld local node indices on rank %d, "
          "but idx_t can represent at most %lld. Rebuild with a wider "
          "SMESH_IDX_TYPE.\n",
          (long)n_local_nodes, rank,
          (long long)std::numeric_limits<idx_t>::max());
      return SMESH_FAILURE;
    }

    idx_t **small_elements =
        (idx_t **)SMESH_ALLOC(nnodesxelem * sizeof(idx_t *));

    for (int d = 0; d < nnodesxelem; ++d) {
      small_elements[d] =
          (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));

      for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
        small_elements[d][i] = (idx_t)elements[d][i];
      }

      SMESH_FREE(elements[d]);
    }
    SMESH_FREE(elements);

    *elements_out = small_elements;

    const ptrdiff_t n_import_nodes = *n_ghost_nodes_out + *n_aura_nodes_out;
    *ghosts_out = (idx_t *)SMESH_ALLOC(n_import_nodes * sizeof(idx_t));
    for (ptrdiff_t i = 0; i < n_import_nodes; ++i) {
      (*ghosts_out)[i] = (idx_t)ghosts[i];
    }

    SMESH_FREE(ghosts);

    return SMESH_SUCCESS;
  }
}

template <typename idx_t, typename geom_t, typename large_idx_t,
          typename Ordering>
int mesh_from_folder_multiblock(
    const MPI_Comm comm, const Path &folder,
    const std::vector<std::string> &block_names,
    const std::vector<enum ElemType> &element_types,
    ptrdiff_t *n_global_elements_out, ptrdiff_t *n_owned_elements_out,
    ptrdiff_t *n_shared_elements_out, ptrdiff_t *n_ghost_elements_out,
    large_idx_t **element_mapping_out, large_idx_t **aura_element_mapping_out,
    int *nnodesxelem_out, ptrdiff_t **n_local_elements_per_block_out,
    idx_t ****elements_per_block_out, int *spatial_dim_out,
    ptrdiff_t *n_global_nodes_out, ptrdiff_t *n_owned_nodes_out,
    ptrdiff_t *n_shared_nodes_out, ptrdiff_t *n_ghost_nodes_out,
    ptrdiff_t *n_aura_nodes_out, large_idx_t **node_mapping_out,
    geom_t ***points_out, int **node_owner_out, ptrdiff_t **node_offsets_out,
    idx_t **ghosts_out, const bool use_sfc, Ordering ordering) {
  (void)element_types;
  const block_idx_t n_blocks = static_cast<block_idx_t>(block_names.size());
  if (n_blocks == 0 || block_names.size() != element_types.size() ||
      !n_local_elements_per_block_out || !elements_per_block_out ||
      !nnodesxelem_out) {
    return SMESH_FAILURE;
  }

  int comm_rank = 0;
  int comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  int spatial_dim = 0;
  geom_t **points = nullptr;
  ptrdiff_t n_local_nodes_file = 0;
  ptrdiff_t n_global_nodes = 0;
  if (mesh_coordinates_from_folder(comm, folder, &spatial_dim, &points,
                                   &n_local_nodes_file,
                                   &n_global_nodes) != SMESH_SUCCESS) {
    SMESH_ERROR("mesh_from_folder_multiblock: failed to read coordinates\n");
    return SMESH_FAILURE;
  }

  std::vector<int> nxe((size_t)n_blocks);
  std::vector<ptrdiff_t> n_local_e((size_t)n_blocks);
  std::vector<ptrdiff_t> n_global_e((size_t)n_blocks);
  std::vector<idx_t **> elems((size_t)n_blocks, nullptr);

  for (block_idx_t b = 0; b < n_blocks; ++b) {
    Path block_folder = Path(folder) / "blocks" / block_names[(size_t)b];
    if (mesh_block_from_folder<idx_t>(comm, block_folder, &nxe[(size_t)b],
                                      &elems[(size_t)b], &n_local_e[(size_t)b],
                                      &n_global_e[(size_t)b]) != SMESH_SUCCESS) {
      for (block_idx_t bb = 0; bb < b; ++bb) {
        for (int d = 0; d < nxe[(size_t)bb]; ++d) {
          SMESH_FREE(elems[(size_t)bb][d]);
        }
        SMESH_FREE(elems[(size_t)bb]);
      }
      for (int d = 0; d < spatial_dim; ++d) {
        SMESH_FREE(points[d]);
      }
      SMESH_FREE(points);
      SMESH_ERROR("mesh_from_folder_multiblock: failed to read block %s\n",
                  block_names[(size_t)b].c_str());
      return SMESH_FAILURE;
    }
  }

  ptrdiff_t n_assigned = 0;
  smesh::large_idx_t *sorted_concat_ids = nullptr;
  idx_t **assigned_elements = nullptr;
  if (distributed_assign_elements_sfc_multiblock<idx_t, geom_t, Ordering>(
          comm, n_blocks, nxe.data(), n_local_e.data(), n_global_e.data(),
          elems.data(), n_global_nodes, points, &n_assigned, &sorted_concat_ids,
          &assigned_elements, use_sfc, ordering) != SMESH_SUCCESS) {
    for (block_idx_t b = 0; b < n_blocks; ++b) {
      for (int d = 0; d < nxe[(size_t)b]; ++d) {
        SMESH_FREE(elems[(size_t)b][d]);
      }
      SMESH_FREE(elems[(size_t)b]);
    }
    for (int d = 0; d < spatial_dim; ++d) {
      SMESH_FREE(points[d]);
    }
    SMESH_FREE(points);
    return SMESH_FAILURE;
  }

  for (block_idx_t b = 0; b < n_blocks; ++b) {
    for (int d = 0; d < nxe[(size_t)b]; ++d) {
      SMESH_FREE(elems[(size_t)b][d]);
    }
    SMESH_FREE(elems[(size_t)b]);
  }

  ptrdiff_t n_global_total = 0;
  std::vector<ptrdiff_t> concat_offset((size_t)n_blocks + 1, 0);
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    concat_offset[(size_t)b] = n_global_total;
    n_global_total += n_global_e[(size_t)b];
  }
  concat_offset[(size_t)n_blocks] = n_global_total;

  const int nxe0 = nxe[0];
  idx_t **concat_elements = nullptr;
  // Pass nullptr mapping; remap to concat_ids after rearrange (same pattern as
  // mesh_from_folder_reordered_basic). Avoids large_idx_t vs smesh::large_idx_t
  // template mismatches.
  const int ret = mesh_create_parallel<idx_t, geom_t, large_idx_t>(
      comm, comm_size, comm_rank, nxe0, assigned_elements, n_assigned,
      n_global_total, spatial_dim, points, n_local_nodes_file, n_global_nodes,
      nullptr, nnodesxelem_out, n_global_elements_out, n_owned_elements_out,
      n_shared_elements_out, n_ghost_elements_out, element_mapping_out,
      aura_element_mapping_out, &concat_elements, spatial_dim_out,
      n_global_nodes_out, n_owned_nodes_out, n_shared_nodes_out,
      n_ghost_nodes_out, n_aura_nodes_out, node_mapping_out, points_out,
      node_owner_out, node_offsets_out, ghosts_out);

  if (ret != SMESH_SUCCESS) {
    SMESH_FREE(sorted_concat_ids);
    return SMESH_FAILURE;
  }

  {
    large_idx_t *const element_mapping = *element_mapping_out;
    const ptrdiff_t elements_start =
        rank_start(n_global_total, comm_size, comm_rank);
    for (ptrdiff_t i = 0; i < *n_owned_elements_out; ++i) {
      const ptrdiff_t local_id =
          static_cast<ptrdiff_t>(element_mapping[i]) - elements_start;
      element_mapping[i] =
          static_cast<large_idx_t>(sorted_concat_ids[local_id]);
    }

    const ptrdiff_t n_aura_elements = *n_ghost_elements_out;
    large_idx_t *const aura_element_mapping = *aura_element_mapping_out;
    smesh::large_idx_t *aura_element_mapping_large =
        n_aura_elements > 0
            ? (smesh::large_idx_t *)SMESH_ALLOC(
                  (size_t)n_aura_elements * sizeof(smesh::large_idx_t))
            : nullptr;
    gather_mapped_field(comm, n_aura_elements, n_global_total,
                        aura_element_mapping,
                        smesh::mpi_type<smesh::large_idx_t>(), sorted_concat_ids,
                        aura_element_mapping_large);
    for (ptrdiff_t i = 0; i < n_aura_elements; ++i) {
      aura_element_mapping[i] =
          static_cast<large_idx_t>(aura_element_mapping_large[i]);
    }
    if (aura_element_mapping_large) {
      SMESH_FREE(aura_element_mapping_large);
    }
  }
  SMESH_FREE(sorted_concat_ids);
  sorted_concat_ids = nullptr;

  const ptrdiff_t n_owned = *n_owned_elements_out;
  const ptrdiff_t n_aura = *n_ghost_elements_out;
  const large_idx_t *const element_mapping = *element_mapping_out;
  const large_idx_t *const aura_element_mapping = *aura_element_mapping_out;

  auto block_of_concat = [&](const large_idx_t cid) -> block_idx_t {
    for (block_idx_t b = 0; b < n_blocks; ++b) {
      if (cid >= static_cast<large_idx_t>(concat_offset[(size_t)b]) &&
          cid < static_cast<large_idx_t>(concat_offset[(size_t)(b + 1)])) {
        return b;
      }
    }
    return n_blocks;
  };

  ptrdiff_t *n_local_per_block =
      (ptrdiff_t *)SMESH_CALLOC((size_t)n_blocks, sizeof(ptrdiff_t));
  for (ptrdiff_t i = 0; i < n_owned; ++i) {
    const block_idx_t b = block_of_concat(element_mapping[i]);
    if (b >= n_blocks) {
      SMESH_ERROR("mesh_from_folder_multiblock: owned concat_id %lld out of "
                  "range\n",
                  (long long)element_mapping[i]);
      for (int d = 0; d < nxe0; ++d) {
        SMESH_FREE(concat_elements[d]);
      }
      SMESH_FREE(concat_elements);
      SMESH_FREE(n_local_per_block);
      return SMESH_FAILURE;
    }
    n_local_per_block[b]++;
  }
  for (ptrdiff_t i = 0; i < n_aura; ++i) {
    const block_idx_t b = block_of_concat(aura_element_mapping[i]);
    if (b >= n_blocks) {
      SMESH_ERROR("mesh_from_folder_multiblock: aura concat_id %lld out of "
                  "range\n",
                  (long long)aura_element_mapping[i]);
      for (int d = 0; d < nxe0; ++d) {
        SMESH_FREE(concat_elements[d]);
      }
      SMESH_FREE(concat_elements);
      SMESH_FREE(n_local_per_block);
      return SMESH_FAILURE;
    }
    n_local_per_block[b]++;
  }

  idx_t ***elements_per_block =
      (idx_t ***)SMESH_ALLOC((size_t)n_blocks * sizeof(idx_t **));
  std::vector<ptrdiff_t> write_cursor((size_t)n_blocks, 0);
  for (block_idx_t b = 0; b < n_blocks; ++b) {
    elements_per_block[b] = (idx_t **)SMESH_ALLOC((size_t)nxe0 * sizeof(idx_t *));
    const ptrdiff_t n_b =
        std::max<ptrdiff_t>(n_local_per_block[b], 1);
    for (int d = 0; d < nxe0; ++d) {
      elements_per_block[b][d] =
          (idx_t *)SMESH_ALLOC((size_t)n_b * sizeof(idx_t));
    }
  }

  auto append_elem = [&](const block_idx_t b, const ptrdiff_t src) {
    const ptrdiff_t dst = write_cursor[(size_t)b]++;
    for (int d = 0; d < nxe0; ++d) {
      elements_per_block[b][d][dst] = concat_elements[d][src];
    }
  };

  for (ptrdiff_t i = 0; i < n_owned; ++i) {
    append_elem(block_of_concat(element_mapping[i]), i);
  }
  for (ptrdiff_t i = 0; i < n_aura; ++i) {
    append_elem(block_of_concat(aura_element_mapping[i]), n_owned + i);
  }

  for (int d = 0; d < nxe0; ++d) {
    SMESH_FREE(concat_elements[d]);
  }
  SMESH_FREE(concat_elements);

  *n_local_elements_per_block_out = n_local_per_block;
  *elements_per_block_out = elements_per_block;
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_DISTRIBUTED_REORDER_IMPL_HPP
