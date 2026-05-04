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

template <typename geom_t>
int Hilbert3ElementOrdering<geom_t>::operator()(
    const ptrdiff_t n_points, const geom_t *const SMESH_RESTRICT x,
    const geom_t *const SMESH_RESTRICT y, const geom_t *const SMESH_RESTRICT z,
    const geom_t x_min, const geom_t x_max, const geom_t y_min,
    const geom_t y_max, const geom_t z_min, const geom_t z_max,
    u32 *const SMESH_RESTRICT encoding) const {
  return encode_hilbert3(n_points, x, y, z, x_min, x_max, y_min, y_max, z_min,
                         z_max, encoding);
}

template <typename idx_t, typename geom_t, typename Ordering>
int distributed_reorder_elements(
    MPI_Comm comm, const int nnodesxelem, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    large_idx_t *const SMESH_RESTRICT sorted_ids,
    Ordering ordering) {
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

  std::vector<idx_t> element_nodes((size_t)n_local_elements *
                                   (size_t)nnodesxelem);
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    for (int d = 0; d < nnodesxelem; ++d) {
      element_nodes[(size_t)e * (size_t)nnodesxelem + (size_t)d] =
          elements[d][e];
    }
  }

  std::vector<geom_t> element_x(element_nodes.size());
  std::vector<geom_t> element_y(element_nodes.size());
  std::vector<geom_t> element_z(element_nodes.size());
  if (gather_mapped_field(comm, (ptrdiff_t)element_nodes.size(), n_global_nodes,
                          element_nodes.data(), mpi_type<geom_t>(), points[0],
                          element_x.data()) != SMESH_SUCCESS ||
      gather_mapped_field(comm, (ptrdiff_t)element_nodes.size(), n_global_nodes,
                          element_nodes.data(), mpi_type<geom_t>(), points[1],
                          element_y.data()) != SMESH_SUCCESS ||
      gather_mapped_field(comm, (ptrdiff_t)element_nodes.size(), n_global_nodes,
                          element_nodes.data(), mpi_type<geom_t>(), points[2],
                          element_z.data()) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  std::vector<geom_t> cx((size_t)n_local_elements);
  std::vector<geom_t> cy((size_t)n_local_elements);
  std::vector<geom_t> cz((size_t)n_local_elements);
  geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max()};
  geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest()};

  const geom_t inv_nnodesxelem = geom_t(1) / geom_t(nnodesxelem);
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    geom_t x = 0;
    geom_t y = 0;
    geom_t z = 0;
    const size_t element_offset = (size_t)e * (size_t)nnodesxelem;
    for (int d = 0; d < nnodesxelem; ++d) {
      const size_t node = element_offset + (size_t)d;
      x += element_x[node];
      y += element_y[node];
      z += element_z[node];
    }

    x *= inv_nnodesxelem;
    y *= inv_nnodesxelem;
    z *= inv_nnodesxelem;
    cx[(size_t)e] = x;
    cy[(size_t)e] = y;
    cz[(size_t)e] = z;

    local_min[0] = std::min(local_min[0], x);
    local_min[1] = std::min(local_min[1], y);
    local_min[2] = std::min(local_min[2], z);
    local_max[0] = std::max(local_max[0], x);
    local_max[1] = std::max(local_max[1], y);
    local_max[2] = std::max(local_max[2], z);
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
                     mpi_type<large_idx_t>(), sorted_keys.data(),
                     sorted_ids, static_cast<int>(n_sorted_elements),
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
      nullptr, nnodesxelem_out, n_global_elements_out,
      n_owned_elements_out, n_shared_elements_out, n_ghost_elements_out,
      element_mapping_out, aura_element_mapping_out, elements_out,
      spatial_dim_out, n_global_nodes_out, n_owned_nodes_out,
      n_shared_nodes_out, n_ghost_nodes_out, n_aura_nodes_out,
      node_mapping_out, points_out, node_owner_out, node_offsets_out,
      ghosts_out);
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
        n_aura_elements > 0
            ? (smesh::large_idx_t *)SMESH_ALLOC(
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

} // namespace smesh

#endif // SMESH_DISTRIBUTED_REORDER_IMPL_HPP
