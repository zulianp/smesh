#ifndef SMESH_DISTRIBUTED_REORDER_IMPL_HPP
#define SMESH_DISTRIBUTED_REORDER_IMPL_HPP

#include "smesh_distributed_reorder.hpp"

#include "smesh_decompose.hpp"
#include "smesh_distributed_aura.hpp"
#include "smesh_sfc.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
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
    Ordering ordering) {
  if (n_local_elements == 0) {
    return SMESH_SUCCESS;
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
  std::vector<large_idx_t> sorted_ids((size_t)n_sorted_elements);
  if (MPI_Sort_bykey(send_keys.data(), send_ids.data(),
                     static_cast<int>(n_local_elements), mpi_type<u32>(),
                     mpi_type<large_idx_t>(), sorted_keys.data(),
                     sorted_ids.data(), static_cast<int>(n_sorted_elements),
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
                            sorted_ids.data(), mpi_type<idx_t>(), elements[d],
                            sorted_elements.data()) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    std::memcpy(elements[d], sorted_elements.data(),
                (size_t)n_sorted_elements * sizeof(idx_t));
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_DISTRIBUTED_REORDER_IMPL_HPP
