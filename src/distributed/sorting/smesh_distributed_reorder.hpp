#ifndef SMESH_DISTRIBUTED_REORDER_HPP
#define SMESH_DISTRIBUTED_REORDER_HPP

#include "smesh_distributed_base.hpp"
#include "smesh_types.hpp"

#include <cstddef>
#include <mpi.h>

namespace smesh {

template <typename geom_t> struct Hilbert3ElementOrdering {
  int operator()(const ptrdiff_t n_points, const geom_t *const SMESH_RESTRICT x,
                 const geom_t *const SMESH_RESTRICT y,
                 const geom_t *const SMESH_RESTRICT z,
                 const geom_t x_min, const geom_t x_max, const geom_t y_min,
                 const geom_t y_max, const geom_t z_min, const geom_t z_max,
                 u32 *const SMESH_RESTRICT encoding) const;
};

template <typename idx_t, typename geom_t,
          typename Ordering = Hilbert3ElementOrdering<geom_t>>
int distributed_sort_elements(
    MPI_Comm comm, const int nnodesxelem, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    Ordering ordering = Ordering());

template <typename idx_t, typename geom_t,
          typename Ordering = Hilbert3ElementOrdering<geom_t>>
int distributed_reorder_elements(
    MPI_Comm comm, const int nnodesxelem, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    Ordering ordering = Ordering());

} // namespace smesh

#endif // SMESH_DISTRIBUTED_REORDER_HPP
