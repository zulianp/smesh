#ifndef SMESH_DISTRIBUTED_REORDER_HPP
#define SMESH_DISTRIBUTED_REORDER_HPP

#include "smesh_distributed_base.hpp"
#include "smesh_path.hpp"
#include "smesh_sfc.hpp"
#include "smesh_types.hpp"

#include <cstddef>
#include <mpi.h>

namespace smesh {

template <typename geom_t>
using ElementOrdering = int (*)(const ptrdiff_t n_points,
                                const geom_t *const SMESH_RESTRICT x,
                                const geom_t *const SMESH_RESTRICT y,
                                const geom_t *const SMESH_RESTRICT z,
                                const geom_t x_min, const geom_t x_max,
                                const geom_t y_min, const geom_t y_max,
                                const geom_t z_min, const geom_t z_max,
                                u32 *const SMESH_RESTRICT encoding);

template <typename idx_t, typename geom_t,
          typename Ordering = ElementOrdering<geom_t>>
int distributed_reorder_elements(
    MPI_Comm comm, const int nnodesxelem, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    large_idx_t *const SMESH_RESTRICT sorted_ids,
    Ordering ordering = encode_hilbert3<geom_t>);

template <typename idx_t, typename geom_t, typename global_idx_t,
          typename Ordering = ElementOrdering<geom_t>>
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
    Ordering ordering = encode_hilbert3<geom_t>);
} // namespace smesh

#endif // SMESH_DISTRIBUTED_REORDER_HPP
