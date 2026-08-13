#ifndef SMESH_DISTRIBUTED_REORDER_HPP
#define SMESH_DISTRIBUTED_REORDER_HPP

#include "smesh_distributed_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"
#include "smesh_sfc.hpp"
#include "smesh_types.hpp"

#include <cstddef>
#include <mpi.h>
#include <string>
#include <vector>

namespace smesh {

template <typename geom_t>
using OrderEncoder = int (*)(const ptrdiff_t n_points,
                                const geom_t *const SMESH_RESTRICT x,
                                const geom_t *const SMESH_RESTRICT y,
                                const geom_t *const SMESH_RESTRICT z,
                                const geom_t x_min, const geom_t x_max,
                                const geom_t y_min, const geom_t y_max,
                                const geom_t z_min, const geom_t z_max,
                                u32 *const SMESH_RESTRICT encoding);

template <typename idx_t, typename geom_t,
          typename Ordering = OrderEncoder<geom_t>>
int distributed_reorder_elements(
    MPI_Comm comm, const int nnodesxelem, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_global_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    large_idx_t *const SMESH_RESTRICT sorted_ids,
    Ordering ordering = encode_hilbert3<geom_t>);

/// Assign ownership of the union of block elements via SFC (or identity keys).
///
/// On entry, each block has file-order `rank_split(n_global_b)` local elements.
/// On exit, this rank owns `rank_split(sum n_global_b)` elements identified by
/// `sorted_concat_ids` (exclusive-scan concat ids). Connectivity is gathered
/// into a CSR e2n (`e2n_ptr` / `e2n_idx`) so blocks may have different nxe.
///
/// @param use_sfc  If false, sort keys are concat ids (file/concat order).
template <typename idx_t, typename geom_t,
          typename Ordering = OrderEncoder<geom_t>>
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
    ptrdiff_t **const SMESH_RESTRICT e2n_ptr_out,
    idx_t **const SMESH_RESTRICT e2n_idx_out, const bool use_sfc = true,
    Ordering ordering = encode_hilbert3<geom_t>);

template <typename idx_t, typename geom_t, typename global_idx_t,
          typename Ordering = OrderEncoder<geom_t>>
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

/// MPI read of a serial multi-block folder (`meta.yaml` + `blocks/<name>/`).
/// Heterogeneous nxe is allowed. Ownership uses an SFC over the union of
/// elements unless @p use_sfc is false (concat-id partition).
template <typename idx_t, typename geom_t, typename large_idx_t,
          typename Ordering = OrderEncoder<geom_t>>
int mesh_from_folder_multiblock(
    const MPI_Comm comm, const Path &folder,
    const std::vector<std::string> &block_names,
    const std::vector<enum ElemType> &element_types,
    ptrdiff_t *n_global_elements_out, ptrdiff_t *n_owned_elements_out,
    ptrdiff_t *n_shared_elements_out, ptrdiff_t *n_ghost_elements_out,
    large_idx_t **element_mapping_out, large_idx_t **aura_element_mapping_out,
    int **nxe_per_block_out, ptrdiff_t **n_local_elements_per_block_out,
    idx_t ****elements_per_block_out, int *spatial_dim_out,
    ptrdiff_t *n_global_nodes_out, ptrdiff_t *n_owned_nodes_out,
    ptrdiff_t *n_shared_nodes_out, ptrdiff_t *n_ghost_nodes_out,
    ptrdiff_t *n_aura_nodes_out, large_idx_t **node_mapping_out,
    geom_t ***points_out, int **node_owner_out, ptrdiff_t **node_offsets_out,
    idx_t **ghosts_out, const bool use_sfc = true,
    Ordering ordering = encode_hilbert3<geom_t>);
} // namespace smesh

#endif // SMESH_DISTRIBUTED_REORDER_HPP
