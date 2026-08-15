#ifndef SMESH_SSMIXED_GRAPH_HPP
#define SMESH_SSMIXED_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

/// Mesh-wide unique SS numbering for mixed HEX8+TET4 blocks.
/// Shared edge interiors get the same ids in both families. Face and volume
/// interiors stay family-local. Does not modify the homogeneous HEX/TET kernels.
template <typename idx_t>
int ssmixed_hex_tet_generate_elements_blocks(
        const int                                                                     L,
        const ptrdiff_t                                                               n_blocks,
        const enum ElemType *const                                                    block_types,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               m_nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        const element_idx_t *const *const                                             hft,
        const block_idx_t *const *const                                               hnbb,
        ptrdiff_t                                                                    *n_unique_nodes_out,
        ptrdiff_t                                                                    *interior_start_out);

/// One node_mapping over HEX and TET lattices. Interleaves families per
/// hierarchical layer so shared edge nodes keep a single id. Do not call
/// sshex8 then sstet4 hierarchical_renumbering_blocks independently.
template <typename idx_t>
int ssmixed_hex_tet_hierarchical_renumbering_blocks(
        const int                                                                     L,
        const int                                                                     nlevels,
        int *const                                                                    levels,
        const ptrdiff_t                                                               n_blocks,
        const enum ElemType *const                                                    block_types,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               nnodes,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        idx_t *const SMESH_RESTRICT                                                   node_mapping,
        const bool                                                                    preserve_corner_ordering);

}  // namespace smesh

#endif
