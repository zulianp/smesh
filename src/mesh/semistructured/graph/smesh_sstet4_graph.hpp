#ifndef SMESH_SSTET4_GRAPH_HPP
#define SMESH_SSTET4_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

template <typename idx_t>
int sstet4_generate_elements(const int                                               L,
                             const ptrdiff_t                                         m_nelements,
                             const ptrdiff_t                                         m_nnodes,
                             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
                             idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
                             ptrdiff_t                                              *n_unique_nodes_out,
                             ptrdiff_t                                              *interior_start_out);

template <typename idx_t>
int sstet4_generate_elements_blocks(
        const int                                                                     L,
        const ptrdiff_t                                                               n_blocks,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               m_nnodes,
        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        const element_idx_t *const *const                                             hft,
        const block_idx_t *const *const                                               hnbb,
        ptrdiff_t                                                                    *n_unique_nodes_out,
        ptrdiff_t                                                                    *interior_start_out);

template <typename idx_t>
int sstet4_hierarchical_renumbering(const int                                         L,
                                    const int                                         nlevels,
                                    int *const                                        levels,
                                    const ptrdiff_t                                   nelements,
                                    const ptrdiff_t                                   nnodes,
                                    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                    idx_t *const SMESH_RESTRICT                       node_mapping,
                                    const bool                                        preserve_corner_ordering);

template <typename idx_t>
int sstet4_hierarchical_renumbering_blocks(
        const int                                                                     L,
        const int                                                                     nlevels,
        int *const                                                                    levels,
        const ptrdiff_t                                                               n_blocks,
        const ptrdiff_t *const                                                        n_elements,
        const ptrdiff_t                                                               nnodes,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT       elements,
        idx_t *const SMESH_RESTRICT                                                   node_mapping,
        const bool                                                                    preserve_corner_ordering);

}  // namespace smesh

#endif  // SMESH_SSTET4_GRAPH_HPP
