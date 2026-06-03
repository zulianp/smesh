#ifndef SMESH_SSQUAD4_GRAPH_HPP
#define SMESH_SSQUAD4_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int quad4_build_edge_graph_from_n2e(
    const ptrdiff_t nelements, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex, count_t **out_rowptr,
    idx_t **out_colidx);

} // namespace smesh

#endif
