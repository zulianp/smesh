#ifndef SMESH_RESTRICTION_CUH
#define SMESH_RESTRICTION_CUH

#include <stddef.h>

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t>
int cu_macrotet4_to_tet4_restriction_elemental(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const enum PrimitiveType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum PrimitiveType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream);

} // namespace smesh

#endif // SMESH_RESTRICTION_CUH
