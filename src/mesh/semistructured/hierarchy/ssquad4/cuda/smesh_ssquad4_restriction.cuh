#ifndef SMESH_SSQUAD4_RESTRICTION_CUH
#define SMESH_SSQUAD4_RESTRICTION_CUH

#include "smesh_base.hpp"
#include "smesh_types.hpp"
#include <stddef.h>

namespace smesh {

int cu_ssquad4_hierarchical_restriction(
    const int level, const ptrdiff_t nelements, const ptrdiff_t stride,
    const idx_t *const SMESH_RESTRICT elements,
    const uint16_t *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const enum PrimitiveType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum PrimitiveType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream);

int cu_ssquad4_restrict(
    const ptrdiff_t nelements, const int from_level,
    const int from_level_stride, idx_t **const SMESH_RESTRICT from_elements,
    const uint16_t *const SMESH_RESTRICT from_element_to_node_incidence_count,
    const int to_level, const int to_level_stride,
    idx_t **const SMESH_RESTRICT to_elements, const int vec_size,
    const enum PrimitiveType from_type, const ptrdiff_t from_stride,
    const void *const SMESH_RESTRICT from, const enum PrimitiveType to_type,
    const ptrdiff_t to_stride, void *const SMESH_RESTRICT to, void *stream);

} // namespace smesh

#endif // CU_SSQUAD4_INTERPOLATE_H