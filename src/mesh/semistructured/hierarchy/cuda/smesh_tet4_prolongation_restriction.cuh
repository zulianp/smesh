#ifndef SMESH_TET4_PROLONGATION_RESTRICTION_CUH
#define SMESH_TET4_PROLONGATION_RESTRICTION_CUH

#include <stddef.h>

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

int cu_tet4_to_macrotet4_prolongation(
    const ptrdiff_t coarse_nnodes,
    const count_t *const SMESH_RESTRICT coarse_rowptr,
    const idx_t *const SMESH_RESTRICT coarse_colidx,
    const idx_t *const SMESH_RESTRICT fine_node_map, const int vec_size,
    const enum PrimitiveType from_type, const void *const SMESH_RESTRICT from,
    const enum PrimitiveType to_type, void *const SMESH_RESTRICT to,
    void *stream);

int cu_macrotet4_to_tet4_restriction(
    const ptrdiff_t coarse_nnodes,
    const count_t *const SMESH_RESTRICT coarse_rowptr,
    const idx_t *const SMESH_RESTRICT coarse_colidx,
    const idx_t *const SMESH_RESTRICT fine_node_map, const int vec_size,
    const enum PrimitiveType from_type, const void *const SMESH_RESTRICT from,
    const enum PrimitiveType to_type, void *const SMESH_RESTRICT to,
    void *stream);

// Element-based (more generic)

int cu_macrotet4_to_tet4_prolongation_element_based(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int vec_size, const enum PrimitiveType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum PrimitiveType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream);

int cu_macrotet4_to_tet4_restriction_element_based(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const uint16_t *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const enum PrimitiveType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum PrimitiveType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream);

} // namespace smesh

#endif // CUT_TET4_PROLONGATION_RESTRICTION_H
