#ifndef CU_SSHEX8_PROLONGATION_CUH
#define CU_SSHEX8_PROLONGATION_CUH

#include "smesh_cuda_base.cuh"
#include "smesh_types.hpp"
#include <stddef.h>

namespace smesh {

int cu_sshex8_hierarchical_prolongation(
    const int level, const ptrdiff_t nelements,
    idx_t **const SMESH_RESTRICT elements, const int vec_size,
    const enum PrimitiveType from_type, const ptrdiff_t from_stride,
    const void *const SMESH_RESTRICT from, const enum PrimitiveType to_type,
    const ptrdiff_t to_stride, void *const SMESH_RESTRICT to, void *stream);

int cu_sshex8_prolongate(
    const ptrdiff_t nelements, const int from_level,
    const int from_level_stride, idx_t **const SMESH_RESTRICT from_elements,
    const int to_level, const int to_level_stride,
    idx_t **const SMESH_RESTRICT to_elements, const int vec_size,
    const enum PrimitiveType from_type, const ptrdiff_t from_stride,
    const void *const SMESH_RESTRICT from, const enum PrimitiveType to_type,
    const ptrdiff_t to_stride, void *const SMESH_RESTRICT to, void *stream);

} // namespace smesh

#endif // CU_SSHEX8_PROLONGATION_CUH
