#ifndef CU_SSHEX8_RESTRICTION_CUH
#define CU_SSHEX8_RESTRICTION_CUH

#include <stddef.h>
#include "smesh_cuda_base.cuh"
#include "smesh_types.hpp"

namespace smesh {

    int cu_sshex8_hierarchical_restriction(const int                                               level,
                                           const ptrdiff_t                                         nelements,
                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                           const uint16_t *const SMESH_RESTRICT element_to_node_incidence_count,
                                           const int                            vec_size,
                                           const enum PrimitiveType             from_type,
                                           const ptrdiff_t                      from_stride,
                                           const void *const SMESH_RESTRICT     from,
                                           const enum PrimitiveType             to_type,
                                           const ptrdiff_t                      to_stride,
                                           void *const SMESH_RESTRICT           to,
                                           void                                *stream);

    int cu_sshex8_restrict(const ptrdiff_t                                         nelements,
                           const int                                               from_level,
                           const int                                               from_level_stride,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                           const uint16_t *const SMESH_RESTRICT                    from_element_to_node_incidence_count,
                           const int                                               to_level,
                           const int                                               to_level_stride,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                           const int                                               vec_size,
                           const enum PrimitiveType                                from_type,
                           const ptrdiff_t                                         from_stride,
                           const void *const SMESH_RESTRICT                        from,
                           const enum PrimitiveType                                to_type,
                           const ptrdiff_t                                         to_stride,
                           void *const SMESH_RESTRICT                              to,
                           void                                                   *stream);

}  // namespace smesh
#endif  // CU_SSHEX8_RESTRICTION_CUH
