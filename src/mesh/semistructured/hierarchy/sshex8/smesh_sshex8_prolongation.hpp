#ifndef SMESH_SSHEX8_PROLONGATION_HPP
#define SMESH_SSHEX8_PROLONGATION_HPP

#include <stddef.h>
#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

    template <typename idx_t, typename T>
    int sshex8_hierarchical_prolongation(int                                                     level,
                                         const ptrdiff_t                                         nelements,
                                         const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         const int                                               vec_size,
                                         const T *const SMESH_RESTRICT                           from,
                                         T *const SMESH_RESTRICT                                 to);

    template <typename idx_t, typename T>
    int sshex8_prolongate(const ptrdiff_t                                         nelements,
                          const int                                               from_level,
                          const int                                               from_level_stride,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                          const int                                               to_level,
                          const int                                               to_level_stride,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                          const int                                               vec_size,
                          const T *const SMESH_RESTRICT                           from,
                          T *const SMESH_RESTRICT                                 to);

}  // namespace smesh

#endif  // SMESH_SSHEX8_PROLONGATION_HPP
