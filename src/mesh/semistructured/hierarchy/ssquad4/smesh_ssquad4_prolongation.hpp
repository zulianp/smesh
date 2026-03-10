#ifndef SMESH_SSQUAD4_PROLONGATION_HPP
#define SMESH_SSQUAD4_PROLONGATION_HPP

#include <stddef.h>
#include "smesh_base.hpp"

namespace smesh {

    template <typename idx_t, typename T>
    int ssquad4_hierarchical_prolongation(int                                                     level,
                                          const ptrdiff_t                                         nelements,
                                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                          const int                                               vec_size,
                                          const T *const SMESH_RESTRICT                           from,
                                          T *const SMESH_RESTRICT                                 to);

    template <typename idx_t, typename T>
    int ssquad4_prolongate(const ptrdiff_t                                         nelements,
                           const int                                               from_level,
                           const int                                               from_level_stride,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                           const int                                               to_level,
                           const int                                               to_level_stride,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                           const int                                               vec_size,
                           const T *const SMESH_RESTRICT                           from,
                           T *const SMESH_RESTRICT                                 to);

    template <typename idx_t, typename count_t>
    int ssquad4_prolongation_crs_nnz(const int                                               level,
                                     const ptrdiff_t                                         nelements,
                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                     const ptrdiff_t                                         to_nnodes,
                                     count_t *const SMESH_RESTRICT                           rowptr);

    template <typename idx_t, typename count_t, typename T>
    int ssquad4_prolongation_crs_fill(const int                                               level,
                                      const ptrdiff_t                                         nelements,
                                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                      const ptrdiff_t                                         to_nnodes,
                                      count_t *const SMESH_RESTRICT                           rowptr,
                                      idx_t *const SMESH_RESTRICT                             colidx,
                                      T *const SMESH_RESTRICT                                 values);

}  // namespace smesh
#endif  // SMESH_SSQUAD4_INTERPOLATE_H
