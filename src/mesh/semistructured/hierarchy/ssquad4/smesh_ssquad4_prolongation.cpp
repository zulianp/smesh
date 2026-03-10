#include "smesh_ssquad4_prolongation.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION(IDX_T, T)                                                               \
    template int ssquad4_prolongate<IDX_T, T>(const ptrdiff_t                                         nelements,                \
                                              const int                                               from_level,               \
                                              const int                                               from_level_stride,        \
                                              const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,            \
                                              const int                                               to_level,                 \
                                              const int                                               to_level_stride,          \
                                              const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,              \
                                              const int                                               vec_size,                 \
                                              const T *const SMESH_RESTRICT                           from,                     \
                                              T *const SMESH_RESTRICT                                 to);                                                      \
    template int ssquad4_hierarchical_prolongation<IDX_T, T>(int                                                     level,     \
                                                             const ptrdiff_t                                         nelements, \
                                                             const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,  \
                                                             const int                                               vec_size,  \
                                                             const T *const SMESH_RESTRICT                           from,      \
                                                             T *const SMESH_RESTRICT                                 to)

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_NNZ(IDX_T, COUNT_T)                                                  \
    template int ssquad4_prolongation_crs_nnz<IDX_T, COUNT_T>(const int                                               level,     \
                                                              const ptrdiff_t                                         nelements, \
                                                              const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,  \
                                                              const ptrdiff_t                                         to_nnodes, \
                                                              COUNT_T *const SMESH_RESTRICT                           rowptr)

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_FILL(IDX_T, COUNT_T, T)                                               \
    template int ssquad4_prolongation_crs_fill<IDX_T, COUNT_T, T>(const int                                               level,     \
                                                                   const ptrdiff_t                                         nelements, \
                                                                   const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,  \
                                                                   const ptrdiff_t                                         to_nnodes, \
                                                                   COUNT_T *const SMESH_RESTRICT                           rowptr,    \
                                                                   IDX_T *const SMESH_RESTRICT                             colidx,    \
                                                                   T *const SMESH_RESTRICT                                 values)

namespace smesh {
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION(i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION(i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION(i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION(i64, f64);

    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_NNZ(i32, i32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_NNZ(i64, i64);

    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_FILL(i32, i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_FILL(i32, i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_FILL(i64, i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_FILL(i64, i64, f64);
}  // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION
#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_NNZ
#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_CRS_FILL
