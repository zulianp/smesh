#include "smesh_sstet4_prolongation.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSTET4_PROLONGATION(IDX_T, T)                                                               \
    template int sstet4_prolongate<IDX_T, T>(const ptrdiff_t                                         nelements,                \
                                             const int                                               from_level,               \
                                             const int                                               from_level_stride,        \
                                             const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,            \
                                             const int                                               to_level,                 \
                                             const int                                               to_level_stride,          \
                                             const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,              \
                                             const int                                               vec_size,                 \
                                             const T *const SMESH_RESTRICT                           from,                     \
                                             T *const SMESH_RESTRICT                                 to);                      \
    template int sstet4_hierarchical_prolongation<IDX_T, T>(int                                                     level,     \
                                                            const ptrdiff_t                                         nelements, \
                                                            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,  \
                                                            const int                                               vec_size,  \
                                                            const T *const SMESH_RESTRICT                           from,      \
                                                            T *const SMESH_RESTRICT                                 to)

namespace smesh {
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_PROLONGATION(i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_PROLONGATION(i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_PROLONGATION(i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_PROLONGATION(i64, f64);
}  // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSTET4_PROLONGATION
