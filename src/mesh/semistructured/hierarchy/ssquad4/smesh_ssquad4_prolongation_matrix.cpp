#include "smesh_ssquad4_prolongation_matrix.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_SYMB(IDX_T, ROWPTR_T, COLIDX_T, T)                                \
    template int ssquad4_prolongation_matrix_symb<IDX_T, ROWPTR_T, COLIDX_T, T>(                                                  \
            const int                                               level,                                                        \
            const ptrdiff_t                                         nelements,                                                    \
            const ptrdiff_t                                         nnodes,                                                       \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,                                                     \
            ROWPTR_T *const SMESH_RESTRICT                          rowptr)

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_FILL(IDX_T, ROWPTR_T, COLIDX_T, T)                                \
    template int ssquad4_prolongation_matrix_fill<IDX_T, ROWPTR_T, COLIDX_T, T>(                                                  \
            const int                                               level,                                                        \
            const ptrdiff_t                                         nelements,                                                    \
            const ptrdiff_t                                         nnodes,                                                       \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,                                                     \
            const ROWPTR_T *const SMESH_RESTRICT                    rowptr,                                                       \
            COLIDX_T *const SMESH_RESTRICT                          colidx,                                                       \
            T *const SMESH_RESTRICT                                 values)

namespace smesh {
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_SYMB(i32, i32, i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_SYMB(i32, i32, i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_SYMB(i64, i64, i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_SYMB(i64, i64, i64, f64);

    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_FILL(i32, i32, i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_FILL(i32, i32, i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_FILL(i64, i64, i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_FILL(i64, i64, i64, f64);
}  // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_SYMB
#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_PROLONGATION_MATRIX_FILL
