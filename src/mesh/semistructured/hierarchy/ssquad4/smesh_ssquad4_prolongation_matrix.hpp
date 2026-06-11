#ifndef SMESH_SSQUAD4_PROLONGATION_MATRIX_HPP
#define SMESH_SSQUAD4_PROLONGATION_MATRIX_HPP

#include "smesh_base.hpp"

#include <stddef.h>

namespace smesh {

    template <typename IDXType, typename RowPtrType, typename ColIdxType, typename T>
    int ssquad4_prolongation_matrix_symb(const int                                                 level,
                                         const ptrdiff_t                                           nelements,
                                         const ptrdiff_t                                           nnodes,
                                         const IDXType *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         RowPtrType *const SMESH_RESTRICT                          rowptr);

    template <typename IDXType, typename RowPtrType, typename ColIdxType, typename T>
    int ssquad4_prolongation_matrix_fill(const int                                                 level,
                                         const ptrdiff_t                                           nelements,
                                         const ptrdiff_t                                           nnodes,
                                         const IDXType *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         const RowPtrType *const SMESH_RESTRICT                    rowptr,
                                         ColIdxType *const SMESH_RESTRICT                          colidx,
                                         T *const SMESH_RESTRICT                                   values);

}  // namespace smesh

#endif  // SMESH_SSQUAD4_PROLONGATION_MATRIX_HPP
