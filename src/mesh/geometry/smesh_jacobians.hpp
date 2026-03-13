#ifndef SMESH_JACOBIANS_HPP
#define SMESH_JACOBIANS_HPP

#include "smesh_types.hpp"

namespace smesh {

// TODO: once we need 2D
// template <typename AdjugateType, typename DeterminantType>
// int tri3_adjugate_fill(
//     const ptrdiff_t nelements,
//     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
//     const ptrdiff_t stride,
//     AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
//     DeterminantType *const SMESH_RESTRICT determinant);

// template <typename AdjugateType, typename DeterminantType>
// int quad4_adjugate_fill(
//     const ptrdiff_t nelements,
//     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
//     const ptrdiff_t stride,
//     AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
//     DeterminantType *const SMESH_RESTRICT determinant);

template <typename AdjugateType, typename DeterminantType>
int tet4_adjugate_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t stride,
    AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
    DeterminantType *const SMESH_RESTRICT determinant);

template <typename AdjugateType, typename DeterminantType>
int hex8_adjugate_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const geom_t qx, const geom_t qy, const geom_t qz, const ptrdiff_t stride,
    AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
    DeterminantType *const SMESH_RESTRICT determinant);
} // namespace smesh

#endif // SMESH_JACOBIANS_HPP