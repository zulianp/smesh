#ifndef SMESH_JACOBIANS_HPP
#define SMESH_JACOBIANS_HPP

#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

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

template <typename AdjugateType, typename DeterminantType>
int adjugate_fill(
    const enum ElemType element_type, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t stride,
    AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
    DeterminantType *const SMESH_RESTRICT determinant) {
  switch (element_type) {
  case TET4:
    return tet4_adjugate_fill(nelements, elements, points, stride, adjugate,
                              determinant);
  case HEX8:
    return hex8_adjugate_fill(nelements, elements, points, 0.5, 0.5, 0.5,
                              stride, adjugate, determinant);
  default:
    SMESH_ERROR("adjugate_fill: Unsupported element type: %s\n",
                type_to_string(element_type));
    return SMESH_FAILURE;
  }
}
} // namespace smesh

#endif // SMESH_JACOBIANS_HPP