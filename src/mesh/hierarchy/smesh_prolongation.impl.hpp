#ifndef SMESH_PROLONGATION_IMPL_HPP
#define SMESH_PROLONGATION_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {
template <typename idx_t, typename T>
int hierarchical_prolongation(
    const enum ElemType from_element, const enum ElemType to_element,
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int vec_size, const T *const SMESH_RESTRICT from,
    T *const SMESH_RESTRICT to) {
  if (from_element == TET4 &&
      (to_element == TET10 || to_element == MACRO_TET4)) {
#pragma omp parallel for
    for (ptrdiff_t e = 0; e < nelements; e++) {
      // P1
      const idx_t i0 = elements[0][e];
      const idx_t i1 = elements[1][e];
      const idx_t i2 = elements[2][e];
      const idx_t i3 = elements[3][e];

      // P2
      const idx_t i4 = elements[4][e];
      const idx_t i5 = elements[5][e];
      const idx_t i6 = elements[6][e];
      const idx_t i7 = elements[7][e];
      const idx_t i8 = elements[8][e];
      const idx_t i9 = elements[9][e];

      SMESH_ASSERT(i0 != i4);

      for (int v = 0; v < vec_size; v++) {
        to[i0 * vec_size + v] = from[i0 * vec_size + v];
        to[i1 * vec_size + v] = from[i1 * vec_size + v];
        to[i2 * vec_size + v] = from[i2 * vec_size + v];
        to[i3 * vec_size + v] = from[i3 * vec_size + v];

        to[i4 * vec_size + v] =
            0.5 * (from[i0 * vec_size + v] + from[i1 * vec_size + v]);
        to[i5 * vec_size + v] =
            0.5 * (from[i1 * vec_size + v] + from[i2 * vec_size + v]);
        to[i6 * vec_size + v] =
            0.5 * (from[i0 * vec_size + v] + from[i2 * vec_size + v]);
        to[i7 * vec_size + v] =
            0.5 * (from[i0 * vec_size + v] + from[i3 * vec_size + v]);
        to[i8 * vec_size + v] =
            0.5 * (from[i1 * vec_size + v] + from[i3 * vec_size + v]);
        to[i9 * vec_size + v] =
            0.5 * (from[i2 * vec_size + v] + from[i3 * vec_size + v]);
      }
    }
  } else if (from_element == TRI3 &&
             (to_element == TRI6 || to_element == MACRO_TRI3)) {
#pragma omp parallel for
    for (ptrdiff_t e = 0; e < nelements; e++) {
      // P1
      const idx_t i0 = elements[0][e];
      const idx_t i1 = elements[1][e];
      const idx_t i2 = elements[2][e];

      // P2
      const idx_t i3 = elements[3][e];
      const idx_t i4 = elements[4][e];
      const idx_t i5 = elements[5][e];

      for (int v = 0; v < vec_size; v++) {
        to[i0 * vec_size + v] = from[i0 * vec_size + v];
        to[i1 * vec_size + v] = from[i1 * vec_size + v];
        to[i2 * vec_size + v] = from[i2 * vec_size + v];

        to[i3 * vec_size + v] =
            0.5 * (from[i0 * vec_size + v] + from[i1 * vec_size + v]);
        to[i4 * vec_size + v] =
            0.5 * (from[i1 * vec_size + v] + from[i2 * vec_size + v]);
        to[i5 * vec_size + v] =
            0.5 * (from[i0 * vec_size + v] + from[i2 * vec_size + v]);
      }
    }
  } else {
    SMESH_ERROR(
        "Unsupported element pair for hierarchical_prolongation %d, %d\n",
        from_element, to_element);
    return SMESH_FAILURE;
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_PROLONGATION_IMPL_HPP
