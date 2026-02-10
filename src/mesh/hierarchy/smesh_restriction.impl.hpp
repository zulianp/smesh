#ifndef SMESH_RESTRICTION_IMPL_HPP
#define SMESH_RESTRICTION_IMPL_HPP

#include "smesh_restriction.hpp"

namespace smesh {

template <typename idx_t, typename T>
int hierarchical_restriction(
    const enum ElemType from_element, const enum ElemType to_element,
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT e2n_count, const int vec_size,
    const T *const SMESH_RESTRICT from, T *const SMESH_RESTRICT to) {
  if (to_element == TET4 &&
      (from_element == TET10 || from_element == MACRO_TET4)) {
#pragma omp parallel for
    for (ptrdiff_t e = 0; e < nelements; e++) {
      // P1
      const ptrdiff_t i0 = elements[0][e];
      const ptrdiff_t i1 = elements[1][e];
      const ptrdiff_t i2 = elements[2][e];
      const ptrdiff_t i3 = elements[3][e];

      // P2
      const ptrdiff_t i4 = elements[4][e];
      const ptrdiff_t i5 = elements[5][e];
      const ptrdiff_t i6 = elements[6][e];
      const ptrdiff_t i7 = elements[7][e];
      const ptrdiff_t i8 = elements[8][e];
      const ptrdiff_t i9 = elements[9][e];

      for (int v = 0; v < vec_size; v++) {
        const ptrdiff_t ii0 = i0 * vec_size + v;
        const ptrdiff_t ii1 = i1 * vec_size + v;
        const ptrdiff_t ii2 = i2 * vec_size + v;
        const ptrdiff_t ii3 = i3 * vec_size + v;
        const ptrdiff_t ii4 = i4 * vec_size + v;
        const ptrdiff_t ii5 = i5 * vec_size + v;
        const ptrdiff_t ii6 = i6 * vec_size + v;
        const ptrdiff_t ii7 = i7 * vec_size + v;
        const ptrdiff_t ii8 = i8 * vec_size + v;
        const ptrdiff_t ii9 = i9 * vec_size + v;

        const T to0 = from[ii0] / e2n_count[i0] +
                      from[ii4] * (0.5 / e2n_count[i4]) +
                      from[ii6] * (0.5 / e2n_count[i6]) +
                      from[ii7] * (0.5 / e2n_count[i7]);

        const T to1 = from[ii1] / e2n_count[i1] +
                      from[ii5] * (0.5 / e2n_count[i5]) +
                      from[ii8] * (0.5 / e2n_count[i8]) +
                      from[ii4] * (0.5 / e2n_count[i4]);

        const T to2 = from[ii2] / e2n_count[i2] +
                      from[ii9] * (0.5 / e2n_count[i9]) +
                      from[ii5] * (0.5 / e2n_count[i5]) +
                      from[ii6] * (0.5 / e2n_count[i6]);

        const T to3 = from[ii3] / e2n_count[i3] +
                      from[ii7] * (0.5 / e2n_count[i7]) +
                      from[ii8] * (0.5 / e2n_count[i8]) +
                      from[ii9] * (0.5 / e2n_count[i9]);

#pragma omp atomic update
        to[ii0] += to0;

#pragma omp atomic update
        to[ii1] += to1;

#pragma omp atomic update
        to[ii2] += to2;

#pragma omp atomic update
        to[ii3] += to3;
      }
    }
  } else if (to_element == TRI3 &&
             (from_element == TRI6 || from_element == MACRO_TRI3)) {
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
        const ptrdiff_t ii0 = i0 * vec_size + v;
        const ptrdiff_t ii1 = i1 * vec_size + v;
        const ptrdiff_t ii2 = i2 * vec_size + v;
        const ptrdiff_t ii3 = i3 * vec_size + v;
        const ptrdiff_t ii4 = i4 * vec_size + v;
        const ptrdiff_t ii5 = i5 * vec_size + v;

        const T to0 = from[ii0] / e2n_count[i0] +
                      from[ii3] * (0.5 / e2n_count[i3]) +
                      from[ii5] * (0.5 / e2n_count[i5]);
        const T to1 = from[ii1] / e2n_count[i1] +
                      from[ii3] * (0.5 / e2n_count[i3]) +
                      from[ii4] * (0.5 / e2n_count[i4]);
        const T to2 = from[ii2] / e2n_count[i2] +
                      from[ii4] * (0.5 / e2n_count[i4]) +
                      from[ii5] * (0.5 / e2n_count[i5]);

#pragma omp atomic update
        to[ii0] += to0;

#pragma omp atomic update
        to[ii1] += to1;

#pragma omp atomic update
        to[ii2] += to2;
      }
    }
  } else {

    SMESH_ERROR(
        "Unsupported element pair for hierarchical_restriction %d, %d\n",
        from_element, to_element);

    return SMESH_FAILURE;
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_RESTRICTION_IMPL_HPP