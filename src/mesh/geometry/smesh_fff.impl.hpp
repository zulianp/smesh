#ifndef SMESH_FFF_IMPL_HPP
#define SMESH_FFF_IMPL_HPP

#include "smesh_jacobians.hpp"

#include "smesh_hex8_inline.hpp"
#include "smesh_tet4_inline.hpp"

namespace smesh {

// TET4

template <typename FFFType>
int tet4_fff_fill_simd(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff) {

  const geom_t *const SMESH_RESTRICT x = points[0];
  const geom_t *const SMESH_RESTRICT y = points[1];
  const geom_t *const SMESH_RESTRICT z = points[2];

#pragma omp parallel
  {
    static const int VEC_SIZE = 32;
    geom_t x0[VEC_SIZE];
    geom_t y0[VEC_SIZE];
    geom_t z0[VEC_SIZE];
    geom_t x1[VEC_SIZE];
    geom_t y1[VEC_SIZE];
    geom_t z1[VEC_SIZE];
    geom_t x2[VEC_SIZE];
    geom_t y2[VEC_SIZE];
    geom_t z2[VEC_SIZE];
    geom_t x3[VEC_SIZE];
    geom_t y3[VEC_SIZE];
    geom_t z3[VEC_SIZE];

#pragma omp for schedule(static)
    for (ptrdiff_t batch = 0; batch < nelements; batch += VEC_SIZE) {
      const ptrdiff_t len = MIN(VEC_SIZE, nelements - batch);

      for (ptrdiff_t l = 0; l < len; l++) {
        const ptrdiff_t e = batch + l;
        const ptrdiff_t i0 = elements[0][e];
        const ptrdiff_t i1 = elements[1][e];
        const ptrdiff_t i2 = elements[2][e];
        const ptrdiff_t i3 = elements[3][e];

        x0[l] = x[i0];
        y0[l] = y[i0];
        z0[l] = z[i0];

        x1[l] = x[i1];
        y1[l] = y[i1];
        z1[l] = z[i1];

        x2[l] = x[i2];
        y2[l] = y[i2];
        z2[l] = z[i2];

        x3[l] = x[i3];
        y3[l] = y[i3];
        z3[l] = z[i3];
      }

      FFFType *const SMESH_RESTRICT fff_0 = &fff[0][batch];
      FFFType *const SMESH_RESTRICT fff_1 = &fff[1][batch];
      FFFType *const SMESH_RESTRICT fff_2 = &fff[2][batch];
      FFFType *const SMESH_RESTRICT fff_3 = &fff[3][batch];
      FFFType *const SMESH_RESTRICT fff_4 = &fff[4][batch];
      FFFType *const SMESH_RESTRICT fff_5 = &fff[5][batch];

#pragma omp simd
      for (ptrdiff_t l = 0; l < len; l++) {
        tet4_fff(x0[l], x1[l], x2[l], x3[l], y0[l], y1[l], y2[l], y3[l], z0[l],
                 z1[l], z2[l], z3[l], &fff_0[l], &fff_1[l], &fff_2[l],
                 &fff_3[l], &fff_4[l], &fff_5[l]);
      }
    }
  }

  return SMESH_SUCCESS;
}

template <typename FFFType>
int tet4_fff_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff) {

  if (stride == 1) {
    return tet4_fff_fill_simd(nelements, elements, points, fff);
  }

  const geom_t *const SMESH_RESTRICT x = points[0];
  const geom_t *const SMESH_RESTRICT y = points[1];
  const geom_t *const SMESH_RESTRICT z = points[2];

#pragma omp parallel
  {

#pragma omp for schedule(static)
    for (ptrdiff_t e = 0; e < nelements; e++) {
      const ptrdiff_t i0 = elements[0][e];
      const ptrdiff_t i1 = elements[1][e];
      const ptrdiff_t i2 = elements[2][e];
      const ptrdiff_t i3 = elements[3][e];

      const ptrdiff_t idx = e * stride;

      tet4_fff(x[i0], x[i1], x[i2], x[i3], y[i0], y[i1], y[i2], y[i3], z[i0],
               z[i1], z[i2], z[i3], &fff[0][idx], &fff[1][idx], &fff[2][idx],
               &fff[3][idx], &fff[4][idx], &fff[5][idx]);
    }
  }
  return SMESH_SUCCESS;
}

// HEX8
template <typename FFFType>
int hex8_fff_fill_simd(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const geom_t qx, const geom_t qy, const geom_t qz,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff) {
  const geom_t *const SMESH_RESTRICT x = points[0];
  const geom_t *const SMESH_RESTRICT y = points[1];
  const geom_t *const SMESH_RESTRICT z = points[2];

#pragma omp parallel
  {
    static const int VEC_SIZE = 32;
    geom_t x0[VEC_SIZE];
    geom_t y0[VEC_SIZE];
    geom_t z0[VEC_SIZE];
    geom_t x1[VEC_SIZE];
    geom_t y1[VEC_SIZE];
    geom_t z1[VEC_SIZE];
    geom_t x2[VEC_SIZE];
    geom_t y2[VEC_SIZE];
    geom_t z2[VEC_SIZE];
    geom_t x3[VEC_SIZE];
    geom_t y3[VEC_SIZE];
    geom_t z3[VEC_SIZE];
    geom_t x4[VEC_SIZE];
    geom_t y4[VEC_SIZE];
    geom_t z4[VEC_SIZE];
    geom_t x5[VEC_SIZE];
    geom_t y5[VEC_SIZE];
    geom_t z5[VEC_SIZE];
    geom_t x6[VEC_SIZE];
    geom_t y6[VEC_SIZE];
    geom_t z6[VEC_SIZE];
    geom_t x7[VEC_SIZE];
    geom_t y7[VEC_SIZE];
    geom_t z7[VEC_SIZE];

#pragma omp for schedule(static)
    for (ptrdiff_t batch = 0; batch < nelements; batch += VEC_SIZE) {
      const ptrdiff_t len = MIN(VEC_SIZE, nelements - batch);

      for (ptrdiff_t l = 0; l < len; l++) {
        const ptrdiff_t e = batch + l;
        const ptrdiff_t i0 = elements[0][e];
        const ptrdiff_t i1 = elements[1][e];
        const ptrdiff_t i2 = elements[2][e];
        const ptrdiff_t i3 = elements[3][e];
        const ptrdiff_t i4 = elements[4][e];
        const ptrdiff_t i5 = elements[5][e];
        const ptrdiff_t i6 = elements[6][e];
        const ptrdiff_t i7 = elements[7][e];

        x0[l] = x[i0];
        y0[l] = y[i0];
        z0[l] = z[i0];

        x1[l] = x[i1];
        y1[l] = y[i1];
        z1[l] = z[i1];

        x2[l] = x[i2];
        y2[l] = y[i2];
        z2[l] = z[i2];

        x3[l] = x[i3];
        y3[l] = y[i3];
        z3[l] = z[i3];

        x4[l] = x[i4];
        y4[l] = y[i4];
        z4[l] = z[i4];

        x5[l] = x[i5];
        y5[l] = y[i5];
        z5[l] = z[i5];

        x6[l] = x[i6];
        y6[l] = y[i6];
        z6[l] = z[i6];

        x7[l] = x[i7];
        y7[l] = y[i7];
        z7[l] = z[i7];
      }

      FFFType *const SMESH_RESTRICT fff_0 = &fff[0][batch];
      FFFType *const SMESH_RESTRICT fff_1 = &fff[1][batch];
      FFFType *const SMESH_RESTRICT fff_2 = &fff[2][batch];
      FFFType *const SMESH_RESTRICT fff_3 = &fff[3][batch];
      FFFType *const SMESH_RESTRICT fff_4 = &fff[4][batch];
      FFFType *const SMESH_RESTRICT fff_5 = &fff[5][batch];

#pragma omp simd
      for (ptrdiff_t l = 0; l < len; l++) {
        hex8_fff(x0[l], x1[l], x2[l], x3[l], x4[l], x5[l], x6[l], x7[l],
                 //
                 y0[l], y1[l], y2[l], y3[l], y4[l], y5[l], y6[l], y7[l],
                 //
                 z0[l], z1[l], z2[l], z3[l], z4[l], z5[l], z6[l], z7[l],
                 //
                 qx, qy, qz,
                 //
                 &fff_0[l], &fff_1[l], &fff_2[l], &fff_3[l], &fff_4[l],
                 &fff_5[l]);
      }
    }
  }
  return SMESH_SUCCESS;
}

template <typename FFFType>
int hex8_fff_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const geom_t qx, const geom_t qy, const geom_t qz, const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff) {
  if (stride == 1) {
    return hex8_fff_fill_simd(nelements, elements, points, qx, qy, qz, fff);
  }

  const geom_t *const SMESH_RESTRICT x = points[0];
  const geom_t *const SMESH_RESTRICT y = points[1];
  const geom_t *const SMESH_RESTRICT z = points[2];

#pragma omp parallel
  {
#pragma omp for schedule(static)
    for (ptrdiff_t e = 0; e < nelements; e++) {
      const ptrdiff_t i0 = elements[0][e];
      const ptrdiff_t i1 = elements[1][e];
      const ptrdiff_t i2 = elements[2][e];
      const ptrdiff_t i3 = elements[3][e];
      const ptrdiff_t i4 = elements[4][e];
      const ptrdiff_t i5 = elements[5][e];
      const ptrdiff_t i6 = elements[6][e];
      const ptrdiff_t i7 = elements[7][e];

      const ptrdiff_t idx = e * stride;

      hex8_fff(x[i0], x[i1], x[i2], x[i3], x[i4], x[i5], x[i6], x[i7],
               //
               y[i0], y[i1], y[i2], y[i3], y[i4], y[i5], y[i6], y[i7],
               //
               z[i0], z[i1], z[i2], z[i3], z[i4], z[i5], z[i6], z[i7],
               //
               qx, qy, qz,
               //
               &fff[0][idx], &fff[1][idx], &fff[2][idx], &fff[3][idx],
               &fff[4][idx], &fff[5][idx]);
    }
  }
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_FFF_IMPL_HPP
