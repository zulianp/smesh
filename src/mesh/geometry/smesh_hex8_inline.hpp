#ifndef SMESH_HEX8_INLINE_HPP
#define SMESH_HEX8_INLINE_HPP

#include "smesh_base.hpp"

#ifndef POW2
#define POW2(a) ((a) * (a))
#endif

namespace smesh {

template <typename T, typename FFFType>
static SMESH_INLINE SMESH_HOST_DEVICE void hex8_fff( // X
    const T px0, const T px1, const T px2, const T px3, const T px4,
    const T px5, const T px6, const T px7,
    // Y
    const T py0, const T py1, const T py2, const T py3, const T py4,
    const T py5, const T py6, const T py7,
    // Z
    const T pz0, const T pz1, const T pz2, const T pz3, const T pz4,
    const T pz5, const T pz6, const T pz7,
    // QPs
    const T qx, const T qy, const T qz,
    // Out
    FFFType *const SMESH_RESTRICT fff0, FFFType *const SMESH_RESTRICT fff1,
    FFFType *const SMESH_RESTRICT fff2, FFFType *const SMESH_RESTRICT fff3,
    FFFType *const SMESH_RESTRICT fff4, FFFType *const SMESH_RESTRICT fff5) {
  const T x0 = qy * qz;
  const T x1 = 1 - qz;
  const T x2 = qy * x1;
  const T x3 = 1 - qy;
  const T x4 = qz * x3;
  const T x5 = x1 * x3;
  const T x6 = x0 * pz6 - x0 * pz7 + x2 * pz2 - x2 * pz3 - x4 * pz4 + x4 * pz5 -
               x5 * pz0 + x5 * pz1;
  const T x7 = qx * qy;
  const T x8 = qx * x3;
  const T x9 = 1 - qx;
  const T x10 = qy * x9;
  const T x11 = x3 * x9;
  const T x12 = qx * qy * px6 + qx * x3 * px5 + qy * x9 * px7 - x10 * px3 -
                x11 * px0 + x3 * x9 * px4 - x7 * px2 - x8 * px1;
  const T x13 = qx * qz;
  const T x14 = qx * x1;
  const T x15 = qz * x9;
  const T x16 = x1 * x9;
  const T x17 = qx * qz * py6 + qx * x1 * py2 + qz * x9 * py7 + x1 * x9 * py3 -
                x13 * py5 - x14 * py1 - x15 * py4 - x16 * py0;
  const T x18 = x12 * x17;
  const T x19 = x0 * px6 - x0 * px7 + x2 * px2 - x2 * px3 - x4 * px4 +
                x4 * px5 - x5 * px0 + x5 * px1;
  const T x20 = qx * qy * py6 + qx * x3 * py5 + qy * x9 * py7 - x10 * py3 -
                x11 * py0 + x3 * x9 * py4 - x7 * py2 - x8 * py1;
  const T x21 = qx * qz * pz6 + qx * x1 * pz2 + qz * x9 * pz7 + x1 * x9 * pz3 -
                x13 * pz5 - x14 * pz1 - x15 * pz4 - x16 * pz0;
  const T x22 = x20 * x21;
  const T x23 = x0 * py6 - x0 * py7 + x2 * py2 - x2 * py3 - x4 * py4 +
                x4 * py5 - x5 * py0 + x5 * py1;
  const T x24 = qx * qy * pz6 + qx * x3 * pz5 + qy * x9 * pz7 - x10 * pz3 -
                x11 * pz0 + x3 * x9 * pz4 - x7 * pz2 - x8 * pz1;
  const T x25 = qx * qz * px6 + qx * x1 * px2 + qz * x9 * px7 + x1 * x9 * px3 -
                x13 * px5 - x14 * px1 - x15 * px4 - x16 * px0;
  const T x26 = x24 * x25;
  const T x27 = x12 * x21 * x23 + x17 * x19 * x24 - x18 * x6 - x19 * x22 +
                x20 * x25 * x6 - x23 * x26;
  const T x28 = -x18 + x20 * x25;
  const T x29 = (1 / POW2(x27));
  const T x30 = x12 * x21 - x26;
  const T x31 = x17 * x24 - x22;
  const T x32 = x12 * x23 - x19 * x20;
  const T x33 = x28 * x29;
  const T x34 = -x12 * x6 + x19 * x24;
  const T x35 = x29 * x30;
  const T x36 = x20 * x6 - x23 * x24;
  const T x37 = x29 * x31;
  const T x38 = x17 * x19 - x23 * x25;
  const T x39 = -x19 * x21 + x25 * x6;
  const T x40 = -x17 * x6 + x21 * x23;
  *fff0 =
      (FFFType)(x27 * (POW2(x28) * x29 + x29 * POW2(x30) + x29 * POW2(x31)));
  *fff1 = (FFFType)(x27 * (x32 * x33 + x34 * x35 + x36 * x37));
  *fff2 = (FFFType)(x27 * (x33 * x38 + x35 * x39 + x37 * x40));
  *fff3 =
      (FFFType)(x27 * (x29 * POW2(x32) + x29 * POW2(x34) + x29 * POW2(x36)));
  *fff4 =
      (FFFType)(x27 * (x29 * x32 * x38 + x29 * x34 * x39 + x29 * x36 * x40));
  *fff5 =
      (FFFType)(x27 * (x29 * POW2(x38) + x29 * POW2(x39) + x29 * POW2(x40)));
}

template <typename T>
static SMESH_INLINE SMESH_HOST_DEVICE void
aahex8_jac_diag(const T px0, const T px6, const T py0, const T py6, const T pz0,
                const T pz6, T *const SMESH_RESTRICT jac_diag0,
                T *const SMESH_RESTRICT jac_diag1,
                T *const SMESH_RESTRICT jac_diag2) {
  const T x0 = -py0 + py6;
  const T x1 = -pz0 + pz6;
  const T x2 = -px0 + px6;
  *jac_diag0 = x0 * x1 / x2;
  *jac_diag1 = x1 * x2 / x0;
  *jac_diag2 = x0 * x2 / x1;
}

template <typename T, typename AdjugateType, typename DeterminantType>
static SMESH_INLINE SMESH_HOST_DEVICE void hex8_adjugate_and_det(
    // X
    const T px0, const T px1, const T px2, const T px3, const T px4,
    const T px5, const T px6, const T px7,
    // Y
    const T py0, const T py1, const T py2, const T py3, const T py4,
    const T py5, const T py6, const T py7,
    // Z
    const T pz0, const T pz1, const T pz2, const T pz3, const T pz4,
    const T pz5, const T pz6, const T pz7,
    // QPs
    const T qx, const T qy, const T qz,
    // Out
    AdjugateType *const SMESH_RESTRICT adjugate0,
    AdjugateType *const SMESH_RESTRICT adjugate1,
    AdjugateType *const SMESH_RESTRICT adjugate2,
    AdjugateType *const SMESH_RESTRICT adjugate3,
    AdjugateType *const SMESH_RESTRICT adjugate4,
    AdjugateType *const SMESH_RESTRICT adjugate5,
    AdjugateType *const SMESH_RESTRICT adjugate6,
    AdjugateType *const SMESH_RESTRICT adjugate7,
    AdjugateType *const SMESH_RESTRICT adjugate8,
    DeterminantType *const SMESH_RESTRICT determinant) {
  T jacobian[9];
  {
    const T x0 = qy * qz;
    const T x1 = 1 - qz;
    const T x2 = qy * x1;
    const T x3 = 1 - qy;
    const T x4 = qz * x3;
    const T x5 = x1 * x3;
    const T x6 = qx * qz;
    const T x7 = qx * x1;
    const T x8 = 1 - qx;
    const T x9 = qz * x8;
    const T x10 = x1 * x8;
    const T x11 = qx * qy;
    const T x12 = qx * x3;
    const T x13 = qy * x8;
    const T x14 = x3 * x8;

    jacobian[0] = x0 * px6 - x0 * px7 + x2 * px2 - x2 * px3 - x4 * px4 +
                  x4 * px5 - x5 * px0 + x5 * px1;
    jacobian[1] = qx * qz * px6 + qx * x1 * px2 + qz * x8 * px7 +
                  x1 * x8 * px3 - x10 * px0 - x6 * px5 - x7 * px1 - x9 * px4;
    jacobian[2] = qx * qy * px6 + qx * x3 * px5 + qy * x8 * px7 - x11 * px2 -
                  x12 * px1 - x13 * px3 - x14 * px0 + x3 * x8 * px4;
    jacobian[3] = x0 * py6 - x0 * py7 + x2 * py2 - x2 * py3 - x4 * py4 +
                  x4 * py5 - x5 * py0 + x5 * py1;
    jacobian[4] = qx * qz * py6 + qx * x1 * py2 + qz * x8 * py7 +
                  x1 * x8 * py3 - x10 * py0 - x6 * py5 - x7 * py1 - x9 * py4;
    jacobian[5] = qx * qy * py6 + qx * x3 * py5 + qy * x8 * py7 - x11 * py2 -
                  x12 * py1 - x13 * py3 - x14 * py0 + x3 * x8 * py4;
    jacobian[6] = x0 * pz6 - x0 * pz7 + x2 * pz2 - x2 * pz3 - x4 * pz4 +
                  x4 * pz5 - x5 * pz0 + x5 * pz1;
    jacobian[7] = qx * qz * pz6 + qx * x1 * pz2 + qz * x8 * pz7 +
                  x1 * x8 * pz3 - x10 * pz0 - x6 * pz5 - x7 * pz1 - x9 * pz4;
    jacobian[8] = qx * qy * pz6 + qx * x3 * pz5 + qy * x8 * pz7 - x11 * pz2 -
                  x12 * pz1 - x13 * pz3 - x14 * pz0 + x3 * x8 * pz4;
  }

  const T x0 = jacobian[4] * jacobian[8];
  const T x1 = jacobian[5] * jacobian[7];
  const T x2 = jacobian[1] * jacobian[8];
  const T x3 = jacobian[1] * jacobian[5];
  const T x4 = jacobian[2] * jacobian[4];

  *adjugate0 = (AdjugateType)(x0 - x1);
  *adjugate1 = (AdjugateType)(jacobian[2] * jacobian[7] - x2);
  *adjugate2 = (AdjugateType)(x3 - x4);
  *adjugate3 =
      (AdjugateType)(-jacobian[3] * jacobian[8] + jacobian[5] * jacobian[6]);
  *adjugate4 =
      (AdjugateType)(jacobian[0] * jacobian[8] - jacobian[2] * jacobian[6]);
  *adjugate5 =
      (AdjugateType)(-jacobian[0] * jacobian[5] + jacobian[2] * jacobian[3]);
  *adjugate6 =
      (AdjugateType)(jacobian[3] * jacobian[7] - jacobian[4] * jacobian[6]);
  *adjugate7 =
      (AdjugateType)(-jacobian[0] * jacobian[7] + jacobian[1] * jacobian[6]);
  *adjugate8 =
      (AdjugateType)(jacobian[0] * jacobian[4] - jacobian[1] * jacobian[3]);
  *determinant =
      (DeterminantType)(jacobian[0] * x0 - jacobian[0] * x1 +
                        jacobian[2] * jacobian[3] * jacobian[7] -
                        jacobian[3] * x2 + jacobian[6] * x3 - jacobian[6] * x4);
}

template <typename T, typename AdjugateType>
static SMESH_INLINE SMESH_HOST_DEVICE void hex8_adjugate(
    // X
    const T px0, const T px1, const T px2, const T px3, const T px4,
    const T px5, const T px6, const T px7,
    // Y
    const T py0, const T py1, const T py2, const T py3, const T py4,
    const T py5, const T py6, const T py7,
    // Z
    const T pz0, const T pz1, const T pz2, const T pz3, const T pz4,
    const T pz5, const T pz6, const T pz7,
    // QPs
    const T qx, const T qy, const T qz,
    // Out
    AdjugateType *const SMESH_RESTRICT adjugate0,
    AdjugateType *const SMESH_RESTRICT adjugate1,
    AdjugateType *const SMESH_RESTRICT adjugate2,
    AdjugateType *const SMESH_RESTRICT adjugate3,
    AdjugateType *const SMESH_RESTRICT adjugate4,
    AdjugateType *const SMESH_RESTRICT adjugate5,
    AdjugateType *const SMESH_RESTRICT adjugate6,
    AdjugateType *const SMESH_RESTRICT adjugate7,
    AdjugateType *const SMESH_RESTRICT adjugate8) {
  T jacobian[9];
  {
    const T x0 = qy * qz;
    const T x1 = 1 - qz;
    const T x2 = qy * x1;
    const T x3 = 1 - qy;
    const T x4 = qz * x3;
    const T x5 = x1 * x3;
    const T x6 = qx * qz;
    const T x7 = qx * x1;
    const T x8 = 1 - qx;
    const T x9 = qz * x8;
    const T x10 = x1 * x8;
    const T x11 = qx * qy;
    const T x12 = qx * x3;
    const T x13 = qy * x8;
    const T x14 = x3 * x8;

    jacobian[0] = x0 * px6 - x0 * px7 + x2 * px2 - x2 * px3 - x4 * px4 +
                  x4 * px5 - x5 * px0 + x5 * px1;
    jacobian[1] = qx * qz * px6 + qx * x1 * px2 + qz * x8 * px7 +
                  x1 * x8 * px3 - x10 * px0 - x6 * px5 - x7 * px1 - x9 * px4;
    jacobian[2] = qx * qy * px6 + qx * x3 * px5 + qy * x8 * px7 - x11 * px2 -
                  x12 * px1 - x13 * px3 - x14 * px0 + x3 * x8 * px4;
    jacobian[3] = x0 * py6 - x0 * py7 + x2 * py2 - x2 * py3 - x4 * py4 +
                  x4 * py5 - x5 * py0 + x5 * py1;
    jacobian[4] = qx * qz * py6 + qx * x1 * py2 + qz * x8 * py7 +
                  x1 * x8 * py3 - x10 * py0 - x6 * py5 - x7 * py1 - x9 * py4;
    jacobian[5] = qx * qy * py6 + qx * x3 * py5 + qy * x8 * py7 - x11 * py2 -
                  x12 * py1 - x13 * py3 - x14 * py0 + x3 * x8 * py4;
    jacobian[6] = x0 * pz6 - x0 * pz7 + x2 * pz2 - x2 * pz3 - x4 * pz4 +
                  x4 * pz5 - x5 * pz0 + x5 * pz1;
    jacobian[7] = qx * qz * pz6 + qx * x1 * pz2 + qz * x8 * pz7 +
                  x1 * x8 * pz3 - x10 * pz0 - x6 * pz5 - x7 * pz1 - x9 * pz4;
    jacobian[8] = qx * qy * pz6 + qx * x3 * pz5 + qy * x8 * pz7 - x11 * pz2 -
                  x12 * pz1 - x13 * pz3 - x14 * pz0 + x3 * x8 * pz4;
  }

  const T x0 = jacobian[4] * jacobian[8];
  const T x1 = jacobian[5] * jacobian[7];
  const T x2 = jacobian[1] * jacobian[8];
  const T x3 = jacobian[1] * jacobian[5];
  const T x4 = jacobian[2] * jacobian[4];

  *adjugate0 = x0 - x1;
  *adjugate1 = jacobian[2] * jacobian[7] - x2;
  *adjugate2 = x3 - x4;
  *adjugate3 = -jacobian[3] * jacobian[8] + jacobian[5] * jacobian[6];
  *adjugate4 = jacobian[0] * jacobian[8] - jacobian[2] * jacobian[6];
  *adjugate5 = -jacobian[0] * jacobian[5] + jacobian[2] * jacobian[3];
  *adjugate6 = jacobian[3] * jacobian[7] - jacobian[4] * jacobian[6];
  *adjugate7 = -jacobian[0] * jacobian[7] + jacobian[1] * jacobian[6];
  *adjugate8 = jacobian[0] * jacobian[4] - jacobian[1] * jacobian[3];
}

template <typename T, typename OutType>
static SMESH_INLINE SMESH_HOST_DEVICE void hex8_l2_project(
    const T jacobian_determinant, const T qx, const T qy, const T qz,
    const T qw, const T v, OutType *const SMESH_RESTRICT out0,
    OutType *const SMESH_RESTRICT out1, OutType *const SMESH_RESTRICT out2,
    OutType *const SMESH_RESTRICT out3, OutType *const SMESH_RESTRICT out4,
    OutType *const SMESH_RESTRICT out5, OutType *const SMESH_RESTRICT out6,
    OutType *const SMESH_RESTRICT out7) {
  const T x0 = 1 - qx;
  const T x1 = 1 - qy;
  const T x2 = (qw * jacobian_determinant) * v;
  const T x3 = x2 * (1 - qz);
  const T x4 = x1 * x3;
  const T x5 = qy * x3;
  const T x6 = qz * x2;
  const T x7 = x1 * x6;
  const T x8 = qy * x6;
  *out0 += x0 * x4;
  *out1 += qx * x4;
  *out2 += qx * x5;
  *out3 += x0 * x5;
  *out4 += x0 * x7;
  *out5 += qx * x7;
  *out6 += qx * x8;
  *out7 += x0 * x8;
}

} // namespace smesh

#endif // SMESH_HEX8_INLINE_HPP
