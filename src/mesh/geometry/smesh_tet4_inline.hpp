#ifndef SMESH_TET4_INLINE_HPP
#define SMESH_TET4_INLINE_HPP

#include "smesh_base.hpp"

#ifndef POW2
#define POW2(a) ((a) * (a))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

namespace smesh {

template <typename T, typename FFFType>
static SMESH_INLINE SMESH_HOST_DEVICE void
tet4_fff(const T px0, const T px1, const T px2, const T px3, const T py0,
         const T py1, const T py2, const T py3, const T pz0, const T pz1,
         const T pz2, const T pz3, FFFType *const SMESH_RESTRICT fff0,
         FFFType *const SMESH_RESTRICT fff1, FFFType *const SMESH_RESTRICT fff2,
         FFFType *const SMESH_RESTRICT fff3, FFFType *const SMESH_RESTRICT fff4,
         FFFType *const SMESH_RESTRICT fff5) {
  const T x0 = -px0 + px1;
  const T x1 = -py0 + py2;
  const T x2 = -pz0 + pz3;
  const T x3 = x1 * x2;
  const T x4 = x0 * x3;
  const T x5 = -py0 + py3;
  const T x6 = -pz0 + pz2;
  const T x7 = x5 * x6;
  const T x8 = x0 * x7;
  const T x9 = -py0 + py1;
  const T x10 = -px0 + px2;
  const T x11 = x10 * x2;
  const T x12 = x11 * x9;
  const T x13 = -pz0 + pz1;
  const T x14 = x10 * x5;
  const T x15 = x13 * x14;
  const T x16 = -px0 + px3;
  const T x17 = x16 * x6 * x9;
  const T x18 = x1 * x16;
  const T x19 = x13 * x18;
  const T x20 = -T(1.0 / 6.0) * x12 + T(1.0 / 6.0) * x15 + T(1.0 / 6.0) * x17 -
                T(1.0 / 6.0) * x19 + T(1.0 / 6.0) * x4 - T(1.0 / 6.0) * x8;
  const T x21 = x14 - x18;
  const T x22 = T(1.0) / POW2(-x12 + x15 + x17 - x19 + x4 - x8);
  const T x23 = -x11 + x16 * x6;
  const T x24 = x3 - x7;
  const T x25 = -x0 * x5 + x16 * x9;
  const T x26 = x21 * x22;
  const T x27 = x0 * x2 - x13 * x16;
  const T x28 = x22 * x23;
  const T x29 = x13 * x5 - x2 * x9;
  const T x30 = x22 * x24;
  const T x31 = x0 * x1 - x10 * x9;
  const T x32 = -x0 * x6 + x10 * x13;
  const T x33 = -x1 * x13 + x6 * x9;
  *fff0 = x20 * (POW2(x21) * x22 + x22 * POW2(x23) + x22 * POW2(x24));
  *fff1 = x20 * (x25 * x26 + x27 * x28 + x29 * x30);
  *fff2 = x20 * (x26 * x31 + x28 * x32 + x30 * x33);
  *fff3 = x20 * (x22 * POW2(x25) + x22 * POW2(x27) + x22 * POW2(x29));
  *fff4 = x20 * (x22 * x25 * x31 + x22 * x27 * x32 + x22 * x29 * x33);
  *fff5 = x20 * (x22 * POW2(x31) + x22 * POW2(x32) + x22 * POW2(x33));
}

template <typename T>
static SMESH_INLINE SMESH_HOST_DEVICE T tet4_det_fff(const T fff0, const T fff1,
                                                     const T fff2, const T fff3,
                                                     const T fff4,
                                                     const T fff5) {
  return fff0 * fff3 * fff5 - fff0 * POW2(fff4) - POW2(fff1) * fff5 +
         2 * fff1 * fff2 * fff4 - POW2(fff2) * fff3;
}

template <typename T, typename AdjugateType, typename DeterminantType>
static SMESH_INLINE SMESH_HOST_DEVICE void
tet4_adjugate_and_det(const T px0, const T px1, const T px2, const T px3,
                      const T py0, const T py1, const T py2, const T py3,
                      const T pz0, const T pz1, const T pz2, const T pz3,
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
  // Compute jacobian in high precision
  T jacobian[9];
  jacobian[0] = -px0 + px1;
  jacobian[1] = -px0 + px2;
  jacobian[2] = -px0 + px3;
  jacobian[3] = -py0 + py1;
  jacobian[4] = -py0 + py2;
  jacobian[5] = -py0 + py3;
  jacobian[6] = -pz0 + pz1;
  jacobian[7] = -pz0 + pz2;
  jacobian[8] = -pz0 + pz3;

  const T x0 = jacobian[4] * jacobian[8];
  const T x1 = jacobian[5] * jacobian[7];
  const T x2 = jacobian[1] * jacobian[8];
  const T x3 = jacobian[1] * jacobian[5];
  const T x4 = jacobian[2] * jacobian[4];

  // Store adjugate in lower precision
  *adjugate0 = x0 - x1;
  *adjugate1 = jacobian[2] * jacobian[7] - x2;
  *adjugate2 = x3 - x4;
  *adjugate3 = -jacobian[3] * jacobian[8] + jacobian[5] * jacobian[6];
  *adjugate4 = jacobian[0] * jacobian[8] - jacobian[2] * jacobian[6];
  *adjugate5 = -jacobian[0] * jacobian[5] + jacobian[2] * jacobian[3];
  *adjugate6 = jacobian[3] * jacobian[7] - jacobian[4] * jacobian[6];
  *adjugate7 = -jacobian[0] * jacobian[7] + jacobian[1] * jacobian[6];
  *adjugate8 = jacobian[0] * jacobian[4] - jacobian[1] * jacobian[3];

  // Store determinant in lower precision
  *determinant = jacobian[0] * x0 - jacobian[0] * x1 +
                 jacobian[2] * jacobian[3] * jacobian[7] - jacobian[3] * x2 +
                 jacobian[6] * x3 - jacobian[6] * x4;
}

} // namespace smesh

#endif // SMESH_TET4_INLINE_HPP
