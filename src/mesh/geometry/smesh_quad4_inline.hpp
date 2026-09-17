#ifndef SMESH_QUAD4_INLINE_HPP
#define SMESH_QUAD4_INLINE_HPP

#include "smesh_base.hpp"

#ifndef POW2
#define POW2(a) ((a) * (a))
#endif

namespace smesh {

template <typename T, typename AdjugateType, typename DeterminantType>
static SMESH_INLINE SMESH_HOST_DEVICE void
quad4_adjugate_and_det(const T px0, const T px1, const T px2, const T px3, const T py0, const T py1, const T py2,
                       const T py3, const T qx, const T qy, AdjugateType *const SMESH_RESTRICT adjugate0,
                       AdjugateType *const SMESH_RESTRICT adjugate1, AdjugateType *const SMESH_RESTRICT adjugate2,
                       AdjugateType *const SMESH_RESTRICT adjugate3, DeterminantType *const SMESH_RESTRICT determinant) {
    const T oqx = T(1) - qx;
    const T oqy = T(1) - qy;
    const T a   = -oqy * px0 + oqy * px1 + qy * px2 - qy * px3;
    const T b   = -oqx * px0 - qx * px1 + qx * px2 + oqx * px3;
    const T c   = -oqy * py0 + oqy * py1 + qy * py2 - qy * py3;
    const T d   = -oqx * py0 - qx * py1 + qx * py2 + oqx * py3;
    *adjugate0   = (AdjugateType)d;
    *adjugate1   = (AdjugateType)(-b);
    *adjugate2   = (AdjugateType)(-c);
    *adjugate3   = (AdjugateType)a;
    *determinant = (DeterminantType)(a * d - b * c);
}

template <typename T, typename FFFType>
static SMESH_INLINE SMESH_HOST_DEVICE void
quad4_fff(const T px0, const T px1, const T px2, const T px3, const T py0, const T py1, const T py2, const T py3,
          const T qx, const T qy, FFFType *const SMESH_RESTRICT fff0, FFFType *const SMESH_RESTRICT fff1,
          FFFType *const SMESH_RESTRICT fff2) {
    T adj0, adj1, adj2, adj3, det;
    quad4_adjugate_and_det(px0, px1, px2, px3, py0, py1, py2, py3, qx, qy, &adj0, &adj1, &adj2, &adj3, &det);
    const T inv = T(1) / det;
    *fff0       = (FFFType)((POW2(adj0) + POW2(adj1)) * inv);
    *fff1       = (FFFType)((adj0 * adj2 + adj1 * adj3) * inv);
    *fff2       = (FFFType)((POW2(adj2) + POW2(adj3)) * inv);
}

}  // namespace smesh

#endif  // SMESH_QUAD4_INLINE_HPP
