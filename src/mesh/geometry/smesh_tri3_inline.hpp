#ifndef SMESH_TRI3_INLINE_HPP
#define SMESH_TRI3_INLINE_HPP

#include "smesh_base.hpp"

#ifndef POW2
#define POW2(a) ((a) * (a))
#endif

namespace smesh {

template <typename T, typename AdjugateType, typename DeterminantType>
static SMESH_INLINE SMESH_HOST_DEVICE void
tri3_adjugate_and_det(const T px0, const T px1, const T px2, const T py0, const T py1, const T py2,
                      AdjugateType *const SMESH_RESTRICT adjugate0, AdjugateType *const SMESH_RESTRICT adjugate1,
                      AdjugateType *const SMESH_RESTRICT adjugate2, AdjugateType *const SMESH_RESTRICT adjugate3,
                      DeterminantType *const SMESH_RESTRICT determinant) {
    const T a = -px0 + px1;
    const T b = -px0 + px2;
    const T c = -py0 + py1;
    const T d = -py0 + py2;
    *adjugate0   = (AdjugateType)d;
    *adjugate1   = (AdjugateType)(-b);
    *adjugate2   = (AdjugateType)(-c);
    *adjugate3   = (AdjugateType)a;
    *determinant = (DeterminantType)(a * d - b * c);
}

template <typename T, typename FFFType>
static SMESH_INLINE SMESH_HOST_DEVICE void tri3_fff(const T px0, const T px1, const T px2, const T py0, const T py1,
                                                    const T py2, FFFType *const SMESH_RESTRICT fff0,
                                                    FFFType *const SMESH_RESTRICT fff1, FFFType *const SMESH_RESTRICT fff2) {
    T adj0, adj1, adj2, adj3, det;
    tri3_adjugate_and_det(px0, px1, px2, py0, py1, py2, &adj0, &adj1, &adj2, &adj3, &det);
    const T inv = T(1) / (T(2) * det);
    *fff0       = (FFFType)((POW2(adj0) + POW2(adj1)) * inv);
    *fff1       = (FFFType)((adj0 * adj2 + adj1 * adj3) * inv);
    *fff2       = (FFFType)((POW2(adj2) + POW2(adj3)) * inv);
}

}  // namespace smesh

#endif  // SMESH_TRI3_INLINE_HPP
