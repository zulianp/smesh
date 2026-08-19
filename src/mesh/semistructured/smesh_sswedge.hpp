#ifndef SMESH_SSWEDGE_HPP
#define SMESH_SSWEDGE_HPP

#include "smesh_base.hpp"
#include "smesh_sstet4.hpp"

namespace smesh {

/// Lattice nodes on a prism: triangle × 1D. N = N_tri(L) * (L+1).
static SMESH_INLINE SMESH_HOST_DEVICE int sswedge_nxe(const int L) {
    return sstet4_n_tri(L) * (L + 1);
}

static SMESH_INLINE SMESH_HOST_DEVICE int sswedge_nxedge(const int L) { return L > 1 ? (L - 1) : 0; }

static SMESH_INLINE SMESH_HOST_DEVICE int sswedge_nx_tri_face(const int L) { return sstet4_nxface(L); }

static SMESH_INLINE SMESH_HOST_DEVICE int sswedge_nx_quad_face(const int L) {
    return L > 1 ? ((L - 1) * (L - 1)) : 0;
}

static SMESH_INLINE SMESH_HOST_DEVICE int sswedge_nxvol(const int L) {
    return (L > 2) ? (sstet4_nxface(L) * (L - 1)) : 0;
}

/// (x,y,z) with x+y <= L, 0 <= z <= L. z is slowest; (x,y) is the HyTeG triangle.
static SMESH_INLINE SMESH_HOST_DEVICE int sswedge_lidx(const int L, const int x, const int y, const int z) {
    SMESH_ASSERT(x >= 0 && y >= 0 && z >= 0);
    SMESH_ASSERT(x + y <= L);
    SMESH_ASSERT(z <= L);
    const int ret = z * sstet4_n_tri(L) + sstet4_n_tri(L) - sstet4_n_tri(L - y) + x;
    SMESH_ASSERT(ret >= 0);
    SMESH_ASSERT(ret < sswedge_nxe(L));
    return ret;
}

}  // namespace smesh

#endif  // SMESH_SSWEDGE_HPP
