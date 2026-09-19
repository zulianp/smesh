#ifndef SMESH_SSPYRAMID_HPP
#define SMESH_SSPYRAMID_HPP

#include "smesh_base.hpp"
#include "smesh_sstet4.hpp"

namespace smesh {

/// sum_{i=1}^n i^2. Returns 0 for n < 1.
static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_sum_sq(const int n) {
    if (n < 1) {
        return 0;
    }
    return n * (n + 1) * (2 * n + 1) / 6;
}

/// Lattice nodes: shrinking quad. N = (L+1)(L+2)(2L+3)/6.
static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_nxe(const int L) {
    return (L + 1) * (L + 2) * (2 * L + 3) / 6;
}

static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_nxedge(const int L) { return L > 1 ? (L - 1) : 0; }

static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_nx_tri_face(const int L) { return sstet4_nxface(L); }

static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_nx_quad_face(const int L) {
    return L > 1 ? ((L - 1) * (L - 1)) : 0;
}

static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_nxvol(const int L) {
    return (L > 2) ? sspyramid_sum_sq(L - 2) : 0;
}

/// Prefix nodes before layer k: sum_{t=0}^{k-1} (L-t+1)^2 = S(L+1) - S(L-k+1).
static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_prefix(const int L, const int k) {
    return sspyramid_sum_sq(L + 1) - sspyramid_sum_sq(L - k + 1);
}

/// (i,j,k) with 0 <= k <= L and 0 <= i,j <= L-k. k is slowest.
static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_lidx(const int L, const int i, const int j, const int k) {
    SMESH_ASSERT(k >= 0 && k <= L);
    SMESH_ASSERT(i >= 0 && j >= 0);
    SMESH_ASSERT(i <= L - k && j <= L - k);
    const int stride = L - k + 1;
    const int ret    = sspyramid_prefix(L, k) + j * stride + i;
    SMESH_ASSERT(ret >= 0);
    SMESH_ASSERT(ret < sspyramid_nxe(L));
    return ret;
}

/// Child pyramids per macro-pyramid (shrinking-quad layer stencil).
/// = sum_{k=0}^{L-1}[(L-k)^2 + (L-k-1)^2]  = L^2 + L(L-1)(2L-1)/3
static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_n_pyr(const int L) {
    if (L < 1) return 0;
    // L=1 → 1 (identity); L=2 → 6; L=4 → 44
    return L * L + L * (L - 1) * (2 * L - 1) / 3;
}

/// Child tets per macro-pyramid (seam tets between upward pyramids).
/// = sum_{k=0}^{L-1} 2*(L-k)*(L-k-1) = 2L(L^2-1)/3
static SMESH_INLINE SMESH_HOST_DEVICE int sspyramid_n_tet(const int L) {
    if (L < 2) return 0;
    return 2 * L * (L * L - 1) / 3;
}

}  // namespace smesh

#endif  // SMESH_SSPYRAMID_HPP
