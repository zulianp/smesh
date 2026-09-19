#ifndef SMESH_SSTET4_HPP
#define SMESH_SSTET4_HPP

#include "smesh_base.hpp"

namespace smesh {

/// Combinatorial tet / triangle counts. N_tet(n) = 0 for n < 0.
static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_n_tri(const int n) {
    if (n < 0) {
        return 0;
    }
    return (n + 1) * (n + 2) / 2;
}

static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_n_tet(const int n) {
    if (n < 0) {
        return 0;
    }
    return (n + 1) * (n + 2) * (n + 3) / 6;
}

/// Lattice nodes on a tet with L micro-edges: (L+1)(L+2)(L+3)/6.
static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_nxe(const int level) { return sstet4_n_tet(level); }

/// Kuhn / Freudenthal micro-tets per macro tet: L^3 (HyTeG / PROTEUS).
static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_txe(const int level) { return level * level * level; }

static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_nxedge(const int L) { return L > 1 ? (L - 1) : 0; }

static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_nxface(const int L) { return (L > 2) ? ((L - 1) * (L - 2) / 2) : 0; }

static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_nxvol(const int L) {
    return (L > 3) ? ((L - 1) * (L - 2) * (L - 3) / 6) : 0;
}

/// HyTeG vertex-DoF linearization on the integer tet lattice.
/// Barycentric (x,y,z) toward vertices 1,2,3 with x,y,z >= 0 and x+y+z <= L
/// (vertex 0 is (0,0,0)). z is slowest, then y, then x — tet analogue of
/// \c sshex8_lidx.
/// t_L(x,y,z) = N_tet(L) - N_tet(L-z) + N_tri(L-z) - N_tri(L-z-y) + x
static SMESH_INLINE SMESH_HOST_DEVICE int sstet4_lidx(const int L, const int x, const int y, const int z) {
    SMESH_ASSERT(x >= 0 && y >= 0 && z >= 0);
    SMESH_ASSERT(x + y + z <= L);
    const int ret = sstet4_n_tet(L) - sstet4_n_tet(L - z) + sstet4_n_tri(L - z) - sstet4_n_tri(L - z - y) + x;
    SMESH_ASSERT(ret >= 0);
    SMESH_ASSERT(ret < sstet4_nxe(L));
    return ret;
}

/// 2D barycentric index on a tet face lattice packed as \c sstet4_fill_side
/// (same as the z=0 slice of \c sstet4_lidx).
static SMESH_INLINE SMESH_HOST_DEVICE int sstri_lidx(const int L, const int x, const int y) {
    SMESH_ASSERT(x >= 0 && y >= 0);
    SMESH_ASSERT(x + y <= L);
    const int ret = sstet4_n_tri(L) - sstet4_n_tri(L - y) + x;
    SMESH_ASSERT(ret >= 0);
    SMESH_ASSERT(ret < sstet4_n_tri(L));
    return ret;
}

static SMESH_INLINE SMESH_HOST_DEVICE int sstri_nxe(const int L) { return sstet4_n_tri(L); }

static SMESH_INLINE SMESH_HOST_DEVICE int sstri_txe(const int L) { return L * L; }

}  // namespace smesh

#endif  // SMESH_SSTET4_HPP
