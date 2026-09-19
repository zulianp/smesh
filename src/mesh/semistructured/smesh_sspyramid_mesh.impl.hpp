#ifndef SMESH_SSPYRAMID_MESH_IMPL_HPP
#define SMESH_SSPYRAMID_MESH_IMPL_HPP

#include "smesh_sspyramid.hpp"
#include "smesh_sspyramid_mesh.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
int sspyramid_fill_points(const int                                                level,
                          const ptrdiff_t                                          nelements,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                          const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
                          geom_t *SMESH_RESTRICT *const SMESH_RESTRICT             points) {
    const int corners[5] = {sspyramid_lidx(level, 0, 0, 0),
                            sspyramid_lidx(level, level, 0, 0),
                            sspyramid_lidx(level, level, level, 0),
                            sspyramid_lidx(level, 0, level, 0),
                            sspyramid_lidx(level, 0, 0, level)};

    for (int k = 0; k <= level; ++k) {
        const geom_t zeta = geom_t(k) / geom_t(level);
        const int    span = level - k;
        for (int j = 0; j <= span; ++j) {
            for (int i = 0; i <= span; ++i) {
                geom_t xi  = geom_t(0);
                geom_t eta = geom_t(0);
                if (span > 0) {
                    xi  = geom_t(i) / geom_t(span);
                    eta = geom_t(j) / geom_t(span);
                }
                const geom_t n0   = (geom_t(1) - xi) * (geom_t(1) - eta);
                const geom_t n1   = xi * (geom_t(1) - eta);
                const geom_t n2   = xi * eta;
                const geom_t n3   = (geom_t(1) - xi) * eta;
                const int    lidx = sspyramid_lidx(level, i, j, k);
                for (int d = 0; d < 3; d++) {
                    for (ptrdiff_t e = 0; e < nelements; e++) {
                        const geom_t base = n0 * macro_mesh_points[d][elements[corners[0]][e]] +
                                            n1 * macro_mesh_points[d][elements[corners[1]][e]] +
                                            n2 * macro_mesh_points[d][elements[corners[2]][e]] +
                                            n3 * macro_mesh_points[d][elements[corners[3]][e]];
                        const geom_t apex = macro_mesh_points[d][elements[corners[4]][e]];
                        points[d][elements[lidx][e]] = (geom_t(1) - zeta) * base + zeta * apex;
                    }
                }
            }
        }
    }
    return SMESH_SUCCESS;
}

/// Explode a single macro-pyramid SS block into PYRAMID5 + TET4 children.
///
/// Lattice: sspyramid_lidx(L, i, j, k), 0 <= k <= L, 0 <= i,j <= L-k.
///
/// For each layer k (s = L-k rows per side):
///   Upward pyramids (s² total):
///     base: (i,j,k),(i+1,j,k),(i+1,j+1,k),(i,j+1,k), apex: (i,j,k+1)
///     i,j in 0..s-1
///   Inverted pyramids ((s-1)² total, only if s >= 2):
///     base: (i,j,k+1),(i+1,j,k+1),(i+1,j+1,k+1),(i,j+1,k+1), apex: (i+1,j+1,k)
///     i,j in 0..s-2
///   Seam tets between adjacent upward pyramids:
///     horizontal (j boundary shared): (i+1,j,k),(i+1,j,k+1),(i+1,j+1,k),(i+1,j+1,k+1)
///       → tet: (i+1,j,k) (i+1,j,k+1) (i,j,k+1) (i+1,j+1,k)   [NO – use apex overlap]
///     The seam tets share the apex of two adjacent upward pyramids and a base edge.
///     Horizontal neighbor (i, j) vs (i+1, j): shared base edge (i+1,j,k)-(i+1,j+1,k)
///       apexes at (i,j,k+1) and (i+1,j,k+1) → tet: (i+1,j,k),(i+1,j+1,k),(i,j,k+1),(i+1,j,k+1)
///     Vertical neighbor (i, j) vs (i, j+1): shared base edge (i,j+1,k)-(i+1,j+1,k)
///       apexes at (i,j,k+1) and (i,j+1,k+1) → tet: (i,j+1,k),(i+1,j+1,k),(i,j,k+1),(i,j+1,k+1)
///
/// Winding convention: PYRAMID5 base is listed CCW viewed from below (away from apex);
///   for upward pyramids apex is above (+k), so base is at layer k (bottom).
///   VTK PYRAMID5 = {n0,n1,n2,n3, apex} where base is CCW from apex viewpoint.
///
/// Total pyramids per macro = sspyramid_n_pyr(L);  total tets = sspyramid_n_tet(L).
template <typename idx_t>
int sspyramid_to_pyramid5_and_tet4(const int                                               level,
                                   const ptrdiff_t                                         nelements,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                   idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             pyr5_elements,
                                   idx_t *SMESH_RESTRICT *const SMESH_RESTRICT             tet4_elements) {
    const int L        = level;
    const int n_pyr    = sspyramid_n_pyr(L);
    const int n_tet    = sspyramid_n_tet(L);

    // Macro-element loop (SoA: column = macro element, row = local node index).
    // We fill child index as:  child_idx = macro_e * n_pyr + local_child
    // and similarly for tets.
    ptrdiff_t lp = 0;  // local pyramid index within this macro
    ptrdiff_t lt = 0;  // local tet index within this macro

    for (int k = 0; k < L; ++k) {
        const int s = L - k;  // layer size (s x s grid at level k)

        // --- Upward pyramids: s*s per layer ---
        for (int j = 0; j < s; ++j) {
            for (int i = 0; i < s; ++i) {
                const int b0 = sspyramid_lidx(L, i,     j,     k);
                const int b1 = sspyramid_lidx(L, i + 1, j,     k);
                const int b2 = sspyramid_lidx(L, i + 1, j + 1, k);
                const int b3 = sspyramid_lidx(L, i,     j + 1, k);
                const int ap = sspyramid_lidx(L, i,     j,     k + 1);

                SMESH_ASSERT(lp < n_pyr);
                for (ptrdiff_t e = 0; e < nelements; ++e) {
                    // VTK PYRAMID5: base CCW when viewed from below (i.e. from -k side)
                    // apex is at k+1 (above), so base CCW from apex = CCW when looking down.
                    pyr5_elements[0][e * n_pyr + lp] = elements[b0][e];
                    pyr5_elements[1][e * n_pyr + lp] = elements[b1][e];
                    pyr5_elements[2][e * n_pyr + lp] = elements[b2][e];
                    pyr5_elements[3][e * n_pyr + lp] = elements[b3][e];
                    pyr5_elements[4][e * n_pyr + lp] = elements[ap][e];
                }
                ++lp;
            }
        }

        // --- Inverted pyramids: (s-1)*(s-1) per layer (only if s >= 2) ---
        if (s >= 2) {
            for (int j = 0; j < s - 1; ++j) {
                for (int i = 0; i < s - 1; ++i) {
                    // base at upper layer k+1, apex at lower layer k
                    const int u0 = sspyramid_lidx(L, i,     j,     k + 1);
                    const int u1 = sspyramid_lidx(L, i + 1, j,     k + 1);
                    const int u2 = sspyramid_lidx(L, i + 1, j + 1, k + 1);
                    const int u3 = sspyramid_lidx(L, i,     j + 1, k + 1);
                    // apex of inverted pyramid: one step diagonally into lower layer
                    const int ap = sspyramid_lidx(L, i + 1, j + 1, k);

                    SMESH_ASSERT(lp < n_pyr);
                    for (ptrdiff_t e = 0; e < nelements; ++e) {
                        // Inverted: base at k+1 (above the apex), apex below.
                        // For consistent outward normals, reverse base winding:
                        pyr5_elements[0][e * n_pyr + lp] = elements[u0][e];
                        pyr5_elements[1][e * n_pyr + lp] = elements[u3][e];
                        pyr5_elements[2][e * n_pyr + lp] = elements[u2][e];
                        pyr5_elements[3][e * n_pyr + lp] = elements[u1][e];
                        pyr5_elements[4][e * n_pyr + lp] = elements[ap][e];
                    }
                    ++lp;
                }
            }
        }

        // --- Seam tets: 2*s*(s-1) per layer ---
        // Horizontal seams: adjacent i-neighbors share the edge (i+1,j,k)-(i+1,j+1,k)
        //   upward-pyr(i,j) apex=(i,j,k+1), upward-pyr(i+1,j) apex=(i+1,j,k+1)
        for (int j = 0; j < s; ++j) {
            for (int i = 0; i < s - 1; ++i) {
                const int a0 = sspyramid_lidx(L, i,     j,     k);      // shared base edge v0
                const int a1 = sspyramid_lidx(L, i + 1, j,     k);      // shared base edge v1 (between pyramids)
                const int a2 = sspyramid_lidx(L, i + 1, j + 1, k);      // shared base edge v2 (between pyramids)
                const int ap0 = sspyramid_lidx(L, i,     j,     k + 1); // apex of left pyramid
                const int ap1 = sspyramid_lidx(L, i + 1, j,     k + 1); // apex of right pyramid

                // Tet between the two upward pyramids: uses the shared base edge + both apexes
                // tet: a1, a2, ap0, ap1 (positive volume winding)
                SMESH_ASSERT(lt < n_tet);
                (void)a0;
                for (ptrdiff_t e = 0; e < nelements; ++e) {
                    tet4_elements[0][e * n_tet + lt] = elements[a1][e];
                    tet4_elements[1][e * n_tet + lt] = elements[a2][e];
                    tet4_elements[2][e * n_tet + lt] = elements[ap0][e];
                    tet4_elements[3][e * n_tet + lt] = elements[ap1][e];
                }
                ++lt;
            }
        }

        // Vertical seams: adjacent j-neighbors share the edge (i,j+1,k)-(i+1,j+1,k)
        //   upward-pyr(i,j) apex=(i,j,k+1), upward-pyr(i,j+1) apex=(i,j+1,k+1)
        for (int j = 0; j < s - 1; ++j) {
            for (int i = 0; i < s; ++i) {
                const int b0 = sspyramid_lidx(L, i,     j + 1, k);      // shared base edge v0
                const int b1 = sspyramid_lidx(L, i + 1, j + 1, k);      // shared base edge v1
                const int bp0 = sspyramid_lidx(L, i,     j,     k + 1); // apex of front pyramid
                const int bp1 = sspyramid_lidx(L, i,     j + 1, k + 1); // apex of back pyramid

                SMESH_ASSERT(lt < n_tet);
                for (ptrdiff_t e = 0; e < nelements; ++e) {
                    tet4_elements[0][e * n_tet + lt] = elements[b0][e];
                    tet4_elements[1][e * n_tet + lt] = elements[b1][e];
                    tet4_elements[2][e * n_tet + lt] = elements[bp0][e];
                    tet4_elements[3][e * n_tet + lt] = elements[bp1][e];
                }
                ++lt;
            }
        }
    }

    SMESH_ASSERT(lp == n_pyr);
    SMESH_ASSERT(lt == n_tet);
    return SMESH_SUCCESS;
}

}  // namespace smesh

#endif  // SMESH_SSPYRAMID_MESH_IMPL_HPP
