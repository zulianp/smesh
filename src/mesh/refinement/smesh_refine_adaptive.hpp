#ifndef SMESH_REFINE_ADAPTIVE_HPP
#define SMESH_REFINE_ADAPTIVE_HPP

#include "smesh_base.hpp"
#include "smesh_refine.hpp"

namespace smesh {

/// Red-green action for a 6-bit edge mask. `which` is the Exodus edge (0..5)
/// or the local face (0..3). Any other mask is promoted to a Bey 8-split.
///
/// A red child's forest slot is its Bey index, 0..7. Green children are stored
/// in the forest too, so the next red level bisects the transition edge instead
/// of 8-splitting the neighbor. Green-edge slot = 8 + edge * 2 + child.
/// Green-face slot = 20 + face * 4 + child.
enum {
    TET4_ADAPT_KEEP    = 0,
    TET4_ADAPT_EDGE    = 1,
    TET4_ADAPT_FACE    = 2,
    TET4_ADAPT_PROMOTE = 3
};

enum {
    TET4_ADAPT_SLOT_GREEN_EDGE = 8,
    TET4_ADAPT_SLOT_GREEN_FACE = 20,
    /// Irregular closure: slot = this + parent face. Boundary is local face 3.
    TET4_ADAPT_SLOT_FAN = 36
};

struct Tet4AdaptAction {
    u8 kind;
    u8 which;
};

/// Indexed by the bitmask of bisected Exodus edges (`tet4_refine_edges`).
inline constexpr Tet4AdaptAction tet4_adapt_action[64] = {
        {TET4_ADAPT_KEEP, 0},    {TET4_ADAPT_EDGE, 0},    {TET4_ADAPT_EDGE, 1},    {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_EDGE, 2},    {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_FACE, 3},
        {TET4_ADAPT_EDGE, 3},    {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_EDGE, 4},    {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_FACE, 0},    {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_EDGE, 5},    {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_FACE, 2},    {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_FACE, 1},    {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0},
        {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}, {TET4_ADAPT_PROMOTE, 0}};

/// Parent face of red child `c`'s local face, or -1 when that face is interior.
/// Child 7's face on parent face 2 is local face 3 (`tet4_face_child_lfi`).
inline constexpr int tet4_child_face_parent[8][4] = {
        {0, -1, 2, 3}, {0, 1, -1, 3}, {-1, 1, 2, 3}, {0, 1, 2, -1},
        {-1, -1, -1, 3}, {0, -1, -1, -1}, {-1, 1, -1, -1}, {-1, -1, -1, 2}};

/// Parent edge of red child `c`'s local edge, or -1.
inline constexpr int tet4_child_edge_parent[8][6] = {
        {0, -1, 2, 3, -1, -1}, {0, 1, -1, -1, 4, -1}, {-1, 1, 2, -1, -1, 5}, {-1, -1, -1, 3, 4, 5},
        {-1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1}};

/// Two tets for one bisected edge. Node ids are corners 0..3 and Exodus mids 4..9.
/// Orders keep the parent's signed volume.
inline constexpr int tet4_green_edge[6][2][4] = {
        {{0, 4, 2, 3}, {4, 1, 2, 3}}, {{1, 5, 0, 3}, {5, 2, 0, 3}}, {{0, 6, 3, 1}, {6, 2, 3, 1}},
        {{0, 7, 1, 2}, {7, 3, 1, 2}}, {{1, 8, 2, 0}, {8, 3, 2, 0}}, {{2, 9, 0, 1}, {9, 3, 0, 1}}};

/// Four tets for one fully bisected face. Triangles match `tet4_face_child`.
inline constexpr int tet4_green_face[4][4][4] = {
        {{0, 7, 4, 2}, {4, 8, 1, 2}, {7, 3, 8, 2}, {7, 8, 4, 2}},
        {{1, 8, 5, 0}, {5, 9, 2, 0}, {8, 3, 9, 0}, {5, 8, 9, 0}},
        {{0, 6, 7, 1}, {6, 2, 9, 1}, {7, 9, 3, 1}, {7, 6, 9, 1}},
        {{0, 4, 6, 3}, {4, 1, 5, 3}, {6, 5, 2, 3}, {4, 5, 6, 3}}};

/// Green-edge child local face → leaf face, or -1.
inline constexpr int tet4_green_edge_face[6][2][4] = {
        {{0, -1, 2, 3}, {0, 1, -1, 3}}, {{1, -1, 0, 3}, {1, 2, -1, 3}}, {{3, -1, 0, 2}, {3, 1, -1, 2}},
        {{2, -1, 3, 0}, {2, 1, -1, 0}}, {{0, -1, 3, 1}, {0, 2, -1, 1}}, {{1, -1, 3, 2}, {1, 0, -1, 2}}};

/// Green-edge child local edge → leaf edge, or -1.
inline constexpr int tet4_green_edge_edge[6][2][6] = {
        {{0, -1, 2, 3, -1, 5}, {0, 1, -1, -1, 4, 5}}, {{1, -1, 0, 4, -1, 3}, {1, 2, -1, -1, 5, 3}},
        {{2, -1, 3, 0, -1, 4}, {2, 5, -1, -1, 1, 4}}, {{3, -1, 0, 2, -1, 1}, {3, 4, -1, -1, 5, 1}},
        {{4, -1, 1, 0, -1, 2}, {4, 5, -1, -1, 3, 2}}, {{5, -1, 2, 1, -1, 0}, {5, 3, -1, -1, 4, 0}}};

/// Green-face child local face → leaf face, or -1.
inline constexpr int tet4_green_face_face[4][4][4] = {
        {{2, -1, 3, 0}, {-1, 1, 3, 0}, {2, 1, -1, 0}, {-1, -1, -1, 0}},
        {{0, -1, 3, 1}, {-1, 2, 3, 1}, {0, 2, -1, 1}, {-1, -1, -1, 1}},
        {{3, -1, 0, 2}, {3, 1, -1, 2}, {-1, 1, 0, 2}, {-1, -1, -1, 2}},
        {{0, -1, 2, 3}, {0, 1, -1, 3}, {-1, 1, 2, 3}, {-1, -1, -1, 3}}};

/// Green-face child local edge → leaf edge, or -1.
inline constexpr int tet4_green_face_edge[4][4][6] = {
        {{3, -1, 0, 2, -1, -1}, {-1, 4, 0, -1, -1, 1}, {3, 4, -1, -1, 5, -1}, {-1, -1, -1, -1, -1, -1}},
        {{4, -1, 1, 0, -1, -1}, {-1, 5, 1, -1, -1, 2}, {4, 5, -1, -1, 3, -1}, {-1, -1, -1, -1, -1, -1}},
        {{2, -1, 3, 0, -1, -1}, {2, 5, -1, -1, 1, -1}, {-1, 5, 3, -1, -1, 4}, {-1, -1, -1, -1, -1, -1}},
        {{0, -1, 2, 3, -1, -1}, {0, 1, -1, -1, 4, -1}, {-1, 1, 2, -1, -1, 5}, {-1, -1, -1, -1, -1, -1}}};

/// Serial TET4 red-green refinement. `level` is the target red depth of every
/// marked element when \p mark_level is null. A non-null \p mark_level gives
/// each mark its own depth; \p level is then the maximum and sizes the closure
/// bound. An element whose target is L raises every face neighbor to at least
/// L-1, and that rule repeats, so the halo is L-1 elements wide. Neighbors are
/// not raised to L. The kernel allocates every output buffer with `SMESH_ALLOC`.
/// Face CRS has `4 * n_elements` rows (coarse element-major, then local face).
/// Edge CRS has `6 * n_elements` rows. Midpoint triples are in creation order.
/// When \p out_root is non-null it receives the coarse element index of each output tet.
template <typename idx_t, typename count_t, typename geom_t>
int mesh_refine_adaptive(
        const ptrdiff_t n_elements, const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
        const int spatial_dim, const ptrdiff_t n_nodes, const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
        const idx_t *const SMESH_RESTRICT marked, const ptrdiff_t n_marked, const int level,
        idx_t ***out_elements, ptrdiff_t *n_out_elements, geom_t ***out_points, ptrdiff_t *n_out_nodes,
        count_t **face_ptr, idx_t **face_elem, idx_t **face_lfi, count_t **edge_ptr, idx_t **edge_elem,
        idx_t **edge_lei, idx_t **mid_lo, idx_t **mid_hi, idx_t **mid_id, ptrdiff_t *n_mids, idx_t **out_root,
        const int *mark_level);

}  // namespace smesh

#endif  // SMESH_REFINE_ADAPTIVE_HPP
