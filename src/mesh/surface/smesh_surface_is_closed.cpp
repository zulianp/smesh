#include "smesh_surface_is_closed.hpp"

#include "smesh_elem_type.hpp"

namespace smesh {

namespace {

static constexpr int tri3_side_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
static constexpr int quad4_side_nodes[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

template <int NSIDES>
SMESH_FORCE_INLINE bool has_side(
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const element_idx_t e, const idx_t a, const idx_t b,
    const int (&side_nodes)[NSIDES][2]) {
  for (int s = 0; s < NSIDES; ++s) {
    const idx_t c = elems[side_nodes[s][0]][e];
    const idx_t d = elems[side_nodes[s][1]][e];
    if ((c == a && d == b) || (c == b && d == a)) {
      return true;
    }
  }

  return false;
}

template <int NSIDES>
bool surface_is_closed_edges(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx,
    const int (&side_nodes)[NSIDES][2]) {
  bool closed = true;

#pragma omp parallel for reduction(&& : closed)
  for (ptrdiff_t e = 0; e < n_elements; ++e) {
    bool element_closed = true;

    for (int s = 0; s < NSIDES; ++s) {
      const idx_t a = elems[side_nodes[s][0]][e];
      const idx_t b = elems[side_nodes[s][1]][e];
      const count_t begin = n2e_ptr[a];
      const count_t end = n2e_ptr[a + 1];
      bool found = false;

      for (count_t it = begin; it < end; ++it) {
        const element_idx_t other = n2e_idx[it];
        if (other != static_cast<element_idx_t>(e) &&
            has_side<NSIDES>(elems, other, a, b, side_nodes)) {
          found = true;
          break;
        }
      }

      if (!found) {
        element_closed = false;
        break;
      }
    }

    closed = closed && element_closed;
  }

  return closed;
}

}  // namespace

bool surface_is_closed(
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_nodes, const count_t *const SMESH_RESTRICT n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx) {
  SMESH_UNUSED(n_nodes);

  switch (element_type) {
    case TRI3:
    case TRISHELL3:
      return surface_is_closed_edges(n_elements, elems, n2e_ptr, n2e_idx,
                                     tri3_side_nodes);
    case QUAD4:
    case QUADSHELL4:
      return surface_is_closed_edges(n_elements, elems, n2e_ptr, n2e_idx,
                                     quad4_side_nodes);
    default:
      return false;
  }
}

} // namespace smesh
