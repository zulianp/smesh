#ifndef SMESH_NODE_HEIGHT_HPP
#define SMESH_NODE_HEIGHT_HPP

// Local mesh size per node, for error bounds on coordinates (compression bound mode h):
// h_n = the smallest height of the elements incident to node n.
//
// Why heights: moving a node by dx changes an element's geometry (Jacobian, FFF) by about
// |dx| / height, not |dx| / edge length. Flat elements can have edges far longer than their
// smallest height (graded tet meshes: up to 1e4x), so an edge-based bound does not bound the
// operator error; a height-based one does, and a perturbation below height / 4 per component can
// not invert a simplex.
//
//   TET*, first 4 nodes        3 V / (largest face area)
//   HEX8 / HEX27, corners      the smallest height of the 8 corner tets (corner + its 3 edge
//                              neighbours); 1/sqrt(3) for the unit cube
//   WEDGE6, PYRAMID5           the smallest height of their corner tets
//   TRI*, TRISHELL*            2 A / (longest edge)
//   QUAD*, QUADSHELL*          the smallest height of the 4 corner triangles
//   EDGE*, BEAM2, EDGESHELL*   the length
// All nodes of an element (mid-edge and face nodes too) get the element's height. Computed in
// double; volumes and areas are taken in absolute value, so inverted elements still get a height.
// Other element types (macro, semi-structured, proteus) are not supported.

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

#include <cstddef>

namespace smesh {

/// Whether element_min_height / accumulate_node_heights handle `type`.
bool node_height_supported(enum ElemType type);

/// Smallest height of one element; p[k] are the coordinates (x, y, z) of its nodes in the
/// element's node order (only the corner nodes are used).
double element_min_height(enum ElemType type, const double (*p)[3]);

/// h[n] = min(h[n], height of every element containing n), for the n_elements elements of one
/// block (elements[k][e]: node k of element e). h must be initialised (e.g. to +infinity).
/// Returns SMESH_FAILURE for unsupported element types.
template <typename idx_t, typename geom_t>
int accumulate_node_heights(enum ElemType type, ptrdiff_t n_elements, const idx_t *const *elements, int spatial_dim,
                            const geom_t *const *points, double *h);

}  // namespace smesh

#endif  // SMESH_NODE_HEIGHT_HPP
