#ifndef SMESH_SSQUAD4_MESH_IMPL_HPP
#define SMESH_SSQUAD4_MESH_IMPL_HPP

#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_mesh.hpp"

namespace smesh {

template <typename idx_t>
int ssquad4_to_standard_quad4_mesh(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements) {
  const int txe = ssquad4_txe(level);

  int lnode[4];
  for (int yi = 0; yi < level; yi++) {
    for (int xi = 0; xi < level; xi++) {
      lnode[0] = ssquad4_lidx(level, xi, yi);
      lnode[1] = ssquad4_lidx(level, xi + 1, yi);
      lnode[2] = ssquad4_lidx(level, xi + 1, yi + 1);
      lnode[3] = ssquad4_lidx(level, xi, yi + 1);

      int le = yi * level + xi;
      SMESH_ASSERT(le < txe);

      for (int l = 0; l < 4; l++) {
        for (ptrdiff_t e = 0; e < nelements; e++) {
          idx_t node = elements[lnode[l]][e];
          quad4_elements[l][e * txe + le] = node;
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_SSQUAD4_MESH_IMPL_HPP