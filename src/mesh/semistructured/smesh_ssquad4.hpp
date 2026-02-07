#ifndef SMESH_SSQUAD4_HPP
#define SMESH_SSQUAD4_HPP

#include "smesh_base.hpp"

namespace smesh {

static SMESH_INLINE int ssquad4_lidx(const int L, const int x, const int y) {
  int Lp1 = L + 1;
  int ret = y * Lp1 + x;

  SMESH_ASSERT(ret < Lp1 * Lp1);
  SMESH_ASSERT(ret >= 0);
  return ret;
}

static SMESH_INLINE int ssquad4_txe(int level) { return level * level; }

static SMESH_INLINE int ssquad4_nxe(int level) {
  const int corners = 4;
  const int edge_nodes = 4 * (level - 1);
  const int area_nodes = (level - 1) * (level - 1);
  return corners + edge_nodes + area_nodes;
}

} // namespace smesh

#endif // SMESH_SSQUAD4_HPP