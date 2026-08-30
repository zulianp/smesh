#ifndef SMESH_SSQUAD4_HPP
#define SMESH_SSQUAD4_HPP

#include "smesh_base.hpp"

namespace smesh {

inline int ssquad4_lidx(const int L, const int x, const int y) {
  int Lp1 = L + 1;
  int ret = y * Lp1 + x;

  SMESH_ASSERT(ret < Lp1 * Lp1);
  SMESH_ASSERT(ret >= 0);
  return ret;
}

inline int ssquad4_txe(int level) { return level * level; }

inline int ssquad4_nxedge(int level) { return level > 1 ? (level - 1) : 0; }

inline int ssquad4_nxface(int level) { return level > 1 ? (level - 1) * (level - 1) : 0; }

inline int ssquad4_nxe(int level) {
  const int corners = 4;
  const int edge_nodes = 4 * ssquad4_nxedge(level);
  const int area_nodes = ssquad4_nxface(level);
  return corners + edge_nodes + area_nodes;
}

inline int ssedge_lidx(const int L, const int x) {
  SMESH_ASSERT(x >= 0 && x <= L);
  return x;
}

inline int ssedge_nxe(const int L) { return L + 1; }

inline int ssedge_txe(const int L) { return L; }

} // namespace smesh

#endif // SMESH_SSQUAD4_HPP