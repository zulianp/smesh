#include "smesh_ssquad4_graph.impl.hpp"
#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(T)                           \
  template int quad4_build_edge_graph_from_n2e<T, T, element_idx_t>(           \
      ptrdiff_t, ptrdiff_t,                                                    \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT,                     \
      const T *const SMESH_RESTRICT,                                           \
      const element_idx_t *const SMESH_RESTRICT, T **, T **);

SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(i32);
SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(i64);
SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(i16);

#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH

} // namespace smesh
