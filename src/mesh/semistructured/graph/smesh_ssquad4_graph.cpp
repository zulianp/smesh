#include "smesh_ssquad4_graph.hpp"
#include "smesh_ssquad4_graph.impl.hpp"
#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(T)                           \
  template int quad4_build_edge_graph_from_n2e<T, T, element_idx_t>(           \
      ptrdiff_t, ptrdiff_t,                                                    \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT,                     \
      const T *const SMESH_RESTRICT,                                           \
      const element_idx_t *const SMESH_RESTRICT, T **, T **);

#define SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_SINGLE(IDX_T)                                                    \
    template int ssquad4_generate_elements<IDX_T>(const int,                                                \
                                                  const ptrdiff_t,                                          \
                                                  const ptrdiff_t,                                          \
                                                  const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,  \
                                                  IDX_T *const SMESH_RESTRICT *const       SMESH_RESTRICT,  \
                                                  ptrdiff_t *,                                              \
                                                  ptrdiff_t *);                                             \
    template int ssquad4_generate_elements_blocks<IDX_T>(                                                   \
            const int,                                                                                      \
            const ptrdiff_t,                                                                                \
            const ptrdiff_t *const,                                                                         \
            const ptrdiff_t,                                                                                \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                  \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                        \
            ptrdiff_t *,                                                                                    \
            ptrdiff_t *);                                                                                   \
    template int ssquad4_hierarchical_renumbering<IDX_T>(const int,                                         \
                                                         const int,                                         \
                                                         int *const,                                        \
                                                         const ptrdiff_t,                                   \
                                                         const ptrdiff_t,                                   \
                                                         IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                         IDX_T *const                       SMESH_RESTRICT, \
                                                         const bool);                                       \
    template int ssquad4_hierarchical_renumbering_blocks<IDX_T>(                                            \
            const int,                                                                                      \
            const int,                                                                                      \
            int *const,                                                                                     \
            const ptrdiff_t,                                                                                \
            const ptrdiff_t *const,                                                                         \
            const ptrdiff_t,                                                                                \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                        \
            IDX_T *const SMESH_RESTRICT,                                                                    \
            const bool);                                                                                    \
    template int ssquad4_extract_nodeset_from_sideset<IDX_T, IDX_T>(                                        \
            const int,                                                                                      \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                        \
            const ptrdiff_t,                                                                                \
            const IDX_T *const SMESH_RESTRICT,                                                              \
            const i16 *const SMESH_RESTRICT,                                                                \
            ptrdiff_t *,                                                                                    \
            IDX_T **SMESH_RESTRICT)

SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(i32);
SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(i64);
SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH(i16);
SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_SINGLE(i32);
SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_SINGLE(i64);

#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_GRAPH
#undef SMESH_EXPLICIT_INSTANTIATE_SSQUAD4_SINGLE

}  // namespace smesh

