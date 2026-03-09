#include "smesh_sshex8_graph.hpp"
#include "smesh_sshex8_graph.impl.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SINGLE(IDX_T)                           \
  template int sshex8_generate_elements<IDX_T>(                                    \
      const int, const ptrdiff_t, const ptrdiff_t,                                \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                    \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, ptrdiff_t *, ptrdiff_t *); \
  template int sshex8_hierarchical_renumbering<IDX_T>(                             \
      const int, const int, int *const, const ptrdiff_t, const ptrdiff_t,         \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT); \
  template int ssquad4_hierarchical_remapping<IDX_T>(                              \
      const int, const int, int *const, const ptrdiff_t, const ptrdiff_t,         \
    IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                    \
      IDX_T **, ptrdiff_t *)

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_GRAPH(IDX_T, ELEM_IDX_T, COUNT_T)       \
  template int sshex8_crs_graph<ELEM_IDX_T, COUNT_T, IDX_T>(                       \
      const int, const ptrdiff_t, const ptrdiff_t,                                 \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                     \
      COUNT_T **, IDX_T **)

#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SIDESET(IDX_T, ELEM_IDX_T)               \
  template int sshex8_extract_surface_from_sideset<IDX_T, ELEM_IDX_T>(             \
      const int, const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,          \
      const ptrdiff_t, const ELEM_IDX_T *const SMESH_RESTRICT,                     \
      const i16 *const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT);              \
  template int sshex8_extract_nodeset_from_sideset<IDX_T, ELEM_IDX_T>(             \
      const int, const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,          \
      const ptrdiff_t, const ELEM_IDX_T *const SMESH_RESTRICT,                     \
      const i16 *const SMESH_RESTRICT, ptrdiff_t *, IDX_T **SMESH_RESTRICT)

SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SINGLE(i32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SINGLE(i64);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_GRAPH(i32, i32, i32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_GRAPH(i64, i64, i64);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SIDESET(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SIDESET(i64, i64);

#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SIDESET
#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_GRAPH
#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_SINGLE

} // namespace smesh