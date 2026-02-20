#include "smesh_decompose.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(IDX_T, COUNT_T, ELEM_IDX_T)       \
  template int create_n2e<IDX_T, COUNT_T, ELEM_IDX_T>(                         \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const ptrdiff_t, const int,                                              \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, COUNT_T **,     \
      ELEM_IDX_T **);                                                          \
  template int redistribute_n2e<IDX_T, COUNT_T, ELEM_IDX_T>(                   \
      MPI_Comm, const int, const int, const ptrdiff_t, const ptrdiff_t,        \
      const ptrdiff_t, const COUNT_T *const SMESH_RESTRICT,                    \
      const ELEM_IDX_T *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT, \
      IDX_T **const SMESH_RESTRICT, COUNT_T **const SMESH_RESTRICT,            \
      ELEM_IDX_T **const SMESH_RESTRICT);                                      \
  template int localize_element_indices<IDX_T, COUNT_T, ELEM_IDX_T>(           \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      IDX_T *const *const SMESH_RESTRICT, const ptrdiff_t,                     \
      const COUNT_T *const SMESH_RESTRICT,                                     \
      const ELEM_IDX_T *const SMESH_RESTRICT,                                  \
      const IDX_T *const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT);        \
  template int rearrange_local_nodes<IDX_T, COUNT_T, ELEM_IDX_T>(              \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      const ptrdiff_t, COUNT_T *const SMESH_RESTRICT,                          \
      ELEM_IDX_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT,           \
      IDX_T **const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,           \
      ptrdiff_t *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT);       \
  template int rearrange_local_elements<IDX_T, COUNT_T, ELEM_IDX_T>(           \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      const ptrdiff_t, COUNT_T *const SMESH_RESTRICT,                          \
      ELEM_IDX_T *const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT,          \
      const ptrdiff_t, ptrdiff_t *const SMESH_RESTRICT,                        \
      ELEM_IDX_T *const SMESH_RESTRICT)

#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(IDX_T)                        \
  template int prepare_node_renumbering<IDX_T>(                                \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const IDX_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT);         \
  template int stitch_aura_elements<IDX_T>(                                    \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const IDX_T *const SMESH_RESTRICT, const int, const ptrdiff_t,           \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const ptrdiff_t,      \
      IDX_T **const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT,              \
      ptrdiff_t *const SMESH_RESTRICT);                                        \
  template int collect_ghost_and_aura_import_indices<IDX_T>(                   \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const ptrdiff_t, const IDX_T *const SMESH_RESTRICT,                      \
      const IDX_T *const SMESH_RESTRICT,                                       \
      const ptrdiff_t *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT);     \
  template int determine_ownership<IDX_T>(                                     \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t, \
      const IDX_T *const SMESH_RESTRICT,                                       \
      const ptrdiff_t *const SMESH_RESTRICT, int *const SMESH_RESTRICT);

// NOTE: `expand_aura_elements_inconsistent` currently requires `idx_t` and
// `element_idx_t` to match due to its output pointer type.
#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(IDX_T, COUNT_T)              \
  template int expand_aura_elements_inconsistent<IDX_T, COUNT_T, IDX_T>(       \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const int,                   \
      COUNT_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT,              \
      const IDX_T *const SMESH_RESTRICT,                                       \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const IDX_T *const SMESH_RESTRICT, const ptrdiff_t, const ptrdiff_t,     \
      IDX_T **const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT,              \
      ptrdiff_t *const SMESH_RESTRICT)

namespace smesh {

SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i32, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i32, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i64, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i64, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i32, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i32, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i64, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i64, i64);

SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(i64);

SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i32, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i64, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i64, i64);

} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA
#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX
#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE
