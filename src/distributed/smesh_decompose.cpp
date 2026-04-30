#include "smesh_decompose.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(IDX_T, COUNT_T, ELEM_IDX_T, LARGE_IDX_T)       \
  template int create_n2e<IDX_T, COUNT_T, ELEM_IDX_T>(                         \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const ptrdiff_t, const int,                                              \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, COUNT_T **,     \
      ELEM_IDX_T **);                                                          \
  template int localize_element_indices<IDX_T, COUNT_T, ELEM_IDX_T, LARGE_IDX_T>( \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      IDX_T *const *const SMESH_RESTRICT, const ptrdiff_t,                     \
      const COUNT_T *const SMESH_RESTRICT,                                     \
      const ELEM_IDX_T *const SMESH_RESTRICT,                                  \
      const LARGE_IDX_T *const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT);  \
  template int rearrange_local_nodes<IDX_T, COUNT_T, ELEM_IDX_T, LARGE_IDX_T>( \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      const ptrdiff_t, COUNT_T *const SMESH_RESTRICT,                          \
      ELEM_IDX_T *const SMESH_RESTRICT, LARGE_IDX_T *const SMESH_RESTRICT,     \
      IDX_T **const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,           \
      ptrdiff_t *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT);       \
  template int rearrange_local_elements<IDX_T, COUNT_T, ELEM_IDX_T, LARGE_IDX_T>(           \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      const ptrdiff_t, COUNT_T *const SMESH_RESTRICT,                          \
      ELEM_IDX_T *const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT,          \
      const ptrdiff_t, ptrdiff_t *const SMESH_RESTRICT,                        \
      LARGE_IDX_T *const SMESH_RESTRICT,                                       \
      const LARGE_IDX_T *const SMESH_RESTRICT)

#define SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(COUNT_T, ELEM_IDX_T, LARGE_IDX_T) \
  template int redistribute_n2e<COUNT_T, ELEM_IDX_T, LARGE_IDX_T>(               \
      MPI_Comm, const int, const int, const ptrdiff_t, const ptrdiff_t,          \
      const ptrdiff_t, const COUNT_T *const SMESH_RESTRICT,                      \
      const ELEM_IDX_T *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,   \
      LARGE_IDX_T **const SMESH_RESTRICT, COUNT_T **const SMESH_RESTRICT,        \
      ELEM_IDX_T **const SMESH_RESTRICT)

#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(IDX_T, L2G_T)                 \
  template int prepare_node_renumbering<IDX_T, L2G_T>(                         \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const L2G_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT);         \
  template int stich_aura_elements<IDX_T, L2G_T>(                             \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const L2G_T *const SMESH_RESTRICT, const int, const ptrdiff_t,           \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const ptrdiff_t,      \
      IDX_T **const SMESH_RESTRICT, L2G_T **const SMESH_RESTRICT,              \
      ptrdiff_t *const SMESH_RESTRICT);                                        \
  template int collect_ghost_and_aura_import_indices<IDX_T, L2G_T>(            \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const ptrdiff_t, const L2G_T *const SMESH_RESTRICT,                      \
      const IDX_T *const SMESH_RESTRICT,                                       \
      const ptrdiff_t *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT);     \
  template int group_ghost_and_aura_by_rank<IDX_T, L2G_T>(                     \
      const int, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      L2G_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT, int *const SMESH_RESTRICT, \
      const int, const ptrdiff_t, const ptrdiff_t, IDX_T **const SMESH_RESTRICT)


#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(IDX_T, COUNT_T, LARGE_IDX_T) \
  template int expand_aura_elements_inconsistent<IDX_T, COUNT_T, IDX_T, LARGE_IDX_T>( \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const int,                   \
      COUNT_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT,              \
      const LARGE_IDX_T *const SMESH_RESTRICT,                                 \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const LARGE_IDX_T *const SMESH_RESTRICT, const ptrdiff_t, const ptrdiff_t, \
      LARGE_IDX_T **const SMESH_RESTRICT, IDX_T **const SMESH_RESTRICT,        \
      ptrdiff_t *const SMESH_RESTRICT)

namespace smesh {

SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i32, i32, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i32, i64, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i64, i32, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32, i64, i64, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i32, i32, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i32, i64, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i64, i32, large_idx_t);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64, i64, i64, large_idx_t);

SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(i32, i32, i32);
SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(i32, i32, i64);
SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(i32, i64, i64);
SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(i64, i64, i32);
SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(i64, i32, i64);
SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE(i64, i64, i64);

template int determine_ownership<i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,
    const i32 *const SMESH_RESTRICT, const ptrdiff_t *const SMESH_RESTRICT,
    int *const SMESH_RESTRICT);
template int determine_ownership<i64>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,
    const i64 *const SMESH_RESTRICT, const ptrdiff_t *const SMESH_RESTRICT,
    int *const SMESH_RESTRICT);

template int rearrange_local_elements<i32, i32, i32, i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const int,
    const ptrdiff_t, i32 *const SMESH_RESTRICT, i32 *const SMESH_RESTRICT,
    i32 **const SMESH_RESTRICT, const ptrdiff_t, ptrdiff_t *const SMESH_RESTRICT,
    i32 *const SMESH_RESTRICT, const i32 *const SMESH_RESTRICT);
template int localize_element_indices<i32, i32, i32, i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const int,
    i32 *const *const SMESH_RESTRICT, const ptrdiff_t,
    const i32 *const SMESH_RESTRICT, const i32 *const SMESH_RESTRICT,
    const i32 *const SMESH_RESTRICT, i32 **const SMESH_RESTRICT);
template int rearrange_local_nodes<i32, i32, i32, i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const int,
    const ptrdiff_t, i32 *const SMESH_RESTRICT, i32 *const SMESH_RESTRICT,
    i32 *const SMESH_RESTRICT, i32 **const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT);
template int localize_element_indices<i32, i64, i64, i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const int,
    i32 *const *const SMESH_RESTRICT, const ptrdiff_t,
    const i64 *const SMESH_RESTRICT, const i64 *const SMESH_RESTRICT,
    const i32 *const SMESH_RESTRICT, i32 **const SMESH_RESTRICT);
template int rearrange_local_nodes<i32, i64, i64, i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const int,
    const ptrdiff_t, i64 *const SMESH_RESTRICT, i64 *const SMESH_RESTRICT,
    i32 *const SMESH_RESTRICT, i32 **const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT);
template int rearrange_local_elements<i32, i64, i64, i32>(
    const int, const int, const ptrdiff_t, const ptrdiff_t, const int,
    const ptrdiff_t, i64 *const SMESH_RESTRICT, i64 *const SMESH_RESTRICT,
    i32 **const SMESH_RESTRICT, const ptrdiff_t, ptrdiff_t *const SMESH_RESTRICT,
    i32 *const SMESH_RESTRICT, const i32 *const SMESH_RESTRICT);

SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(i32, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX(i64, i64);

SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i32, i32, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i32, i64, i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i32, i32, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i32, i64, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i64, i32, i64);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA(i64, i64, i64);
template int expand_aura_elements_inconsistent<i64, i32, i32, i64>(
    MPI_Comm, const ptrdiff_t, const ptrdiff_t, const int,
    i32 *const SMESH_RESTRICT, i32 *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT, const ptrdiff_t, const ptrdiff_t,
    i64 **const SMESH_RESTRICT, i64 **const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT);
template int expand_aura_elements_inconsistent<i64, i64, i32, i64>(
    MPI_Comm, const ptrdiff_t, const ptrdiff_t, const int,
    i64 *const SMESH_RESTRICT, i32 *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT, const ptrdiff_t, const ptrdiff_t,
    i64 **const SMESH_RESTRICT, i64 **const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT);
template int expand_aura_elements_inconsistent<i32, i64, i64, i32>(
    MPI_Comm, const ptrdiff_t, const ptrdiff_t, const int,
    i64 *const SMESH_RESTRICT, i64 *const SMESH_RESTRICT,
    const i32 *const SMESH_RESTRICT,
    const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    const i32 *const SMESH_RESTRICT, const ptrdiff_t, const ptrdiff_t,
    i32 **const SMESH_RESTRICT, i32 **const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT);
template int expand_aura_elements_inconsistent<i32, i64, i64, i64>(
    MPI_Comm, const ptrdiff_t, const ptrdiff_t, const int,
    i64 *const SMESH_RESTRICT, i64 *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT,
    const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    const i64 *const SMESH_RESTRICT, const ptrdiff_t, const ptrdiff_t,
    i64 **const SMESH_RESTRICT, i32 **const SMESH_RESTRICT,
    ptrdiff_t *const SMESH_RESTRICT);

} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_AURA
#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE_IDX
#undef SMESH_EXPLICIT_INSTANTIATE_REDISTRIBUTE
#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE
