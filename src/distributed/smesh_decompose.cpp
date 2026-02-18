#include "smesh_decompose.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(T)                                \
  template int create_n2e<T, T, T>(                                            \
      MPI_Comm, const ptrdiff_t, const ptrdiff_t, const ptrdiff_t,             \
      const ptrdiff_t, const int,                                              \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT, T **, T **);        \
  template int redistribute_n2e<T, T, T>(                                      \
      MPI_Comm, const int, const int, const ptrdiff_t, const ptrdiff_t,        \
      const ptrdiff_t, const T *const SMESH_RESTRICT,                          \
      const T *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,          \
      T **const SMESH_RESTRICT, T **const SMESH_RESTRICT,                      \
      T **const SMESH_RESTRICT);                                               \
  template int localize_element_indices<T, T, T>(                              \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      T *const *const SMESH_RESTRICT, const ptrdiff_t,                         \
      const T *const SMESH_RESTRICT, const T *const SMESH_RESTRICT,            \
      const T *const SMESH_RESTRICT, T **const SMESH_RESTRICT);                \
  template int rearrange_local_nodes<T, T, T>(                                 \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      const ptrdiff_t, T *const SMESH_RESTRICT, T *const SMESH_RESTRICT,       \
      T *const SMESH_RESTRICT, T **const SMESH_RESTRICT,                       \
      ptrdiff_t *const SMESH_RESTRICT, ptrdiff_t *const SMESH_RESTRICT,        \
      ptrdiff_t *const SMESH_RESTRICT);                                        \
  template int rearrange_local_elements<T, T, T>(                              \
      const int, const int, const ptrdiff_t, const ptrdiff_t, const int,       \
      const ptrdiff_t, T *const SMESH_RESTRICT, T *const SMESH_RESTRICT,       \
      T **const SMESH_RESTRICT, const ptrdiff_t,                               \
      ptrdiff_t *const SMESH_RESTRICT)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i32);
SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE(i64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_DECOMPOSE