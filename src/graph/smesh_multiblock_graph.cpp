
#include <cstddef>
#include <cstdint>

#include "smesh_multiblock_graph.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(T)                         \
  template int create_multiblock_n2e<T, T, T>(                                 \
      block_idx_t, const enum ElemType *, const ptrdiff_t *,                   \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT *, ptrdiff_t,        \
      block_idx_t **, T **, T **);                                             \
  template int create_multiblock_crs_graph_from_n2e<T, T, T>(                  \
      block_idx_t, const enum ElemType *, const ptrdiff_t *, ptrdiff_t,        \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT *,                   \
      const T *const SMESH_RESTRICT, const T *const SMESH_RESTRICT,            \
      const block_idx_t *const SMESH_RESTRICT, T **, T **);                    \
  template int create_multiblock_crs_graph<T, T, T>(                           \
      block_idx_t, const enum ElemType *, const ptrdiff_t *,                   \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT *, ptrdiff_t, T **,  \
      T **);                                                                   \
  template int create_multiblock_crs_graph_upper_triangular_from_n2e<T, T, T>( \
      block_idx_t, const enum ElemType *, const ptrdiff_t *, ptrdiff_t,        \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT *,                   \
      const T *const SMESH_RESTRICT, const T *const SMESH_RESTRICT,            \
      const block_idx_t *const SMESH_RESTRICT, T **, T **);                    \
  template int create_multiblock_crs_graph_upper_triangular<T, T, T>(          \
      block_idx_t, const enum ElemType *, const ptrdiff_t *,                   \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT *, ptrdiff_t, T **,  \
      T **)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(i32);
SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(i64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH
