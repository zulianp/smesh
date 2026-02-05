#include <cstddef>
#include <cstdint>

#include "smesh_crs_graph.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_CRS_GRAPH(T)                                \
  template T find_idx<T>(T, const T *, T);                                     \
  template int create_n2e<T, T, T>(                                            \
      ptrdiff_t, ptrdiff_t, int,                                               \
      T *const SMESH_RESTRICT * const SMESH_RESTRICT, T **, T **);             \
  template int create_crs_graph_for_elem_type<T, T>(                            \
      enum ElemType, ptrdiff_t, ptrdiff_t,                                     \
      T *const SMESH_RESTRICT * const SMESH_RESTRICT, T **, T **);             \
  template int create_crs_graph<T, T>(                                         \
      ptrdiff_t, ptrdiff_t, T *const SMESH_RESTRICT * const SMESH_RESTRICT,    \
      T **, T **);                                                            \
  template int create_crs_graph_3<T, T>(                                       \
      ptrdiff_t, ptrdiff_t, T *const SMESH_RESTRICT * const SMESH_RESTRICT,    \
      T **, T **);                                                            \
  template int block_crs_to_crs<T, T, T>(ptrdiff_t, int, const T *, const T *, \
                                        const T *, T *, T *, T *);            \
  template int crs_to_coo<T, T>(ptrdiff_t, const T *, T *);                    \
  template int sorted_coo_to_crs<T, T>(T, const T *, ptrdiff_t, T *);          \
  template int create_dual_graph<T, T, T>(ptrdiff_t, ptrdiff_t, enum ElemType, \
                                         T *const SMESH_RESTRICT *            \
                                             const SMESH_RESTRICT,            \
                                         T **, T **);                         \
  template int crs_graph_block_to_scalar<T, T>(ptrdiff_t, int, const T *,      \
                                               const T *, T *, T *);          \
  template int create_crs_graph_from_element<T, T>(                            \
      ptrdiff_t, ptrdiff_t, int,                                               \
      T *const SMESH_RESTRICT * const SMESH_RESTRICT, T **, T **);             \
  template int create_crs_graph_upper_triangular_from_element<T, T>(           \
      ptrdiff_t, ptrdiff_t, int,                                               \
      T *const SMESH_RESTRICT * const SMESH_RESTRICT, T **, T **)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_CRS_GRAPH(std::int32_t);
SMESH_EXPLICIT_INSTANTIATE_CRS_GRAPH(std::int64_t);
SMESH_EXPLICIT_INSTANTIATE_CRS_GRAPH(std::int16_t);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_CRS_GRAPH
