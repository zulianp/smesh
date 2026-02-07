#include <cstddef>
#include <cstdint>

#include "smesh_types.hpp"
#include "smesh_adjaciency.impl.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_ADJACENCY(T)                                 \
  template void create_element_adj_table_from_dual_graph<T, T, T>(              \
      const ptrdiff_t, enum ElemType,                                           \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT, const T *const,            \
      const T *const, T *const SMESH_RESTRICT);                                 \
  template void create_element_adj_table_from_dual_graph_soa<T, T, T>(          \
      const ptrdiff_t, enum ElemType,                                           \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT, const T *const,            \
      const T *const,                                                          \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT);                           \
  template void create_element_adj_table<T, T, T>(                              \
      const ptrdiff_t, const ptrdiff_t, enum ElemType,                          \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT, T **SMESH_RESTRICT);       \
  template void extract_surface_connectivity_with_adj_table<T, T, T>(           \
      const ptrdiff_t, const ptrdiff_t, enum ElemType,                          \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT, ptrdiff_t *,               \
      T **SMESH_RESTRICT, T **SMESH_RESTRICT);                                  \
  template int extract_sideset_from_adj_table<T>(                               \
      const enum ElemType, const ptrdiff_t, const T *const SMESH_RESTRICT,      \
      ptrdiff_t *SMESH_RESTRICT, T **SMESH_RESTRICT, int16_t **SMESH_RESTRICT); \
  template int extract_skin_sideset<T, T, T>(                                   \
      const ptrdiff_t, const ptrdiff_t, const enum ElemType,                    \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT, ptrdiff_t *SMESH_RESTRICT, \
      T **SMESH_RESTRICT, int16_t **SMESH_RESTRICT);                            \
  template int extract_surface_from_sideset<T, T, T>(                           \
      const enum ElemType, T *const SMESH_RESTRICT *const SMESH_RESTRICT,       \
      const ptrdiff_t, const T *const SMESH_RESTRICT,                           \
      const int16_t *const SMESH_RESTRICT,                                      \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT);                           \
  template int extract_nodeset_from_sideset<T, T, T>(                           \
      const enum ElemType, T *const SMESH_RESTRICT *const SMESH_RESTRICT,       \
      const ptrdiff_t, const T *const SMESH_RESTRICT,                           \
      const int16_t *const SMESH_RESTRICT, ptrdiff_t *,                         \
      T **SMESH_RESTRICT);                                                      \
  template int extract_nodeset_from_sidesets<T, T, T>(                          \
      ptrdiff_t, const enum ElemType[], T **const SMESH_RESTRICT[],             \
      const ptrdiff_t[], const T *const SMESH_RESTRICT[],                       \
      const int16_t *const SMESH_RESTRICT[], ptrdiff_t *,                       \
      T **SMESH_RESTRICT)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_ADJACENCY(i16);
SMESH_EXPLICIT_INSTANTIATE_ADJACENCY(i32);
SMESH_EXPLICIT_INSTANTIATE_ADJACENCY(i64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_ADJACENCY
