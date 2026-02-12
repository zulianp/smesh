#include "smesh_reorder.impl.hpp"
#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_REORDER(IDX_T, ELEM_IDX_T)        \
  template int mesh_block_reorder<IDX_T, ELEM_IDX_T>(                           \
      const int, const ptrdiff_t,                                               \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                  \
      const ELEM_IDX_T *const SMESH_RESTRICT,                                   \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT)

#define SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_RENUMBER_NODES(IDX_T)             \
  template int mesh_block_renumber_element_nodes<IDX_T>(                        \
      const int, const ptrdiff_t,                                               \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                        \
      IDX_T *const SMESH_RESTRICT, IDX_T *const SMESH_RESTRICT)

#define SMESH_EXPLICIT_INSTANTIATE_REORDER_ARRAY(IDX_T, T)                      \
  template int reorder_scatter<IDX_T, T>(                                       \
      const ptrdiff_t, const IDX_T *const SMESH_RESTRICT,                       \
      const T *const SMESH_RESTRICT, T *const SMESH_RESTRICT);                  \
  template int reorder_gather<IDX_T, T>(                                        \
      const ptrdiff_t, const IDX_T *const SMESH_RESTRICT,                       \
      const T *const SMESH_RESTRICT, T *const SMESH_RESTRICT)

SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_REORDER(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_REORDER(i64, i32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_REORDER(i64, i64);

SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_RENUMBER_NODES(i32);
SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_RENUMBER_NODES(i64);

SMESH_EXPLICIT_INSTANTIATE_REORDER_ARRAY(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_REORDER_ARRAY(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_REORDER_ARRAY(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_REORDER_ARRAY(i64, f64);

#undef SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_REORDER
#undef SMESH_EXPLICIT_INSTANTIATE_MESH_BLOCK_RENUMBER_NODES
#undef SMESH_EXPLICIT_INSTANTIATE_REORDER_ARRAY

} // namespace smesh

