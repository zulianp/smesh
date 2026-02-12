#include "smesh_promotions.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(IDX_T, COUNT_T, ELEM_IDX_T,      \
                                              GEOM_T)                          \
  template void mesh_tet4_to_tet15<IDX_T, COUNT_T, ELEM_IDX_T>(                \
      const ptrdiff_t, const ptrdiff_t,                                        \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
      const COUNT_T *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,  \
      const ELEM_IDX_T *const SMESH_RESTRICT,                                  \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      ptrdiff_t *const SMESH_RESTRICT);                                        \
  template void mesh_tet4_to_tet15_points<IDX_T, COUNT_T, GEOM_T>(             \
      const ptrdiff_t, const ptrdiff_t,                                        \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                \
      const COUNT_T *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,  \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT);                     \
  template int p1_to_p2<COUNT_T, IDX_T, GEOM_T>(                               \
      const enum ElemType element_type, const ptrdiff_t n_elements,            \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT p1_elements,     \
      const int spatial_dim, const ptrdiff_t p1_n_nodes,                       \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT p1_points,      \
      const COUNT_T *const SMESH_RESTRICT n2n_ptr,                             \
      const IDX_T *const SMESH_RESTRICT n2n_idx,                               \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT p2_elements,           \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT p2_points);           \
  template void quad4_to_hex8_extrude<IDX_T, GEOM_T>(                          \
      const ptrdiff_t nsides, const ptrdiff_t nnodes,                          \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT quad4_elements,  \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const ptrdiff_t nlayers, const GEOM_T height,                            \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements,         \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT extruded_points);     \
  template int tri3_to_wedge6_extrude<IDX_T, GEOM_T>(                          \
      const ptrdiff_t nsides, const ptrdiff_t nnodes,                          \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT tri3_elements,  \
      const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT points,         \
      const ptrdiff_t nlayers, const GEOM_T height,                            \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT wedge6_elements,       \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT extruded_points);

namespace smesh {
// Match the common explicit-instantiation convention used across the repo.
SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(i16, i32, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(i32, i32, i32, f32);
SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS(i64, i32, i32, f32);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_PROMOTIONS