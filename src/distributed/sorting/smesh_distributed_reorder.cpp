#include "smesh_distributed_reorder.impl.hpp"

#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(IDX_T, GEOM_T)          \
  template int distributed_reorder_elements<IDX_T, GEOM_T,                     \
                                            ElementOrdering<GEOM_T>>(          \
      MPI_Comm, const int, const ptrdiff_t, const ptrdiff_t,                   \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const ptrdiff_t,      \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                      \
      large_idx_t *const SMESH_RESTRICT,                                       \
      ElementOrdering<GEOM_T>)

#define SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER_REORDERED(                 \
    IDX_T, GEOM_T, LARGE_IDX_T)                                                \
  template int mesh_from_folder_reordered<                                     \
      IDX_T, GEOM_T, LARGE_IDX_T, ElementOrdering<GEOM_T>>(                    \
      const MPI_Comm comm, const Path &folder, int *nnodesxelem_out,           \
      ptrdiff_t *n_global_elements_out, ptrdiff_t *n_owned_elements_out,       \
      ptrdiff_t *n_shared_elements_out, ptrdiff_t *n_ghost_elements_out,       \
      LARGE_IDX_T **element_mapping_out,                                       \
      LARGE_IDX_T **aura_element_mapping_out, IDX_T ***elements_out,           \
      int *spatial_dim_out, ptrdiff_t *n_global_nodes_out,                     \
      ptrdiff_t *n_owned_nodes_out, ptrdiff_t *n_shared_nodes_out,             \
      ptrdiff_t *n_ghost_nodes_out, ptrdiff_t *n_aura_nodes_out,               \
      LARGE_IDX_T **node_mapping_out, GEOM_T ***points_out,                    \
      int **node_owner_out, ptrdiff_t **node_offsets_out, IDX_T **ghosts_out,  \
      ElementOrdering<GEOM_T> ordering)

SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER(i64, f64);

SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER_REORDERED(i32, f32, i32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER_REORDERED(i32, f32, i64);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER_REORDERED(i64, f32, i64);

#undef SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER_REORDERED
#undef SMESH_EXPLICIT_INSTANTIATE_DISTRIBUTED_REORDER

} // namespace smesh
