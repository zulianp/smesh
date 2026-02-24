#include "smesh_distributed_read.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(IDX_T, GEOM_T,             \
                                                    LARGE_IDX_T)               \
  template int mesh_from_folder<IDX_T, GEOM_T, LARGE_IDX_T>(                   \
      const MPI_Comm comm, const Path &folder, int *nnodesxelem_out,           \
      ptrdiff_t *n_global_elements_out, ptrdiff_t *n_owned_elements_out,       \
      ptrdiff_t *n_shared_elements_out,                                         \
      ptrdiff_t *n_ghost_elements_out, LARGE_IDX_T **element_mapping_out,      \
      IDX_T ***elements_out, int *spatial_dim_out,                             \
      ptrdiff_t *n_global_nodes_out, ptrdiff_t *n_owned_nodes_out,             \
      ptrdiff_t *n_shared_nodes_out,                                            \
      ptrdiff_t *n_ghost_nodes_out, LARGE_IDX_T **node_mapping_out,            \
      GEOM_T ***points_out, int **node_owner_out,                              \
      ptrdiff_t **node_offsets_out, IDX_T **ghosts_out)

namespace smesh {
// Explicit instantiation
template int read_mapped_field<i32>(MPI_Comm comm, const char *input_path,
                                    ptrdiff_t n_local, ptrdiff_t n_global,
                                    const i32 *SMESH_RESTRICT mapping,
                                    MPI_Datatype data_type,
                                    void *SMESH_RESTRICT data_out);

template int read_mapped_field<i64>(MPI_Comm comm, const char *input_path,
                                    ptrdiff_t n_local, ptrdiff_t n_global,
                                    const i64 *SMESH_RESTRICT mapping,
                                    MPI_Datatype data_type,
                                    void *SMESH_RESTRICT data_out);

template int mesh_coordinates_from_folder<f32>(
    MPI_Comm comm, const Path &folder, int *spatial_dim_out, f32 ***points_out,
    ptrdiff_t *n_local_nodes_out, ptrdiff_t *n_global_nodes_out);

template int mesh_block_from_folder(MPI_Comm comm, const Path &folder,
                                    int *nnodesxelem_out, i32 ***const elems,
                                    ptrdiff_t *const n_local_elements_out,
                                    ptrdiff_t *const n_global_elements_out);

SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i32, f32, i32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i32, f32, i64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER
