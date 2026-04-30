#include "smesh_distributed_create.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_HEX8_CUBE_CREATE_DISTRIBUTED(IDX_T, GEOM_T)  \
  template int hex8_cube_create_distributed<IDX_T, GEOM_T>(                    \
      MPI_Comm comm, const ptrdiff_t nx, const ptrdiff_t ny,                   \
      const ptrdiff_t nz, const GEOM_T xmin, const GEOM_T ymin,                \
      const GEOM_T zmin, const GEOM_T xmax, const GEOM_T ymax,                 \
      const GEOM_T zmax, int *nnodesxelem_out,                                 \
      ptrdiff_t *n_local_elements_out, ptrdiff_t *n_global_elements_out,       \
      IDX_T ***elems_out, int *spatial_dim_out,                                \
      ptrdiff_t *n_local_nodes_out, ptrdiff_t *n_global_nodes_out,             \
      GEOM_T ***points_out)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_HEX8_CUBE_CREATE_DISTRIBUTED(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_HEX8_CUBE_CREATE_DISTRIBUTED(i64, f32);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_HEX8_CUBE_CREATE_DISTRIBUTED
