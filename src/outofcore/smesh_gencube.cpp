#include "smesh_gencube.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_GENCUBE(IDX_T, GEOM_T)                      \
  template int mesh_hex8_cube_to_folder<IDX_T, GEOM_T>(                        \
      const Path &folder, const ptrdiff_t nx, const ptrdiff_t ny,              \
      const ptrdiff_t nz, const GEOM_T xmin, const GEOM_T ymin,                \
      const GEOM_T zmin, const GEOM_T xmax, const GEOM_T ymax,                 \
      const GEOM_T zmax, const ptrdiff_t z_chunk_size)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_GENCUBE(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_GENCUBE(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_GENCUBE(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_GENCUBE(i64, f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_GENCUBE
