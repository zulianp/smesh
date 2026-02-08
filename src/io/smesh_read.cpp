#include "smesh_path.hpp"
#include "smesh_read.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(IDX_T, GEOM_T)             \
  template int mesh_from_folder<IDX_T, GEOM_T>(const Path &, int *,            \
                                               ptrdiff_t *, IDX_T ***, int *,  \
                                               ptrdiff_t *, GEOM_T ***);

namespace smesh {
  // SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i64, f32);
} // namespace smesh


#undef SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER