#include "smesh_path.hpp"
#include "smesh_read.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(IDX_T, GEOM_T)             \
  template int mesh_from_folder<IDX_T, GEOM_T>(const Path &, int *,            \
                                               ptrdiff_t *, IDX_T ***, int *,  \
                                               ptrdiff_t *, GEOM_T ***);

#define SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(T)               \
  template int array_read_convert_from_extension<T>(const Path &path,         \
                                                    T **data,                 \
                                                    ptrdiff_t *n_elements);

namespace smesh {
  // SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER(i64, f64);

// Called directly by Mesh::read (smesh_mesh.cpp). Without these that call relied
// on mesh_from_folder above leaving an out-of-line copy behind, which g++ 14 at
// -O3 does not: it inlines the function and nothing links.
template int mesh_block_from_folder<i32>(const Path &, int *, i32 ***, ptrdiff_t *);
template int mesh_block_from_folder<i64>(const Path &, int *, i64 ***, ptrdiff_t *);

SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(f16);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(f32);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(f64);

SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(i8);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(i16);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(i32);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(i64);

SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(u8);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(u16);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(u32);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(u64);

SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(char);

#if defined(__clang__)
SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION(long);
#endif
} // namespace smesh


#undef SMESH_EXPLICIT_INSTANTIATE_MESH_FROM_FOLDER
#undef SMESH_EXPLICIT_INSTANTIATE_ARRAY_READ_FROM_EXTENSION
