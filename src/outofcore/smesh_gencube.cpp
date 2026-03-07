#include "smesh_gencube.impl.hpp"

namespace smesh {
template int mesh_hex8_cube_to_folder<i32, f32>(
    const Path &folder, const ptrdiff_t nx, const ptrdiff_t ny,
    const ptrdiff_t nz, const f32 xmin, const f32 ymin, const f32 zmin,
    const f32 xmax, const f32 ymax, const f32 zmax,
    const ptrdiff_t z_chunk_size);

template int mesh_hex8_cube_to_folder<i64, f32>(
    const Path &folder, const ptrdiff_t nx, const ptrdiff_t ny,
    const ptrdiff_t nz, const f32 xmin, const f32 ymin, const f32 zmin,
    const f32 xmax, const f32 ymax, const f32 zmax,
    const ptrdiff_t z_chunk_size);

} // namespace smesh
