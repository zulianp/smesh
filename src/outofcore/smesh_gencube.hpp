#ifndef SMESH_GENCUBE_HPP
#define SMESH_GENCUBE_HPP

#include "smesh_base.hpp"
#include "smesh_path.hpp"
#include "smesh_types.hpp"

#include <limits>

namespace smesh {

template <typename idx_t, typename geom_t>
int mesh_hex8_cube_to_folder(const Path &folder, const ptrdiff_t nx,
                             const ptrdiff_t ny, const ptrdiff_t nz,
                             const geom_t xmin, const geom_t ymin,
                             const geom_t zmin, const geom_t xmax,
                             const geom_t ymax, const geom_t zmax,
                             const ptrdiff_t z_chunk_size = (ptrdiff_t)(std::numeric_limits<int>::max()));

} // namespace smesh

#endif // SMESH_GENCUBE_HPP