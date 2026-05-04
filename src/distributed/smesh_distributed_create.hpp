#ifndef SMESH_DISTRIBUTED_CREATE_HPP
#define SMESH_DISTRIBUTED_CREATE_HPP

#include "smesh_distributed_base.hpp"

#include <mpi.h>
#include <stddef.h>

namespace smesh {

// Generate a hex8 cube mesh distributed across `comm` in the same per-rank
// layout that `mesh_block_from_folder` + `mesh_coordinates_from_folder`
// produce when reading a serialized hex8 cube. The output is ready to be
// passed directly into `mesh_create_parallel`.
//
// Per-rank chunk convention: matches `rank_split`/`rank_start`
// (which is also matrixio's `array_create_from_file` chunk convention).
//
// Global sizing:
//   n_global_elements = nx * ny * nz
//   n_global_nodes    = (nx+1) * (ny+1) * (nz+1)
//
// Output ownership: caller frees `(*elems_out)[d]` and `(*elems_out)`,
// `(*points_out)[d]` and `(*points_out)`.
template <typename idx_t, typename geom_t>
int hex8_cube_create_distributed(
    MPI_Comm comm, const ptrdiff_t nx, const ptrdiff_t ny, const ptrdiff_t nz,
    const geom_t xmin, const geom_t ymin, const geom_t zmin, const geom_t xmax,
    const geom_t ymax, const geom_t zmax,
    // Elements
    int *nnodesxelem_out, ptrdiff_t *n_local_elements_out,
    ptrdiff_t *n_global_elements_out, idx_t ***elems_out,
    // Nodes
    int *spatial_dim_out, ptrdiff_t *n_local_nodes_out,
    ptrdiff_t *n_global_nodes_out, geom_t ***points_out);

} // namespace smesh

#endif // SMESH_DISTRIBUTED_CREATE_HPP
