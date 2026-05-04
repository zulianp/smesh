#ifndef SMESH_DISTRIBUTED_CREATE_IMPL_HPP
#define SMESH_DISTRIBUTED_CREATE_IMPL_HPP

#include "smesh_alloc.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_create.hpp"
#include "smesh_tracer.hpp"

namespace smesh {

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
    ptrdiff_t *n_global_nodes_out, geom_t ***points_out) {
  SMESH_TRACE_SCOPE("hex8_cube_create_distributed");

  if (nx <= 0 || ny <= 0 || nz <= 0) {
    SMESH_ERROR("hex8_cube_create_distributed: invalid grid sizes "
                "(nx=%ld, ny=%ld, nz=%ld)\n",
                (long)nx, (long)ny, (long)nz);
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerHex = 8;
  constexpr int kSpatialDim = 3;

  const ptrdiff_t n_global_elements = nx * ny * nz;
  const ptrdiff_t nx1 = nx + 1;
  const ptrdiff_t ny1 = ny + 1;
  const ptrdiff_t nz1 = nz + 1;
  const ptrdiff_t n_global_nodes = nx1 * ny1 * nz1;

  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("hex8_cube_create_distributed: cube too small "
                "(n_global_elements=%ld, n_global_nodes=%ld, comm_size=%d)\n",
                (long)n_global_elements, (long)n_global_nodes, comm_size);
    return SMESH_FAILURE;
  }

  const ptrdiff_t e_start = rank_start(n_global_elements, comm_size, comm_rank);
  const ptrdiff_t n_local_elements =
      rank_split(n_global_elements, comm_size, comm_rank);
  const ptrdiff_t n_start = rank_start(n_global_nodes, comm_size, comm_rank);
  const ptrdiff_t n_local_nodes =
      rank_split(n_global_nodes, comm_size, comm_rank);

  // Hex8 corner offsets matching `mesh_hex8_cube_to_folder` ordering:
  //   v: (xi+ix[v], yi+iy[v], zi+iz[v])
  static const ptrdiff_t ix[kNodesPerHex] = {0, 1, 1, 0, 0, 1, 1, 0};
  static const ptrdiff_t iy[kNodesPerHex] = {0, 0, 1, 1, 0, 0, 1, 1};
  static const ptrdiff_t iz[kNodesPerHex] = {0, 0, 0, 0, 1, 1, 1, 1};

  // Allocate output buffers
  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerHex * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerHex; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }

  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  // Connectivity: file ordering is e = zi*(ny*nx) + yi*nx + xi.
  // Node id is iv = xi + yi*(nx+1) + zi*(ny+1)*(nx+1).
  const ptrdiff_t exy = nx * ny;
  const ptrdiff_t nxy1 = nx1 * ny1;
  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t ge = e_start + le;
    const ptrdiff_t zi = ge / exy;
    const ptrdiff_t rem = ge - zi * exy;
    const ptrdiff_t yi = rem / nx;
    const ptrdiff_t xi = rem - yi * nx;

    for (int v = 0; v < kNodesPerHex; ++v) {
      const ptrdiff_t node_id =
          (xi + ix[v]) + (yi + iy[v]) * nx1 + (zi + iz[v]) * nxy1;
      elems[v][le] = static_cast<idx_t>(node_id);
    }
  }

  // Coordinates: file ordering is node = xi + yi*(nx+1) + zi*(ny+1)*(nx+1).
  const geom_t hx = (xmax - xmin) / static_cast<geom_t>(nx);
  const geom_t hy = (ymax - ymin) / static_cast<geom_t>(ny);
  const geom_t hz = (zmax - zmin) / static_cast<geom_t>(nz);

  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    const ptrdiff_t zi = gn / nxy1;
    const ptrdiff_t rem = gn - zi * nxy1;
    const ptrdiff_t yi = rem / nx1;
    const ptrdiff_t xi = rem - yi * nx1;

    points[0][ln] = static_cast<geom_t>(xmin + static_cast<geom_t>(xi) * hx);
    points[1][ln] = static_cast<geom_t>(ymin + static_cast<geom_t>(yi) * hy);
    points[2][ln] = static_cast<geom_t>(zmin + static_cast<geom_t>(zi) * hz);
  }

  *nnodesxelem_out = kNodesPerHex;
  *n_local_elements_out = n_local_elements;
  *n_global_elements_out = n_global_elements;
  *elems_out = elems;

  *spatial_dim_out = kSpatialDim;
  *n_local_nodes_out = n_local_nodes;
  *n_global_nodes_out = n_global_nodes;
  *points_out = points;

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_DISTRIBUTED_CREATE_IMPL_HPP
