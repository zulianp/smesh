#ifndef SMESH_DISTRIBUTED_CREATE_IMPL_HPP
#define SMESH_DISTRIBUTED_CREATE_IMPL_HPP

#include "smesh_alloc.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_create.hpp"
#include "smesh_hex8_nozzle_grid.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <math.h>
#include <vector>

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

  // Match serial `mesh_fill_hex8_cube`: double spacing, then cast to geom_t.
  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;
  const double hz = (zmax - zmin) * 1. / nz;

  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    const ptrdiff_t zi = gn / nxy1;
    const ptrdiff_t rem = gn - zi * nxy1;
    const ptrdiff_t yi = rem / nx1;
    const ptrdiff_t xi = rem - yi * nx1;

    points[0][ln] = (geom_t)(xmin + xi * hx);
    points[1][ln] = (geom_t)(ymin + yi * hy);
    points[2][ln] = (geom_t)(zmin + zi * hz);
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

template <typename idx_t, typename geom_t>
int tet4_cube_create_distributed(
    MPI_Comm comm, const ptrdiff_t nx, const ptrdiff_t ny, const ptrdiff_t nz,
    const geom_t xmin, const geom_t ymin, const geom_t zmin, const geom_t xmax,
    const geom_t ymax, const geom_t zmax, int *nnodesxelem_out,
    ptrdiff_t *n_local_elements_out, ptrdiff_t *n_global_elements_out,
    idx_t ***elems_out, int *spatial_dim_out, ptrdiff_t *n_local_nodes_out,
    ptrdiff_t *n_global_nodes_out, geom_t ***points_out) {
  SMESH_TRACE_SCOPE("tet4_cube_create_distributed");

  if (nx <= 0 || ny <= 0 || nz <= 0) {
    SMESH_ERROR("tet4_cube_create_distributed: invalid grid sizes "
                "(nx=%ld, ny=%ld, nz=%ld)\n",
                (long)nx, (long)ny, (long)nz);
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerTet = 4;
  constexpr int kSpatialDim = 3;

  const ptrdiff_t n_hex = nx * ny * nz;
  const ptrdiff_t n_global_elements = 12 * n_hex;
  const ptrdiff_t nx1 = nx + 1;
  const ptrdiff_t ny1 = ny + 1;
  const ptrdiff_t nz1 = nz + 1;
  const ptrdiff_t nnodes_vertices = nx1 * ny1 * nz1;
  const ptrdiff_t n_global_nodes = nnodes_vertices + n_hex;

  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("tet4_cube_create_distributed: cube too small "
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

  static const int face_nodes[6][4] = {{0, 1, 2, 3}, {4, 7, 6, 5},
                                       {0, 4, 5, 1}, {3, 2, 6, 7},
                                       {0, 3, 7, 4}, {1, 5, 6, 2}};
  static const ptrdiff_t ix[8] = {0, 1, 1, 0, 0, 1, 1, 0};
  static const ptrdiff_t iy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
  static const ptrdiff_t iz[8] = {0, 0, 0, 0, 1, 1, 1, 1};

  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerTet * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerTet; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }

  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  const ptrdiff_t exy = nx * ny;
  const ptrdiff_t nxy1 = nx1 * ny1;
  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t gt = e_start + le;
    const ptrdiff_t hex_e = gt / 12;
    const int tet_in_hex = (int)(gt - hex_e * 12);
    const int face = tet_in_hex / 2;
    const int sub = tet_in_hex - face * 2;

    const ptrdiff_t zi = hex_e / exy;
    const ptrdiff_t rem = hex_e - zi * exy;
    const ptrdiff_t yi = rem / nx;
    const ptrdiff_t xi = rem - yi * nx;

    idx_t cube_nodes[8];
    for (int v = 0; v < 8; ++v) {
      cube_nodes[v] = static_cast<idx_t>((xi + ix[v]) + (yi + iy[v]) * nx1 +
                                         (zi + iz[v]) * nxy1);
    }
    const idx_t center_idx =
        static_cast<idx_t>(nnodes_vertices + hex_e);
    const int *fn = face_nodes[face];
    if (sub == 0) {
      elems[0][le] = cube_nodes[fn[0]];
      elems[1][le] = cube_nodes[fn[1]];
      elems[2][le] = cube_nodes[fn[2]];
      elems[3][le] = center_idx;
    } else {
      elems[0][le] = cube_nodes[fn[0]];
      elems[1][le] = cube_nodes[fn[2]];
      elems[2][le] = cube_nodes[fn[3]];
      elems[3][le] = center_idx;
    }
  }

  const geom_t hx = (xmax - xmin) / static_cast<geom_t>(nx);
  const geom_t hy = (ymax - ymin) / static_cast<geom_t>(ny);
  const geom_t hz = (zmax - zmin) / static_cast<geom_t>(nz);

  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    if (gn < nnodes_vertices) {
      const ptrdiff_t zi = gn / nxy1;
      const ptrdiff_t rem = gn - zi * nxy1;
      const ptrdiff_t yi = rem / nx1;
      const ptrdiff_t xi = rem - yi * nx1;
      points[0][ln] = static_cast<geom_t>(xmin + static_cast<geom_t>(xi) * hx);
      points[1][ln] = static_cast<geom_t>(ymin + static_cast<geom_t>(yi) * hy);
      points[2][ln] = static_cast<geom_t>(zmin + static_cast<geom_t>(zi) * hz);
    } else {
      const ptrdiff_t hex_e = gn - nnodes_vertices;
      const ptrdiff_t zi = hex_e / exy;
      const ptrdiff_t rem = hex_e - zi * exy;
      const ptrdiff_t yi = rem / nx;
      const ptrdiff_t xi = rem - yi * nx;
      points[0][ln] =
          static_cast<geom_t>(xmin + (static_cast<geom_t>(xi) + geom_t(0.5)) * hx);
      points[1][ln] =
          static_cast<geom_t>(ymin + (static_cast<geom_t>(yi) + geom_t(0.5)) * hy);
      points[2][ln] =
          static_cast<geom_t>(zmin + (static_cast<geom_t>(zi) + geom_t(0.5)) * hz);
    }
  }

  *nnodesxelem_out = kNodesPerTet;
  *n_local_elements_out = n_local_elements;
  *n_global_elements_out = n_global_elements;
  *elems_out = elems;

  *spatial_dim_out = kSpatialDim;
  *n_local_nodes_out = n_local_nodes;
  *n_global_nodes_out = n_global_nodes;
  *points_out = points;

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int quad4_square_create_distributed(
    MPI_Comm comm, const ptrdiff_t nx, const ptrdiff_t ny, const geom_t xmin,
    const geom_t ymin, const geom_t xmax, const geom_t ymax,
    int *nnodesxelem_out, ptrdiff_t *n_local_elements_out,
    ptrdiff_t *n_global_elements_out, idx_t ***elems_out, int *spatial_dim_out,
    ptrdiff_t *n_local_nodes_out, ptrdiff_t *n_global_nodes_out,
    geom_t ***points_out) {
  SMESH_TRACE_SCOPE("quad4_square_create_distributed");

  if (nx <= 0 || ny <= 0) {
    SMESH_ERROR("quad4_square_create_distributed: invalid grid sizes "
                "(nx=%ld, ny=%ld)\n",
                (long)nx, (long)ny);
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerQuad = 4;
  constexpr int kSpatialDim = 2;

  const ptrdiff_t n_global_elements = nx * ny;
  const ptrdiff_t nx1 = nx + 1;
  const ptrdiff_t ny1 = ny + 1;
  const ptrdiff_t n_global_nodes = nx1 * ny1;

  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("quad4_square_create_distributed: square too small "
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

  static const ptrdiff_t ix[kNodesPerQuad] = {0, 1, 1, 0};
  static const ptrdiff_t iy[kNodesPerQuad] = {0, 0, 1, 1};

  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerQuad * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerQuad; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }

  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t ge = e_start + le;
    const ptrdiff_t yi = ge / nx;
    const ptrdiff_t xi = ge - yi * nx;
    for (int v = 0; v < kNodesPerQuad; ++v) {
      elems[v][le] =
          static_cast<idx_t>((xi + ix[v]) + (yi + iy[v]) * nx1);
    }
  }

  const geom_t hx = (xmax - xmin) / static_cast<geom_t>(nx);
  const geom_t hy = (ymax - ymin) / static_cast<geom_t>(ny);

  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    const ptrdiff_t yi = gn / nx1;
    const ptrdiff_t xi = gn - yi * nx1;
    points[0][ln] = static_cast<geom_t>(xmin + static_cast<geom_t>(xi) * hx);
    points[1][ln] = static_cast<geom_t>(ymin + static_cast<geom_t>(yi) * hy);
  }

  *nnodesxelem_out = kNodesPerQuad;
  *n_local_elements_out = n_local_elements;
  *n_global_elements_out = n_global_elements;
  *elems_out = elems;

  *spatial_dim_out = kSpatialDim;
  *n_local_nodes_out = n_local_nodes;
  *n_global_nodes_out = n_global_nodes;
  *points_out = points;

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int tri3_square_create_distributed(
    MPI_Comm comm, const ptrdiff_t nx, const ptrdiff_t ny, const geom_t xmin,
    const geom_t ymin, const geom_t xmax, const geom_t ymax,
    int *nnodesxelem_out, ptrdiff_t *n_local_elements_out,
    ptrdiff_t *n_global_elements_out, idx_t ***elems_out, int *spatial_dim_out,
    ptrdiff_t *n_local_nodes_out, ptrdiff_t *n_global_nodes_out,
    geom_t ***points_out) {
  SMESH_TRACE_SCOPE("tri3_square_create_distributed");

  if (nx <= 0 || ny <= 0) {
    SMESH_ERROR("tri3_square_create_distributed: invalid grid sizes "
                "(nx=%ld, ny=%ld)\n",
                (long)nx, (long)ny);
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerTri = 3;
  constexpr int kSpatialDim = 2;

  const ptrdiff_t n_global_elements = 2 * nx * ny;
  const ptrdiff_t nx1 = nx + 1;
  const ptrdiff_t ny1 = ny + 1;
  const ptrdiff_t n_global_nodes = nx1 * ny1;

  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("tri3_square_create_distributed: square too small "
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

  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerTri * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerTri; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }

  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t ge = e_start + le;
    const ptrdiff_t quad = ge >> 1;
    const ptrdiff_t yi = quad / nx;
    const ptrdiff_t xi = quad - yi * nx;
    const idx_t i0 = static_cast<idx_t>(xi + yi * nx1);
    const idx_t i1 = static_cast<idx_t>((xi + 1) + yi * nx1);
    const idx_t i2 = static_cast<idx_t>((xi + 1) + (yi + 1) * nx1);
    const idx_t i3 = static_cast<idx_t>(xi + (yi + 1) * nx1);
    if ((ge & 1) == 0) {
      elems[0][le] = i0;
      elems[1][le] = i1;
      elems[2][le] = i3;
    } else {
      elems[0][le] = i1;
      elems[1][le] = i2;
      elems[2][le] = i3;
    }
  }

  const geom_t hx = (xmax - xmin) / static_cast<geom_t>(nx);
  const geom_t hy = (ymax - ymin) / static_cast<geom_t>(ny);
  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    const ptrdiff_t yi = gn / nx1;
    const ptrdiff_t xi = gn - yi * nx1;
    points[0][ln] = static_cast<geom_t>(xmin + static_cast<geom_t>(xi) * hx);
    points[1][ln] = static_cast<geom_t>(ymin + static_cast<geom_t>(yi) * hy);
  }

  *nnodesxelem_out = kNodesPerTri;
  *n_local_elements_out = n_local_elements;
  *n_global_elements_out = n_global_elements;
  *elems_out = elems;

  *spatial_dim_out = kSpatialDim;
  *n_local_nodes_out = n_local_nodes;
  *n_global_nodes_out = n_global_nodes;
  *points_out = points;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int quad4_ring_create_distributed(
    MPI_Comm comm, const geom_t inner_radius, const geom_t outer_radius,
    const ptrdiff_t nlayers, const ptrdiff_t nelements,
    int *nnodesxelem_out, ptrdiff_t *n_local_elements_out,
    ptrdiff_t *n_global_elements_out, idx_t ***elems_out, int *spatial_dim_out,
    ptrdiff_t *n_local_nodes_out, ptrdiff_t *n_global_nodes_out,
    geom_t ***points_out) {
  SMESH_TRACE_SCOPE("quad4_ring_create_distributed");

  if (nlayers <= 0 || nelements <= 2 || inner_radius <= 0 ||
      outer_radius <= inner_radius) {
    SMESH_ERROR("quad4_ring_create_distributed: invalid parameters\n");
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerQuad = 4;
  constexpr int kSpatialDim = 3;

  const ptrdiff_t n_global_elements = nlayers * nelements;
  const ptrdiff_t n_global_nodes = (nlayers + 1) * nelements;
  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("quad4_ring_create_distributed: ring too small "
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

  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerQuad * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerQuad; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }

  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t ge = e_start + le;
    const ptrdiff_t layer = ge / nelements;
    const ptrdiff_t i = ge - layer * nelements;
    const idx_t n0 = static_cast<idx_t>(layer * nelements + ((i + 1) % nelements));
    const idx_t n1 = static_cast<idx_t>(layer * nelements + i);
    elems[0][le] = n0;
    elems[1][le] = n1;
    elems[2][le] = static_cast<idx_t>(n1 + nelements);
    elems[3][le] = static_cast<idx_t>(n0 + nelements);
  }

  const geom_t dangle = 2 * static_cast<geom_t>(M_PI) / static_cast<geom_t>(nelements);
  const geom_t dh = (outer_radius - inner_radius) / static_cast<geom_t>(nlayers);
  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    const ptrdiff_t layer = gn / nelements;
    const ptrdiff_t i = gn - layer * nelements;
    const geom_t r = inner_radius + dh * static_cast<geom_t>(layer);
    const geom_t a = dangle * static_cast<geom_t>(i);
    points[0][ln] = cos(a) * r;
    points[1][ln] = sin(a) * r;
    points[2][ln] = 0;
  }

  *nnodesxelem_out = kNodesPerQuad;
  *n_local_elements_out = n_local_elements;
  *n_global_elements_out = n_global_elements;
  *elems_out = elems;

  *spatial_dim_out = kSpatialDim;
  *n_local_nodes_out = n_local_nodes;
  *n_global_nodes_out = n_global_nodes;
  *points_out = points;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int hex8_lshape_create_distributed(
    MPI_Comm comm, const ptrdiff_t nx, const ptrdiff_t ny, const ptrdiff_t nz,
    const geom_t xmax, const geom_t ymax, const geom_t zmax,
    const ptrdiff_t nxs, const ptrdiff_t nys, int *nnodesxelem_out,
    ptrdiff_t *n_local_elements_out, ptrdiff_t *n_global_elements_out,
    idx_t ***elems_out, int *spatial_dim_out, ptrdiff_t *n_local_nodes_out,
    ptrdiff_t *n_global_nodes_out, geom_t ***points_out) {
  SMESH_TRACE_SCOPE("hex8_lshape_create_distributed");

  if (nx < 1 || ny < 1 || nz < 1 || nxs < 1 || nys < 1 || nxs >= nx ||
      nys >= ny || !(xmax > 0) || !(ymax > 0) || !(zmax > 0)) {
    SMESH_ERROR("hex8_lshape_create_distributed: invalid sizes "
                "(nx=%ld, ny=%ld, nz=%ld, nxs=%ld, nys=%ld)\n",
                (long)nx, (long)ny, (long)nz, (long)nxs, (long)nys);
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerHex = 8;
  constexpr int kSpatialDim = 3;
  static const ptrdiff_t ix[kNodesPerHex] = {0, 1, 1, 0, 0, 1, 1, 0};
  static const ptrdiff_t iy[kNodesPerHex] = {0, 0, 1, 1, 0, 0, 1, 1};
  static const ptrdiff_t iz[kNodesPerHex] = {0, 0, 0, 0, 1, 1, 1, 1};

  const ptrdiff_t n_kept_per_z = nx * ny - nxs * nys;
  const ptrdiff_t n_global_elements = n_kept_per_z * nz;
  const ptrdiff_t n_global_nodes =
      (nx + 1) * (ny + 1) * (nz + 1) - nxs * nys * (nz + 1);

  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("hex8_lshape_create_distributed: mesh too small "
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

  auto node_id = [&](const ptrdiff_t xi, const ptrdiff_t yi,
                     const ptrdiff_t zi) -> ptrdiff_t {
    const ptrdiff_t nx1 = nx + 1;
    const ptrdiff_t g = xi + yi * nx1 + zi * nx1 * (ny + 1);
    ptrdiff_t unused = zi * nxs * nys;
    if (yi < nys) {
      unused += yi * nxs + std::min(xi, nxs);
    } else {
      unused += nys * nxs;
    }
    return g - unused;
  };

  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerHex * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerHex; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }
  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  const ptrdiff_t n_kept_notch_rows = nys * (nx - nxs);
  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t ge = e_start + le;
    const ptrdiff_t zi = ge / n_kept_per_z;
    const ptrdiff_t rem = ge - zi * n_kept_per_z;
    ptrdiff_t xi, yi;
    if (rem < n_kept_notch_rows) {
      yi = rem / (nx - nxs);
      xi = nxs + rem - yi * (nx - nxs);
    } else {
      const ptrdiff_t rem2 = rem - n_kept_notch_rows;
      yi = nys + rem2 / nx;
      xi = rem2 - (yi - nys) * nx;
    }
    for (int v = 0; v < kNodesPerHex; ++v) {
      elems[v][le] = static_cast<idx_t>(
          node_id(xi + ix[v], yi + iy[v], zi + iz[v]));
    }
  }

  const double hx = (double)xmax / (double)nx;
  const double hy = (double)ymax / (double)ny;
  const double hz = (double)zmax / (double)nz;
  const ptrdiff_t n_nodes_per_z = (nx + 1) * (ny + 1) - nxs * nys;
  const ptrdiff_t n_kept_node_notch_rows = nys * (nx + 1 - nxs);
  const ptrdiff_t nx1 = nx + 1;
  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    const ptrdiff_t zi = gn / n_nodes_per_z;
    const ptrdiff_t rem = gn - zi * n_nodes_per_z;
    ptrdiff_t xi, yi;
    if (rem < n_kept_node_notch_rows) {
      yi = rem / (nx1 - nxs);
      xi = nxs + rem - yi * (nx1 - nxs);
    } else {
      const ptrdiff_t rem2 = rem - n_kept_node_notch_rows;
      yi = nys + rem2 / nx1;
      xi = rem2 - (yi - nys) * nx1;
    }
    points[0][ln] = (geom_t)(xi * hx);
    points[1][ln] = (geom_t)(yi * hy);
    points[2][ln] = (geom_t)(zi * hz);
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

template <typename idx_t, typename geom_t>
int hex8_nozzle_create_distributed(
    MPI_Comm comm, const geom_t *x_breaks, const geom_t *bore_radius,
    const ptrdiff_t *n_axial, const ptrdiff_t n_segments,
    const ptrdiff_t expansion, const geom_t expanded_radius,
    const ptrdiff_t n_core, const ptrdiff_t n_bore, const ptrdiff_t n_outer,
    const geom_t core_fraction, const geom_t radial_grading,
    const geom_t axial_grading, int *nnodesxelem_out,
    ptrdiff_t *n_local_elements_out, ptrdiff_t *n_global_elements_out,
    idx_t ***elems_out, int *spatial_dim_out, ptrdiff_t *n_local_nodes_out,
    ptrdiff_t *n_global_nodes_out, geom_t ***points_out) {
  SMESH_TRACE_SCOPE("hex8_nozzle_create_distributed");

  std::vector<geom_t> x_vec(x_breaks, x_breaks + n_segments + 1);
  std::vector<geom_t> r_vec(bore_radius, bore_radius + n_segments + 1);
  std::vector<ptrdiff_t> na_vec(n_axial, n_axial + n_segments);
  if (!nozzle_arguments_valid(x_vec, r_vec, na_vec, expansion, expanded_radius,
                              n_core, n_bore, n_outer, core_fraction,
                              radial_grading, axial_grading,
                              "hex8_nozzle_create_distributed")) {
    return SMESH_FAILURE;
  }

  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  constexpr int kNodesPerHex = 8;
  constexpr int kSpatialDim = 3;
  const bool has_expansion = expansion >= 0 && expansion < n_segments;
  const NozzleSection sec(n_core, n_bore, n_outer, has_expansion,
                          (double)core_fraction, (double)radial_grading);
  const NozzlePlanes pl(x_vec, r_vec, na_vec, (double)axial_grading);
  std::vector<NozzleQuad> quads;
  nozzle_build_quads(sec, quads);

  const ptrdiff_t n_inner_quads = (ptrdiff_t)sec.n * sec.n + sec.n_bore * sec.n_perim;
  const ptrdiff_t n_quads = (ptrdiff_t)quads.size();
  const ptrdiff_t k_exp = has_expansion
                              ? nozzle_first_downstream_cell(na_vec, expansion)
                              : 0;
  const ptrdiff_t n2d = sec.n2d;
  const ptrdiff_t n2d_inner = has_expansion ? sec.n2d_inner : n2d;
  const ptrdiff_t n_up_e = k_exp * n_inner_quads;
  const ptrdiff_t n_down_quads = n_quads;
  const ptrdiff_t n_global_elements =
      n_up_e + (pl.n_cells_x - k_exp) * n_down_quads;
  const ptrdiff_t n_global_nodes =
      k_exp * n2d_inner + (pl.n_planes - k_exp) * n2d;

  if (n_global_elements < comm_size || n_global_nodes < comm_size) {
    SMESH_ERROR("hex8_nozzle_create_distributed: mesh too small "
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

  auto compacted_node = [&](const ptrdiff_t k, const ptrdiff_t v) -> ptrdiff_t {
    if (k < k_exp) {
      return k * n2d_inner + v;
    }
    return k_exp * n2d_inner + (k - k_exp) * n2d + v;
  };

  idx_t **elems = (idx_t **)SMESH_ALLOC(kNodesPerHex * sizeof(idx_t *));
  for (int d = 0; d < kNodesPerHex; ++d) {
    elems[d] = (idx_t *)SMESH_ALLOC(n_local_elements * sizeof(idx_t));
  }
  geom_t **points = (geom_t **)SMESH_ALLOC(kSpatialDim * sizeof(geom_t *));
  for (int d = 0; d < kSpatialDim; ++d) {
    points[d] = (geom_t *)SMESH_ALLOC(n_local_nodes * sizeof(geom_t));
  }

  for (ptrdiff_t le = 0; le < n_local_elements; ++le) {
    const ptrdiff_t ge = e_start + le;
    ptrdiff_t k, q;
    if (ge < n_up_e) {
      k = ge / n_inner_quads;
      q = ge - k * n_inner_quads;
    } else {
      const ptrdiff_t rem = ge - n_up_e;
      k = k_exp + rem / n_down_quads;
      q = rem - (k - k_exp) * n_down_quads;
    }
    const NozzleQuad &quad = quads[(size_t)q];
    for (int a = 0; a < 4; ++a) {
      elems[a][le] = static_cast<idx_t>(compacted_node(k, quad.v[a]));
      elems[a + 4][le] = static_cast<idx_t>(compacted_node(k + 1, quad.v[a]));
    }
  }

  std::vector<double> by, bz, oy, oz;
  sec.reference(by, bz, oy, oz);
  const double R = has_expansion ? (double)expanded_radius : 0.0;
  const ptrdiff_t n_up_nodes = k_exp * n2d_inner;
  for (ptrdiff_t ln = 0; ln < n_local_nodes; ++ln) {
    const ptrdiff_t gn = n_start + ln;
    ptrdiff_t k, v;
    if (gn < n_up_nodes) {
      k = gn / n2d_inner;
      v = gn - k * n2d_inner;
    } else {
      const ptrdiff_t rem = gn - n_up_nodes;
      k = k_exp + rem / n2d;
      v = rem - (k - k_exp) * n2d;
    }
    const double rb = pl.plane_rb[(size_t)k];
    points[0][ln] = (geom_t)pl.plane_x[(size_t)k];
    points[1][ln] = (geom_t)(rb * by[(size_t)v] + (R - rb) * oy[(size_t)v]);
    points[2][ln] = (geom_t)(rb * bz[(size_t)v] + (R - rb) * oz[(size_t)v]);
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
