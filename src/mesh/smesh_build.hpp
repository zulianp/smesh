#ifndef SMESH_BUILD_HPP
#define SMESH_BUILD_HPP

#include "smesh_base.hpp"

/**
 * @brief Fill a structured HEX8 mesh for an axis-aligned box.
 *
 * Generates a regular \f$nx \times ny \times nz\f$ hexahedral mesh over
 * \f$[xmin,xmax]\times[ymin,ymax]\times[zmin,zmax]\f$.
 *
 * - **Elements (connectivity)** are written in SoA layout:
 *   `elements[local_node][elem]`, where `local_node` in \([0,7]\) and
 *   `elem` in \([0, nx*ny*nz)\).
 * - **Points (coordinates)** are written in SoA layout:
 *   `points[dim][node]`, where `dim` in \([0,2]\) corresponds to \((x,y,z)\) and
 *   `node` in \([0,(nx+1)(ny+1)(nz+1))\).
 *
 * **Node numbering** is x-fastest, then y, then z:
 * `node(xi,yi,zi) = xi + (nx+1)*yi + (nx+1)*(ny+1)*zi`.
 *
 * **HEX8 local node ordering** per element is:
 * \((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)\)
 * in \((xi,yi,zi)\) offsets from the element's \((xi,yi,zi)\) index.
 *
 * @tparam idx_t Index type for connectivity (e.g. `int32_t`, `int64_t`).
 * @tparam geom_t Coordinate scalar type (e.g. `float`, `double`).
 *
 * @param nx Number of elements in x direction (must be > 0).
 * @param ny Number of elements in y direction (must be > 0).
 * @param nz Number of elements in z direction (must be > 0).
 * @param xmin Minimum x coordinate.
 * @param ymin Minimum y coordinate.
 * @param zmin Minimum z coordinate.
 * @param xmax Maximum x coordinate (must be > xmin).
 * @param ymax Maximum y coordinate (must be > ymin).
 * @param zmax Maximum z coordinate (must be > zmin).
 * @param elements Pointer to an array of 8 pointers, each of length `nx*ny*nz`.
 * @param points Pointer to an array of 3 pointers, each of length
 *        `(nx+1)*(ny+1)*(nz+1)`.
 */
template <typename idx_t, typename geom_t>
void mesh_fill_hex8_cube(const int nx, const int ny, const int nz,
                         const geom_t xmin, const geom_t ymin,
                         const geom_t zmin, const geom_t xmax,
                         const geom_t ymax, const geom_t zmax,
                         idx_t *SMESH_RESTRICT *const SMESH_RESTRICT elements,
                         geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points);

/**
 * @brief Fill a structured TRI3 mesh for an axis-aligned rectangle.
 *
 * Generates a regular \f$nx \times ny\f$ quadrilateral grid over
 * \f$[xmin,xmax]\times[ymin,ymax]\f$, and splits each quad into 2 triangles:
 * \((i0,i1,i3)\) and \((i1,i2,i3)\), where
 * `i0=(0,0)`, `i1=(1,0)`, `i2=(1,1)`, `i3=(0,1)` in local quad coordinates.
 *
 * - **Elements (connectivity)** are written in SoA layout:
 *   `elements[local_node][elem]`, where `local_node` in \([0,2]\) and
 *   `elem` in \([0, 2*nx*ny)\).
 * - **Points (coordinates)** are written in SoA layout:
 *   `points[dim][node]`, where `dim` in \([0,1]\) corresponds to \((x,y)\) and
 *   `node` in \([0,(nx+1)(ny+1))\).
 *
 * **Node numbering** is x-fastest, then y:
 * `node(xi,yi) = xi + (nx+1)*yi`.
 *
 * @tparam idx_t Index type for connectivity (e.g. `int32_t`, `int64_t`).
 * @tparam geom_t Coordinate scalar type (e.g. `float`, `double`).
 *
 * @param nx Number of quads in x direction (must be > 0).
 * @param ny Number of quads in y direction (must be > 0).
 * @param xmin Minimum x coordinate.
 * @param ymin Minimum y coordinate.
 * @param xmax Maximum x coordinate (must be > xmin).
 * @param ymax Maximum y coordinate (must be > ymin).
 * @param elements Pointer to an array of 3 pointers, each of length `2*nx*ny`.
 * @param points Pointer to an array of 2 pointers, each of length
 *        `(nx+1)*(ny+1)`.
 */
template <typename idx_t, typename geom_t>
void mesh_fill_tri3_square(const int nx, const int ny, const geom_t xmin,
                           const geom_t ymin, const geom_t xmax,
                           const geom_t ymax,
                           idx_t *SMESH_RESTRICT *const SMESH_RESTRICT elements,
                           geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points);

/**
 * @brief Fill a structured QUAD4 mesh for an axis-aligned rectangle.
 *
 * Generates a regular \f$nx \times ny\f$ quad mesh over
 * \f$[xmin,xmax]\times[ymin,ymax]\f$.
 *
 * - **Elements (connectivity)** are written in SoA layout:
 *   `elements[local_node][elem]`, where `local_node` in \([0,3]\) and
 *   `elem` in \([0, nx*ny)\).
 * - **Points (coordinates)** are written in SoA layout:
 *   `points[dim][node]`, where `dim` in \([0,1]\) corresponds to \((x,y)\) and
 *   `node` in \([0,(nx+1)(ny+1))\).
 *
 * **Node numbering** is x-fastest, then y:
 * `node(xi,yi) = xi + (nx+1)*yi`.
 *
 * **QUAD4 local node ordering** per element is:
 * \((0,0),(1,0),(1,1),(0,1)\) in \((xi,yi)\) offsets from the element's
 * \((xi,yi)\) index.
 *
 * @tparam idx_t Index type for connectivity (e.g. `int32_t`, `int64_t`).
 * @tparam geom_t Coordinate scalar type (e.g. `float`, `double`).
 *
 * @param nx Number of elements in x direction (must be > 0).
 * @param ny Number of elements in y direction (must be > 0).
 * @param xmin Minimum x coordinate.
 * @param ymin Minimum y coordinate.
 * @param xmax Maximum x coordinate (must be > xmin).
 * @param ymax Maximum y coordinate (must be > ymin).
 * @param elements Pointer to an array of 4 pointers, each of length `nx*ny`.
 * @param points Pointer to an array of 2 pointers, each of length
 *        `(nx+1)*(ny+1)`.
 */
template <typename idx_t, typename geom_t>
void mesh_fill_quad4_square(
    const int nx, const int ny, const geom_t xmin, const geom_t ymin,
    const geom_t xmax, const geom_t ymax,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points);

#endif // SMESH_BUILD_HPP
