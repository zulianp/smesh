/**
 * @file smesh_build.hpp
 * @brief Structured mesh generation helpers.
 */
#ifndef SMESH_BUILD_HPP
#define SMESH_BUILD_HPP

#include "smesh_base.hpp"

namespace smesh {

    template <typename idx_t, typename geom_t>
    void mesh_fill_hex8_reference_cube(
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
        geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

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
 *   `points[dim][node]`, where `dim` in \([0,2]\) corresponds to \((x,y,z)\)
 * and `node` in \([0,(nx+1)(ny+1)(nz+1))\).
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
void mesh_fill_hex8_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

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
void mesh_fill_tri3_square(
    const int nx, const int ny, const geom_t xmin, const geom_t ymin,
    const geom_t xmax, const geom_t ymax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

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
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

/**
 * @brief Fill a structured TET4 mesh for an axis-aligned box.
 *
 * Generates a regular \f$nx \times ny \times nz\f$ HEX8 grid over
 * \f$[xmin,xmax]\times[ymin,ymax]\times[zmin,zmax]\f$ and subdivides each HEX8
 * cell into **12 tetrahedra** by introducing a single **cell-center node**.
 *
 * - **Elements (connectivity)** are written in SoA layout:
 *   `elements[local_node][tet]`, where `local_node` in \([0,3]\) and
 *   `tet` in \([0, 12*nx*ny*nz)\).
 * - **Points (coordinates)** are written in SoA layout:
 *   `points[dim][node]`, where `dim` in \([0,2]\) corresponds to \((x,y,z)\)
 * and `node` in \([0, (nx+1)(ny+1)(nz+1) + nx*ny*nz)\).
 *
 * **Vertex node numbering** is x-fastest, then y, then z:
 * `node(xi,yi,zi) = xi + (nx+1)*yi + (nx+1)*(ny+1)*zi`.
 *
 * **Cell-center nodes** are appended after all vertex nodes, one per HEX8 cell,
 * in the same element traversal order:
 * `center(elem) = (nx+1)(ny+1)(nz+1) + elem`, where
 * `elem = xi + nx*yi + nx*ny*zi`.
 *
 * **TET4 generation per cell**: for each of the 6 quad faces, the face is split
 * into 2 triangles \((0,1,2)\) and \((0,2,3)\) (in that face’s local node
 * ordering), and each triangle forms a tetrahedron with the cell center as the
 * 4th node. This yields 12 tetrahedra per cell.
 *
 * @tparam idx_t Index type for connectivity (e.g. `int32_t`, `int64_t`).
 * @tparam geom_t Coordinate scalar type (e.g. `float`, `double`).
 *
 * @param nx Number of HEX cells in x direction (must be > 0).
 * @param ny Number of HEX cells in y direction (must be > 0).
 * @param nz Number of HEX cells in z direction (must be > 0).
 * @param xmin Minimum x coordinate.
 * @param ymin Minimum y coordinate.
 * @param zmin Minimum z coordinate.
 * @param xmax Maximum x coordinate (must be > xmin).
 * @param ymax Maximum y coordinate (must be > ymin).
 * @param zmax Maximum z coordinate (must be > zmin).
 * @param elements Pointer to an array of 4 pointers, each of length
 *        `12*nx*ny*nz`.
 * @param points Pointer to an array of 3 pointers, each of length
 *        `(nx+1)*(ny+1)*(nz+1) + nx*ny*nz`.
 */
template <typename idx_t, typename geom_t>
void mesh_fill_tet4_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

/**
 * @brief Fill a structured HEX8 mesh for an axis-aligned box, split into two
 * subdomains along a chosen axis.
 *
 * Builds the same vertex grid and coordinates as `mesh_fill_hex8_cube`, but
 * partitions the HEX8 elements into **left** and **right** sets according to:
 * `idx3[dim_split] < split_index`, where `idx3 = (xi, yi, zi)` are the element
 * indices with x-fastest traversal.
 *
 * - **Connectivity** is written in SoA layout:
 *   `left_elements[local_node][elem_left]` and
 *   `right_elements[local_node][elem_right]`, where `local_node` in \([0,7]\).
 * - **Points (coordinates)** are written in SoA layout:
 *   `points[dim][node]`, where `dim` in \([0,2]\) corresponds to \((x,y,z)\).
 *
 * **Vertex node numbering** is x-fastest, then y, then z:
 * `node(xi,yi,zi) = xi + (nx+1)*yi + (nx+1)*(ny+1)*zi`.
 *
 * **Required capacities**:
 * - `points[d]` length: `(nx+1)*(ny+1)*(nz+1)`
 * - `left_elements[l]` length: `split_index * ny * nz` if `dim_split==0`,
 *   `nx * split_index * nz` if `dim_split==1`, or `nx * ny * split_index` if
 *   `dim_split==2`
 * - `right_elements[l]` length: `nx*ny*nz - left_count`
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
 * @param dim_split Split axis: 0=x, 1=y, 2=z.
 * @param split_index Element-index threshold along `dim_split` in \([0,n]\).
 * @param left_elements Pointer to an array of 8 pointers for left connectivity.
 * @param right_elements Pointer to an array of 8 pointers for right
 * connectivity.
 * @param points Pointer to an array of 3 pointers for coordinates.
 */
template <typename idx_t, typename geom_t>
void mesh_fill_hex8_bidomain_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax, const int dim_split, const idx_t split_index,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT left_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT right_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

/**
 * @brief Fill a structured HEX8 mesh for an axis-aligned box, split into a
 * 3D checkerboard (parity) partition.
 *
 * Builds the same vertex grid and coordinates as `mesh_fill_hex8_cube`, but
 * partitions the HEX8 elements into **white** and **black** sets by parity of
 * the element indices:
 * - white if \((xi + yi + zi)\ \%\ 2 == 0\)
 * - black otherwise
 *
 * This requires `nx`, `ny`, and `nz` to be **even**, so that each set contains
 * exactly half of the elements.
 *
 * - **Connectivity** is written in SoA layout:
 *   `white_elements[local_node][elem_white]` and
 *   `black_elements[local_node][elem_black]`, where `local_node` in \([0,7]\).
 * - **Points (coordinates)** are written in SoA layout:
 *   `points[dim][node]`, where `dim` in \([0,2]\) corresponds to \((x,y,z)\).
 *
 * **Vertex node numbering** is x-fastest, then y, then z:
 * `node(xi,yi,zi) = xi + (nx+1)*yi + (nx+1)*(ny+1)*zi`.
 *
 * **Required capacities**:
 * - `points[d]` length: `(nx+1)*(ny+1)*(nz+1)`
 * - `white_elements[l]` length: `(nx*ny*nz)/2`
 * - `black_elements[l]` length: `(nx*ny*nz)/2`
 *
 * @tparam idx_t Index type for connectivity (e.g. `int32_t`, `int64_t`).
 * @tparam geom_t Coordinate scalar type (e.g. `float`, `double`).
 *
 * @param nx Number of elements in x direction (must be > 0, and even).
 * @param ny Number of elements in y direction (must be > 0, and even).
 * @param nz Number of elements in z direction (must be > 0, and even).
 * @param xmin Minimum x coordinate.
 * @param ymin Minimum y coordinate.
 * @param zmin Minimum z coordinate.
 * @param xmax Maximum x coordinate (must be > xmin).
 * @param ymax Maximum y coordinate (must be > ymin).
 * @param zmax Maximum z coordinate (must be > zmin).
 * @param white_elements Pointer to an array of 8 pointers for white
 * connectivity.
 * @param black_elements Pointer to an array of 8 pointers for black
 * connectivity.
 * @param points Pointer to an array of 3 pointers for coordinates.
 */
template <typename idx_t, typename geom_t>
void mesh_fill_hex8_checkerboard_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT black_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT white_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

template <typename idx_t, typename geom_t>
void mesh_fill_quad4_ring(
    const geom_t inner_radius, const geom_t outer_radius,
    const ptrdiff_t nlayers, const ptrdiff_t nelements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);
} // namespace smesh

#endif // SMESH_BUILD_HPP
