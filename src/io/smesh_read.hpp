#ifndef SMESH_READ_HPP
#define SMESH_READ_HPP

#include "smesh_base.hpp"
#include <string_view>

namespace smesh {

/**
 * @brief Read a binary array from disk.
 *
 * Reads the entire file at @p path as a contiguous array of elements of type
 * @p T.
 *
 * @tparam T Element type stored in the file.
 * @param[in]  path       Path to the input file.
 * @param[out] data       On success, set to a newly allocated buffer holding
 *                        @p *n_elements elements. Allocation uses `malloc()`;
 *                        the caller owns the buffer and must `free()` it.
 *                        Set to `nullptr` on failure.
 * @param[out] n_elements On success, number of elements read. Set to 0 on
 *                        failure.
 *
 * @return `SMESH_SUCCESS` on success, `SMESH_FAILURE` on error (e.g. file open
 *         or read failure).
 */
template <typename T>
int array_read(const Path &path, T **data, ptrdiff_t *n_elements);

/**
 * @brief Read a binary array from disk and convert element type.
 *
 * Reads the file at @p path as an array of @p FileType elements, converts each
 * value to @p TargetType and returns the converted array.
 *
 * If `std::is_same_v<FileType, TargetType>` this is equivalent to
 * `array_read<TargetType>()`.
 *
 * @tparam FileType   Element type stored in the file.
 * @tparam TargetType Element type to return.
 * @param[in]  path       Path to the input file.
 * @param[out] data       On success, set to a newly allocated buffer holding
 *                        @p *n_elements elements (allocated with `malloc()`).
 *                        Caller must `free()` it. Set to `nullptr` on failure.
 * @param[out] n_elements On success, number of elements read (in FileType, and
 *                        also number of returned TargetType elements). Set to 0
 *                        on failure.
 *
 * @return `SMESH_SUCCESS` on success, `SMESH_FAILURE` on error.
 */
template <typename FileType, typename TargetType>
int array_read_convert(const Path &path, TargetType **data,
                       ptrdiff_t *n_elements);

/**
 * @brief Read a binary array and convert based on file extension.
 *
 * Supported extensions:
 * - `.raw`: interpreted as type @p T (no conversion, trusted)
 * - `.float16`, `.float32`, `.float64`
 * - `.int16`, `.int32`, `.int64`
 *
 * Note: the conversion direction is from the on-disk type implied by the file
 * extension to @p T.
 *
 * @tparam T Element type to return.
 * @param[in]  path       Input file path.
 * @param[out] data       Newly allocated buffer on success (via `malloc()`),
 *                        `nullptr` on failure. Caller must `free()` it.
 * @param[out] n_elements Number of elements read/returned on success, 0 on
 *                        failure.
 *
 * @return `SMESH_SUCCESS` on success, `SMESH_FAILURE` on error (including
 *         unsupported extensions).
 */
template <typename T>
int array_read_convert_from_extension(const Path &path, T **data,
                                      ptrdiff_t *n_elements);

/**
 * @brief Read element connectivity (indices) from a mesh folder.
 *
 * The folder is expected to contain index files matching `i*.*` (e.g. `i0.int32`,
 * `i1.int32`, ...). The numeric suffix after `i` specifies the local node
 * position within an element, and defines the output array slot.
 *
 * Each index file must have identical length (number of elements).
 *
 * @tparam idx_t Index type to return (e.g. `i32`, `i64`).
 * @param[in]  folder         Mesh folder.
 * @param[out] nnodesxelem_out Number of node indices per element.
 * @param[out] nelements_out   Number of elements.
 * @param[out] elems_out       On success, `(*elems_out)[k]` points to the index
 *                             array for local node k, with length
 *                             @p *nelements_out. The outer array and each inner
 *                             array are allocated with `calloc()`/`malloc()`;
 *                             caller must `free()` each inner array and then
 *                             the outer pointer array.
 *
 * @return `SMESH_SUCCESS` on success, `SMESH_FAILURE` on error.
 */
template <typename idx_t>
int mesh_block_from_folder(const Path &folder, int *nnodesxelem_out,
                           ptrdiff_t *nelements_out, idx_t ***elems_out);

/**
 * @brief Read point coordinates from a mesh folder.
 *
 * The folder is expected to contain coordinate files in either of these naming
 * schemes (with a supported extension):
 * - `x.*`, `y.*`, `z.*`
 * - `x0.*`, `x1.*`, `x2.*`
 *
 * The spatial dimension is inferred from which of x/y/z are present.
 *
 * @tparam geom_t Coordinate scalar type to return (e.g. `f32`, `f64`).
 * @param[in]  folder          Mesh folder.
 * @param[out] spatial_dim_out Spatial dimension (1, 2, or 3).
 * @param[out] points_out      On success, an array of length
 *                             @p *spatial_dim_out where each entry points to a
 *                             coordinate array of length @p *nnodes_out.
 *                             All allocations use `malloc()`; caller must
 *                             `free()` each coordinate array and then the outer
 *                             pointer array.
 * @param[out] nnodes_out      Number of nodes.
 *
 * @return `SMESH_SUCCESS` on success, `SMESH_FAILURE` on error.
 */
template <typename geom_t>
int mesh_coordinates_from_folder(const Path &folder, int *spatial_dim_out,
                                 geom_t ***points_out, ptrdiff_t *nnodes_out);

/**
 * @brief Read connectivity and coordinates from a mesh folder.
 *
 * Convenience wrapper around `mesh_block_from_folder()` and
 * `mesh_coordinates_from_folder()`.
 *
 * @tparam idx_t  Index type to return.
 * @tparam geom_t Coordinate scalar type to return.
 * @param[in]  folder Mesh folder.
 * @param[out] nnodesxelem_out See `mesh_block_from_folder()`.
 * @param[out] nelements_out   See `mesh_block_from_folder()`.
 * @param[out] elems_out       See `mesh_block_from_folder()`.
 * @param[out] spatial_dim_out See `mesh_coordinates_from_folder()`.
 * @param[out] nnodes_out      See `mesh_coordinates_from_folder()`.
 * @param[out] points_out      See `mesh_coordinates_from_folder()`.
 *
 * @return `SMESH_SUCCESS` on success, `SMESH_FAILURE` on error.
 */
template <typename idx_t, typename geom_t>
int mesh_from_folder(const Path &folder, int *nnodesxelem_out,
                     ptrdiff_t *nelements_out, idx_t ***elems_out,
                     int *spatial_dim_out, ptrdiff_t *nnodes_out,
                     geom_t ***points_out);

} // namespace smesh

#endif // SMESH_READ_HPP