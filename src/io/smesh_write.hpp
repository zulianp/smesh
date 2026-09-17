#ifndef SMESH_WRITE_HPP
#define SMESH_WRITE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"

#include <vector> 
#include <string>
#include <string_view>

#include "smesh_types.hpp"

namespace smesh {

int array_write(const Path &path, const enum PrimitiveType type, const void *const SMESH_RESTRICT data,
                const ptrdiff_t n_elements);

template <typename T>
int array_write(const Path &path, const T *const SMESH_RESTRICT data,
                const ptrdiff_t n_elements);

template <typename FileType, typename T>
int array_write_convert(const Path &path, const T *const SMESH_RESTRICT data,
                        const ptrdiff_t n_elements);

template <typename T>
int array_write_convert_from_extension(const Path &path,
                                       const T *const SMESH_RESTRICT data,
                                       const ptrdiff_t n_elements);

template <typename idx_t>
int mesh_block_to_folder(
    const Path &folder, int nnodesxelem, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements);

template <typename geom_t>
int mesh_coordinates_to_folder(
    const Path &folder, int spatial_dim, const ptrdiff_t n_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

/// `idx_type` / `geom_type` name the types actually written to disk (as
/// produced by `TypeToString`). They are separate arguments because the writers
/// are templated: the connectivity/coordinates on disk do not have to be the
/// compiled `idx_t`/`geom_t`, and the meta file must describe what was written.
int mesh_write_yaml_basic(const Path &path, enum ElemType element_type,
                          const ptrdiff_t n_elements, const int spatial_dim,
                          const ptrdiff_t n_nodes,
                          const std::string_view idx_type,
                          const std::string_view geom_type,
                          enum GeomMap geom_map = ISOPARAMETRIC);

/// Overload for writers that use the compiled `idx_t`/`geom_t`.
inline int mesh_write_yaml_basic(const Path &path, enum ElemType element_type,
                                 const ptrdiff_t n_elements,
                                 const int spatial_dim,
                                 const ptrdiff_t n_nodes,
                                 enum GeomMap geom_map = ISOPARAMETRIC) {
  return mesh_write_yaml_basic(path, element_type, n_elements, spatial_dim,
                               n_nodes, TypeToString<idx_t>::value(),
                               TypeToString<geom_t>::value(), geom_map);
}

template <typename idx_t, typename geom_t>
int mesh_to_folder(const Path &path, enum ElemType element_type,
                   const ptrdiff_t n_elements,
                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                   const int spatial_dim, const ptrdiff_t n_nodes,
                   const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                   enum GeomMap geom_map = ISOPARAMETRIC);

int mesh_multiblock_write_yaml(const Path &path, const uint16_t n_blocks,
                               const std::vector<std::string> &block_names,
                               const std::vector<enum ElemType> &element_types,
                               const std::vector<ptrdiff_t> &n_elements,
                               const int spatial_dim, const ptrdiff_t n_nodes,
                               const std::string_view idx_type,
                               const std::string_view geom_type,
                               const std::vector<enum GeomMap> &geom_maps = {});

/// Overload for writers that use the compiled `idx_t`/`geom_t`.
inline int mesh_multiblock_write_yaml(
    const Path &path, const uint16_t n_blocks,
    const std::vector<std::string> &block_names,
    const std::vector<enum ElemType> &element_types,
    const std::vector<ptrdiff_t> &n_elements, const int spatial_dim,
    const ptrdiff_t n_nodes,
    const std::vector<enum GeomMap> &geom_maps = {}) {
  return mesh_multiblock_write_yaml(path, n_blocks, block_names, element_types,
                                    n_elements, spatial_dim, n_nodes,
                                    TypeToString<idx_t>::value(),
                                    TypeToString<geom_t>::value(), geom_maps);
}

/// Stream SoA connectivity (`elements[d][e]`) to AoS file `e * nxe + d`.
int mesh_write_soa_to_aos(const Path &path, int n_nodes_x_elem, const ptrdiff_t n_elements,
                          const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements);

/// Same layout, reading existing `i0.<idx_t>`, `i1.<idx_t>`, … in `soa_folder`.
int mesh_write_soa_files_to_aos(const Path &soa_folder, const Path &aos_path, int n_nodes_x_elem,
                                const ptrdiff_t n_elements);

template <typename idx_t, typename geom_t>
int mesh_multiblock_to_folder(const Path &path,
                              const std::vector<std::string> &block_names,
                              const std::vector<enum ElemType> &element_types,
                              const std::vector<ptrdiff_t> &n_elements,
                              const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                                  elements[],
                              const int spatial_dim, const ptrdiff_t n_nodes,
                              const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                                  points,
                              const std::vector<enum GeomMap> &geom_maps = {});

} // namespace smesh

#endif // SMESH_WRITE_HPP
