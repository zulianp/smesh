#ifndef SMESH_WRITE_HPP
#define SMESH_WRITE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"

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

int mesh_write_yaml_basic(const Path &path, enum ElemType element_type,
                          const ptrdiff_t n_elements, const int spatial_dim,
                          const ptrdiff_t n_nodes);

template <typename idx_t, typename geom_t>
int mesh_to_folder(const Path &path, enum ElemType element_type,
                   const ptrdiff_t n_elements,
                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                   const int spatial_dim, const ptrdiff_t n_nodes,
                   const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points);

int mesh_multiblock_write_yaml(const Path &path, const uint16_t n_blocks,
                               const std::vector<std::string> &block_names,
                               const std::vector<enum ElemType> &element_types,
                               const std::vector<ptrdiff_t> &n_elements,
                               const int spatial_dim, const ptrdiff_t n_nodes);

template <typename idx_t, typename geom_t>
int mesh_multiblock_to_folder(const Path &path,
                              const std::vector<std::string> &block_names,
                              const std::vector<enum ElemType> &element_types,
                              const std::vector<ptrdiff_t> &n_elements,
                              const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                                  elements[],
                              const int spatial_dim, const ptrdiff_t n_nodes,
                              const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT
                                  points);

} // namespace smesh

#endif // SMESH_WRITE_HPP
