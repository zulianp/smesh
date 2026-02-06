#ifndef SMESH_WRITE_IMPL_HPP
#define SMESH_WRITE_IMPL_HPP

#include "smesh_glob.hpp"
#include "smesh_write.hpp"

#include <string_view>
#include <vector>

namespace smesh {

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
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements) {
  SMESH_UNUSED(folder);
  SMESH_UNUSED(nnodesxelem);
  SMESH_UNUSED(nelements);
  SMESH_UNUSED(elements);
  return SMESH_FAILURE;
}

template <typename geom_t>
int mesh_coordinates_to_folder(
    const Path &folder, int spatial_dim, const ptrdiff_t n_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  SMESH_UNUSED(folder);
  SMESH_UNUSED(spatial_dim);
  SMESH_UNUSED(n_nodes);
  SMESH_UNUSED(points);
  return SMESH_FAILURE;
}

template <typename idx_t, typename geom_t>
int mesh_to_folder(const Path &path, enum ElemType element_type,
                   const ptrdiff_t n_elements,
                   idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                   const int spatial_dim, const ptrdiff_t n_nodes,
                   geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  SMESH_UNUSED(path);
  SMESH_UNUSED(element_type);
  SMESH_UNUSED(n_elements);
  SMESH_UNUSED(elements);
  SMESH_UNUSED(spatial_dim);
  SMESH_UNUSED(n_nodes);
  SMESH_UNUSED(points);
  return SMESH_FAILURE;
}

template <typename idx_t, typename geom_t>
int mesh_multiblock_to_folder(const std::vector<std::string_view> &block_names,
                              const std::vector<enum ElemType> &element_types,
                              const std::vector<ptrdiff_t> &n_elements,
                              idx_t **const elements[], const int spatial_dim,
                              const ptrdiff_t n_nodes, geom_t **const points) {
  SMESH_UNUSED(block_names);
  SMESH_UNUSED(element_types);
  SMESH_UNUSED(n_elements);
  SMESH_UNUSED(elements);
  SMESH_UNUSED(spatial_dim);
  SMESH_UNUSED(n_nodes);
  SMESH_UNUSED(points);
  return SMESH_FAILURE;
}

} // namespace smesh

#endif // SMESH_WRITE_IMPL_HPP
