#ifndef SMESH_WRITE_IMPL_HPP
#define SMESH_WRITE_IMPL_HPP

#include "smesh_glob.hpp"
#include "smesh_write.hpp"

#include <string_view>
#include <vector>

namespace smesh {

template <typename T>
int array_write(const Path &path, const T *const SMESH_RESTRICT data,
                const ptrdiff_t n_elements) {
  FILE *fp = fopen(path.c_str(), "wb");
  if (!fp) {
    SMESH_ERROR("write_raw_array: Unable to write file %s\n", path.c_str());
    return SMESH_FAILURE;
  }

  fwrite(data, sizeof(T), n_elements, fp);
  fclose(fp);
  return SMESH_SUCCESS;
}

template <typename FileType, typename T>
int array_write_convert(const Path &path, const T *const SMESH_RESTRICT data,
                        const ptrdiff_t n_elements) {
  if (std::is_same_v<FileType, T>) {
    return array_write<T>(path, data, n_elements);
  }
  FILE *fp = fopen(path.c_str(), "wb");
  if (!fp) {
    SMESH_ERROR("write_raw_array: Unable to write file %s\n", path.c_str());
    return SMESH_FAILURE;
  }

  const ptrdiff_t buffer_size = std::min(ptrdiff_t(4096), n_elements);
  FileType *buffer = (FileType *)malloc(buffer_size * sizeof(FileType));

  for (ptrdiff_t i = 0; i < n_elements; i += buffer_size) {
    ptrdiff_t n = std::min(buffer_size, n_elements - i);
    for (ptrdiff_t j = 0; j < n; j++) {
      buffer[j] = static_cast<FileType>(data[i + j]);
    }
    fwrite(buffer, sizeof(FileType), n, fp);
  }

  free(buffer);
  fclose(fp);
  return SMESH_SUCCESS;
}

template <typename T>
int array_write_convert_from_extension(const Path &path,
                                       const T *const SMESH_RESTRICT data,
                                       const ptrdiff_t n_elements) {
  auto ext = path.extension();
  if (ext == ".raw") {
    return array_write<T>(path, data, n_elements);
  } else if (ext == "float16") {
    return array_write_convert<f16, T>(path, data, n_elements);
  } else if (ext == "float32") {
    return array_write_convert<f32, T>(path, data, n_elements);
  } else if (ext == "float64") {
    return array_write_convert<f64, T>(path, data, n_elements);
  } else if (ext == "int16") {
    return array_write_convert<i16, T>(path, data, n_elements);
  } else if (ext == "int32") {
    return array_write_convert<i32, T>(path, data, n_elements);
  } else if (ext == "int64") {
    return array_write_convert<i64, T>(path, data, n_elements);
  } else {
    SMESH_ERROR("Unsupported file extension %s for file %s\n", ext.c_str(),
                path.c_str());
    return SMESH_FAILURE;
  }
}

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
