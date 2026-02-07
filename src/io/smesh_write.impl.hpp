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
    const Path &folder, int n_nodes_x_elem, const ptrdiff_t n_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements) {
  int ret = SMESH_SUCCESS;
  for (int d = 0; d < n_nodes_x_elem; ++d) {
    std::string fname = std::string("i") + std::to_string(d) + "." +
                        std::string(TypeToString<idx_t>::value());
    Path i_path = folder / Path(fname);
    if (array_write_convert_from_extension<idx_t>(
            i_path, elements[d], n_elements) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
  }

  return ret;
}

template <typename geom_t>
int mesh_coordinates_to_folder(
    const Path &folder, int spatial_dim, const ptrdiff_t n_nodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {

  static constexpr std::string_view xyz[3] = {"x", "y", "z"};
  if (create_directory(folder) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  int ret = SMESH_SUCCESS;
  for (int d = 0; d < spatial_dim; ++d) {
    std::string fname =
        std::string(xyz[d]) + "." + std::string(TypeToString<geom_t>::value());
    Path x_path = folder / Path(fname);
    if (array_write_convert_from_extension<geom_t>(x_path, points[d],
                                                   n_nodes) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
  }

  return ret;
}

template <typename idx_t, typename geom_t>
int mesh_to_folder(const Path &path, enum ElemType element_type,
                   const ptrdiff_t n_elements,
                   idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                   const int spatial_dim, const ptrdiff_t n_nodes,
                   geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  int n_nodes_x_elem = elem_num_nodes(element_type);
  if (mesh_block_to_folder(path, n_nodes_x_elem, n_elements, elements) !=
      SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }
  if (mesh_coordinates_to_folder(path, spatial_dim, n_nodes, points) !=
      SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  if (mesh_write_yaml_basic(path, element_type, n_elements, spatial_dim,
                            n_nodes) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int mesh_multiblock_to_folder(const Path &path,
                              const std::vector<std::string_view> &block_names,
                              const std::vector<enum ElemType> &element_types,
                              const std::vector<ptrdiff_t> &n_elements,
                              idx_t **const elements[], const int spatial_dim,
                              const ptrdiff_t n_nodes, geom_t **const points) {

  SMESH_ASSERT(block_names.size() == element_types.size())
  SMESH_ASSERT(block_names.size() == n_elements.size());

  if (create_directory(path) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  Path block_path = path / "blocks";
  if (create_directory(block_path) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  int ret = mesh_coordinates_to_folder(path, spatial_dim, n_nodes, points);
  int n_blocks = block_names.size();
  for (int b = 0; b < n_blocks; ++b) {
    Path output_path = block_path / block_names[b];
    if (create_directory(output_path) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    if (mesh_block_to_folder(output_path, element_types[b], n_elements[b],
                             elements[b]) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
  }

  if (mesh_multiblock_write_yaml(path, n_blocks, block_names, element_types,
                                 n_elements, spatial_dim,
                                 n_nodes) != SMESH_SUCCESS) {
    ret = SMESH_FAILURE;
  }

  return ret;
}

} // namespace smesh

#endif // SMESH_WRITE_IMPL_HPP
