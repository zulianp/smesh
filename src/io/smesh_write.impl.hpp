#ifndef SMESH_WRITE_IMPL_HPP
#define SMESH_WRITE_IMPL_HPP

#include "smesh_glob.hpp"
#include "smesh_alloc.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_compression_options.hpp"
#include "smesh_node_height.hpp"
#include "smesh_write.hpp"

#include <cstdio>
#include <string>
#include <string_view>
#include <vector>

namespace smesh {

static inline void remove_other_typed_files(const Path &folder,
                                            const std::string &stem,
                                            const Path &keep) {
  const std::string keep_s = keep.to_string();
  const std::vector<std::string> files =
      find_files((folder / (stem + ".*")).to_string());
  for (const auto &f : files) {
    if (f != keep_s) {
      std::remove(f.c_str());
    }
  }
}

template <typename T>
int array_write(const Path &path, const T *const SMESH_RESTRICT data,
                const ptrdiff_t n_elements) {
  FILE *fp = fopen(path.c_str(), "wb");
  if (!fp) {
    SMESH_ERROR("write_raw_array: Unable to write file %s\n", path.c_str());
    return SMESH_FAILURE;
  }

  const size_t n = static_cast<size_t>(n_elements);
  const bool ok = fwrite(data, sizeof(T), n, fp) == n;
  if (fclose(fp) != 0 || !ok) {
    fprintf(stderr, "array_write: failed to write %s\n", path.c_str());
    return SMESH_FAILURE;
  }
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
  FileType *buffer = (FileType *)SMESH_ALLOC(buffer_size * sizeof(FileType));

  bool ok = true;
  for (ptrdiff_t i = 0; ok && i < n_elements; i += buffer_size) {
    ptrdiff_t n = std::min(buffer_size, n_elements - i);
    for (ptrdiff_t j = 0; j < n; j++) {
      buffer[j] = static_cast<FileType>(data[i + j]);
    }
    ok = fwrite(buffer, sizeof(FileType), static_cast<size_t>(n), fp) == static_cast<size_t>(n);
  }

  SMESH_FREE(buffer);
  if (fclose(fp) != 0 || !ok) {
    fprintf(stderr, "array_write_convert: failed to write %s\n", path.c_str());
    return SMESH_FAILURE;
  }
  return SMESH_SUCCESS;
}

/// Write `data` to `path` in the type named by its extension. Compressed paths
/// (`<stem>.<dtype>.smz`) are written with `opts` (lossless if off).
template <typename T>
int array_write_convert_from_extension(
    const Path &path, const T *const SMESH_RESTRICT data,
    const ptrdiff_t n_elements, const compression::CompressionOptions &opts) {
  using namespace compression;
  if (!is_compressed_path(path)) {
    return array_write_convert_from_extension<T>(path, data, n_elements);
  }
  const std::string ext = typed_extension(path);
  if (ext == "raw") {
    if constexpr (is_smz_type_v<T>) {
      return smz_write<T>(path, data, n_elements, opts);
    }
  } else if (ext == "float16") {
    return smz_write_convert<f16, T>(path, data, n_elements, opts);
  } else if (ext == "float32") {
    return smz_write_convert<f32, T>(path, data, n_elements, opts);
  } else if (ext == "float64") {
    return smz_write_convert<f64, T>(path, data, n_elements, opts);
  } else if (ext == "int16") {
    return smz_write_convert<i16, T>(path, data, n_elements, opts);
  } else if (ext == "int32") {
    return smz_write_convert<i32, T>(path, data, n_elements, opts);
  } else if (ext == "int64") {
    return smz_write_convert<i64, T>(path, data, n_elements, opts);
  }
  SMESH_ERROR("Unsupported file extension %s for file %s\n", ext.c_str(),
              path.c_str());
  return SMESH_FAILURE;
}

template <typename T>
int array_write_convert_from_extension(const Path &path,
                                       const T *const SMESH_RESTRICT data,
                                       const ptrdiff_t n_elements) {
  if (is_compressed_path(path)) {
    return array_write_convert_from_extension<T>(
        path, data, n_elements,
        compression::CompressionOptions::lossless());
  }
  auto ext = path.extension();
  if (ext == "raw") {
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

/// Coordinates in bound mode h need per-node bounds eb * h_n, the smallest height of the elements
/// at each node (over all blocks). Fills `bound` and returns the options to use: the requested
/// ones, or lossless (with a warning) for element types without a height.
template <typename idx_t, typename geom_t>
static inline compression::CompressionOptions coordinate_bounds(
    const compression::CompressionOptions &requested, const size_t n_blocks,
    const enum ElemType *types, const ptrdiff_t *n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT *elements,
    const int spatial_dim, const ptrdiff_t n_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    std::vector<double> *bound) {
  if (!requested.enabled || requested.mode != compression::BoundMode::NodeHeight) {
    return requested;
  }
  bound->assign(size_t(n_nodes), INFINITY);
  for (size_t b = 0; b < n_blocks; b++) {
    // (drop the restrict qualifiers of the writer's signature)
    if (accumulate_node_heights<idx_t, geom_t>(
            types[b], n_elements[b], const_cast<const idx_t *const *>(elements[b]),
            spatial_dim, const_cast<const geom_t *const *>(points), bound->data()) !=
        SMESH_SUCCESS) {
      fprintf(stderr,
              "compression: no element height for %s: coordinates are written losslessly\n",
              type_to_string(types[b]));
      bound->clear();
      compression::CompressionOptions lossless = requested;
      lossless.mode = compression::BoundMode::Lossless;
      return lossless;
    }
  }
  for (double &v : *bound) v *= requested.eb;  // inf (orphan nodes) stays inf: stored exactly
  return requested;
}

/// Write through the extension dispatch, compressed when `opts` asks for it.
template <typename T>
static inline int array_write_maybe_compressed(
    const Path &path, const T *data, const ptrdiff_t n,
    const compression::CompressionOptions *opts) {
  return (opts && opts->enabled)
             ? array_write_convert_from_extension<T>(path, data, n, *opts)
             : array_write_convert_from_extension<T>(path, data, n);
}

template <typename idx_t>
int mesh_block_to_folder(
    const Path &folder, int n_nodes_x_elem, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const compression::CompressionOptions *opts) {
  int ret = SMESH_SUCCESS;
  for (int d = 0; d < n_nodes_x_elem; ++d) {
    std::string fname = std::string("i") + std::to_string(d) + "." +
                        std::string(TypeToString<idx_t>::value()) +
                        compression::compressed_suffix(opts);
    Path i_path = folder / Path(fname);
    if (array_write_maybe_compressed<idx_t>(i_path, elements[d], n_elements,
                                            opts) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
    remove_other_typed_files(folder, std::string("i") + std::to_string(d), i_path);
  }

  const auto extras =
      detect_files(folder / "i*.*", {"raw", "int16", "int32", "int64"});
  for (const auto &p : extras) {
    const std::string stem = typed_stem(p);
    int ii = -1;
    if (!parse_soa_index_stem(stem, &ii)) {
      continue;
    }
    if (ii >= n_nodes_x_elem) {
      std::remove(p.c_str());
    }
  }

  return ret;
}

template <typename geom_t>
int mesh_coordinates_to_folder(
    const Path &folder, int spatial_dim, const ptrdiff_t n_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const compression::CompressionOptions *opts, const double *bound) {

  static constexpr std::string_view xyz[3] = {"x", "y", "z"};
  if (create_directory(folder) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  int ret = SMESH_SUCCESS;
  for (int d = 0; d < spatial_dim; ++d) {
    std::string fname = std::string(xyz[d]) + "." +
                        std::string(TypeToString<geom_t>::value()) +
                        compression::compressed_suffix(opts);
    Path x_path = folder / Path(fname);
    const int werr =
        (opts && opts->enabled && bound)
            ? compression::smz_write<geom_t>(x_path, points[d], n_nodes, *opts, bound)
            : array_write_maybe_compressed<geom_t>(x_path, points[d], n_nodes, opts);
    if (werr != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
    remove_other_typed_files(folder, std::string(xyz[d]), x_path);
  }

  return ret;
}

template <typename idx_t, typename geom_t>
int mesh_to_folder(const Path &path, enum ElemType element_type,
                   const ptrdiff_t n_elements,
                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                   const int spatial_dim, const ptrdiff_t n_nodes,
                   const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                   enum GeomMap geom_map,
                   const compression::MeshWriteOptions *opts) {
  const compression::CompressionOptions *conn = opts ? &opts->connectivity : nullptr;
  std::vector<double> bound;
  compression::CompressionOptions xyz_opts;
  if (opts) {
    xyz_opts = coordinate_bounds<idx_t, geom_t>(opts->coordinates, 1, &element_type,
                                                &n_elements, &elements, spatial_dim,
                                                n_nodes, points, &bound);
  }
  const compression::CompressionOptions *xyz = opts ? &xyz_opts : nullptr;
  int n_nodes_x_elem = elem_num_nodes(element_type);
  if (mesh_block_to_folder(path, n_nodes_x_elem, n_elements, elements, conn) !=
      SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }
  if (mesh_coordinates_to_folder(path, spatial_dim, n_nodes, points, xyz,
                                 bound.empty() ? nullptr : bound.data()) !=
      SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  if (mesh_write_yaml_basic(path, element_type, n_elements, spatial_dim,
                            n_nodes, TypeToString<idx_t>::value(),
                            TypeToString<geom_t>::value(), geom_map,
                            compression::compressed_suffix(conn),
                            compression::compressed_suffix(xyz)) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  return SMESH_SUCCESS;
}

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
                              const std::vector<enum GeomMap> &geom_maps,
                              const compression::MeshWriteOptions *opts) {
  const compression::CompressionOptions *conn = opts ? &opts->connectivity : nullptr;
  std::vector<double> bound;
  compression::CompressionOptions xyz_opts;
  if (opts) {
    xyz_opts = coordinate_bounds<idx_t, geom_t>(
        opts->coordinates, element_types.size(), element_types.data(),
        n_elements.data(), elements, spatial_dim, n_nodes, points, &bound);
  }
  const compression::CompressionOptions *xyz = opts ? &xyz_opts : nullptr;

  SMESH_ASSERT(block_names.size() == element_types.size());
  SMESH_ASSERT(block_names.size() == n_elements.size());

  if (create_directory(path) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  Path block_path = path / "blocks";
  if (create_directory(block_path) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  int ret = mesh_coordinates_to_folder(path, spatial_dim, n_nodes, points, xyz,
                                       bound.empty() ? nullptr : bound.data());
  const int n_blocks = static_cast<int>(block_names.size());
  for (int b = 0; b < n_blocks; ++b) {
    Path output_path = block_path / block_names[b];
    if (create_directory(output_path) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    const int n_nodes_x_elem = elem_num_nodes(element_types[b]);
    if (mesh_block_to_folder(output_path, n_nodes_x_elem, n_elements[b],
                             elements[b], conn) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
  }

  if (mesh_multiblock_write_yaml(path, n_blocks, block_names, element_types,
                                 n_elements, spatial_dim, n_nodes,
                                 TypeToString<idx_t>::value(),
                                 TypeToString<geom_t>::value(), geom_maps,
                                 compression::compressed_suffix(conn),
                                 compression::compressed_suffix(xyz)) !=
      SMESH_SUCCESS) {
    ret = SMESH_FAILURE;
  }

  const auto leftover_root =
      detect_files(path / "i*.*", {"raw", "int16", "int32", "int64"});
  for (const auto &p : leftover_root) {
    std::remove(p.c_str());
  }

  return ret;
}

} // namespace smesh

#endif // SMESH_WRITE_IMPL_HPP
