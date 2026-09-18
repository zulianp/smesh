#ifndef SMESH_READ_IMPL_HPP
#define SMESH_READ_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_alloc.hpp"
#include "smesh_file_extensions.hpp"

#include <map>
#include <string_view>
#include <vector>

namespace smesh {

template <typename T>
int array_read(const Path &path, T **data, ptrdiff_t *n_elements) {
  FILE *fp = fopen(path.c_str(), "rb");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s\n", path.c_str());
    return SMESH_FAILURE;
  }
  fseek(fp, 0, SEEK_END);
  *n_elements = ftell(fp) / sizeof(T);
  fseek(fp, 0, SEEK_SET);
  *data = (T *)SMESH_ALLOC(*n_elements * sizeof(T));

  int ret = SMESH_SUCCESS;
  const size_t n = static_cast<size_t>(*n_elements);
  if (fread(*data, sizeof(T), n, fp) != n) {
    SMESH_FREE(*data);
    *data = nullptr;
    fprintf(stderr, "Failed to read file %s\n", path.c_str());
    ret = SMESH_FAILURE;
  }

  fclose(fp);
  return ret;
}

template <typename FileType, typename TargetType>
int array_read_convert(const Path &path, TargetType **data,
                       ptrdiff_t *n_elements) {
  if (std::is_same_v<FileType, TargetType>) {
    return array_read<TargetType>(path, data, n_elements);
  }

  FILE *fp = fopen(path.c_str(), "rb");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s\n", path.c_str());
    return SMESH_FAILURE;
  }
  fseek(fp, 0, SEEK_END);
  *n_elements = ftell(fp) / sizeof(FileType);
  fseek(fp, 0, SEEK_SET);

  *data = (TargetType *)SMESH_ALLOC(*n_elements * sizeof(TargetType));

  int ret = SMESH_SUCCESS;
  const size_t n = static_cast<size_t>(*n_elements);
  if (sizeof(FileType) <= sizeof(TargetType)) {
    FileType *temp = (FileType *)*data;
    if (fread(temp, sizeof(FileType), n, fp) != n) {
      fprintf(stderr, "Failed to read file %s\n", path.c_str());
      ret = SMESH_FAILURE;
    } else {
      for (ptrdiff_t i = *n_elements - 1; i >= 0; i--) {
        (*data)[i] = (TargetType)temp[i];
      }
    }
  } else {
    FileType *temp = (FileType *)SMESH_ALLOC(*n_elements * sizeof(FileType));
    if (fread(temp, sizeof(FileType), n, fp) != n) {
      fprintf(stderr, "Failed to read file %s\n", path.c_str());
      ret = SMESH_FAILURE;
    } else {
      for (ptrdiff_t i = 0; i < *n_elements; i++) {
        (*data)[i] = (TargetType)temp[i];
      }
    }
    SMESH_FREE(temp);
  }

  fclose(fp);
  if (ret == SMESH_FAILURE) {
    SMESH_FREE(*data);
    *data = nullptr;
  }

  return ret;
}

template <typename T>
int array_read_convert_from_extension(const Path &path, T **data,
                                      ptrdiff_t *n_elements) {
  auto ext = path.extension();
  if (ext == "raw") {
    // We trust the user that the raw file is of the correct type.
    return array_read<T>(path, data, n_elements);
  } else if (ext == "float16") {
    return array_read_convert<f16, T>(path, data, n_elements);
  } else if (ext == "float32") {
    return array_read_convert<f32, T>(path, data, n_elements);
  } else if (ext == "float64") {
    return array_read_convert<f64, T>(path, data, n_elements);
  } else if (ext == "int16") {
    return array_read_convert<i16, T>(path, data, n_elements);
  } else if (ext == "int32") {
    return array_read_convert<i32, T>(path, data, n_elements);
  } else if (ext == "int64") {
    return array_read_convert<i64, T>(path, data, n_elements);
  } else {
    SMESH_ERROR("Unsupported file extension %s for file %s\n", ext.c_str(),
                path.c_str());
    return SMESH_FAILURE;
  }
}

template <typename idx_t>
int mesh_block_from_folder(const Path &folder, int *nnodesxelem_out,
                           idx_t ***elems_out, ptrdiff_t *nelements_out) {
  ptrdiff_t n_elements = 0;

  std::vector<Path> i_files =
      detect_files(folder / "i*.*", {"raw", "int16", "int32", "int64"});

  std::map<int, std::vector<Path>> by_index;
  int max_ii = -1;
  for (const auto &i_path : i_files) {
    const std::string stem = i_path.file_name();
    int ii = 0;
    if (!parse_soa_index_stem(stem, &ii)) {
      continue;
    }
    by_index[ii].push_back(i_path);
    if (ii > max_ii) {
      max_ii = ii;
    }
  }

  const int nnodesxelem = max_ii + 1;
  if (nnodesxelem <= 0) {
    SMESH_ERROR("No connectivity files found in input folder %s\n", folder.c_str());
    *elems_out = nullptr;
    *nnodesxelem_out = 0;
    *nelements_out = 0;
    return SMESH_FAILURE;
  }

  idx_t **elems = (idx_t **)SMESH_CALLOC(nnodesxelem, sizeof(idx_t *));
  for (int d = 0; d < nnodesxelem; d++) {
    elems[d] = nullptr;
  }

  int ret = SMESH_SUCCESS;
  {
    ptrdiff_t n_elements0 = 0;
    for (int ii = 0; ii < nnodesxelem; ++ii) {
      const auto found = by_index.find(ii);
      if (found == by_index.end() || found->second.empty()) {
        SMESH_ERROR("Missing connectivity file i%d in %s\n", ii, folder.c_str());
        ret = SMESH_FAILURE;
        continue;
      }
      const std::string preferred =
          std::string("i") + std::to_string(ii) + "." +
          std::string(TypeToString<idx_t>::value());
      Path i_path = select_one_typed_file(found->second, preferred);

      idx_t *idx = 0;
      if (array_read_convert_from_extension<idx_t>(i_path, &idx, &n_elements) !=
          SMESH_SUCCESS) {
        SMESH_ERROR("Failed to read index file %s\n", i_path.c_str());
        ret = SMESH_FAILURE;
      }
      elems[ii] = idx;

      if (ii == 0) {
        n_elements0 = n_elements;
      } else {
        assert(n_elements0 == n_elements);

        if (n_elements0 != n_elements) {
          SMESH_ERROR("Inconsistent lenghts in input %ld != %ld\n",
                      (long)n_elements0, (long)n_elements);
          ret = SMESH_FAILURE;
        }
      }
    }
  }

  if (ret == SMESH_FAILURE) {
    for (int d = 0; d < nnodesxelem; d++) {
      SMESH_FREE(elems[d]);
    }
    SMESH_FREE(elems);
    *elems_out = nullptr;
    *nnodesxelem_out = 0;
    *nelements_out = 0;
    return SMESH_FAILURE;
  }

  *nnodesxelem_out = nnodesxelem;
  *nelements_out = n_elements;
  *elems_out = elems;

  return SMESH_SUCCESS;
}

template <typename geom_t>
int mesh_coordinates_from_folder(const Path &folder, int *spatial_dim_out,
                                 geom_t ***points_out, ptrdiff_t *nnodes_out) {
  ptrdiff_t n_nodes = 0;

  std::vector<Path> points_paths = select_coordinate_files(folder);
  const int ndims = static_cast<int>(points_paths.size());

  if (!ndims) {
    SMESH_ERROR("No coordinates found in input folder %s\n", folder.c_str());
    return SMESH_FAILURE;
  }

  geom_t **points = (geom_t **)SMESH_ALLOC(sizeof(geom_t *) * ndims);
  for (int d = 0; d < ndims; d++) {
    points[d] = 0;
  }

  int ret = SMESH_SUCCESS;
  for (int d = 0; d < ndims; ++d) {
    geom_t *points_d = 0;
    if (array_read_convert_from_extension(points_paths[d], &points_d,
                                          &n_nodes) != SMESH_SUCCESS) {
      ret = SMESH_FAILURE;
    }
    points[d] = points_d;
  }

  if (ret == SMESH_FAILURE) {
    for (int d = 0; d < ndims; d++) {
      SMESH_FREE(points[d]);
    }

    SMESH_FREE(points);
    *points_out = nullptr;
    *spatial_dim_out = 0;
    *nnodes_out = 0;
    return SMESH_FAILURE;
  }

  *spatial_dim_out = ndims;
  *nnodes_out = n_nodes;
  *points_out = points;
  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int mesh_from_folder(const Path &folder, int *nnodesxelem_out,
                     ptrdiff_t *nelements_out, idx_t ***elems_out,
                     int *spatial_dim_out, ptrdiff_t *nnodes_out,
                     geom_t ***points_out) {

  if (mesh_block_from_folder(folder, nnodesxelem_out, elems_out, nelements_out
                             ) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  if (mesh_coordinates_from_folder(folder, spatial_dim_out, points_out,
                                   nnodes_out) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_READ_IMPL_HPP