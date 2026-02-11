#ifndef SMESH_READ_IMPL_HPP
#define SMESH_READ_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_file_extensions.hpp"

#include <string_view>

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
  *data = (T *)malloc(*n_elements * sizeof(T));

  int ret = SMESH_SUCCESS;
  const size_t n = static_cast<size_t>(*n_elements);
  if (fread(*data, sizeof(T), n, fp) != n) {
    free(*data);
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

  *data = (TargetType *)malloc(*n_elements * sizeof(TargetType));

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
    FileType *temp = (FileType *)malloc(*n_elements * sizeof(FileType));
    if (fread(temp, sizeof(FileType), n, fp) != n) {
      fprintf(stderr, "Failed to read file %s\n", path.c_str());
      ret = SMESH_FAILURE;
    } else {
      for (ptrdiff_t i = 0; i < *n_elements; i++) {
        (*data)[i] = (TargetType)temp[i];
      }
    }
    free(temp);
  }

  fclose(fp);
  if (ret == SMESH_FAILURE) {
    free(*data);
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
                           ptrdiff_t *nelements_out, idx_t ***elems_out) {
  ptrdiff_t n_elements = 0;

  std::vector<Path> i_files =
      detect_files(folder / "i*.*", {"raw", "int16", "int32", "int64"});

  int nnodesxelem = i_files.size();

  idx_t **elems = (idx_t **)calloc(nnodesxelem, sizeof(idx_t *));
  for (int d = 0; d < nnodesxelem; d++) {
    elems[d] = nullptr;
  }

  int ret = SMESH_SUCCESS;
  {
    ptrdiff_t n_elements0 = 0;
    for (int d = 0; d < nnodesxelem; ++d) {
      Path i_path = i_files[d];
      std::string filename = i_path.file_name();
      int ii = std::stoi(filename.substr(1, filename.find_last_of('.')));

      idx_t *idx = 0;
      if (array_read_convert_from_extension<idx_t>(i_path, &idx, &n_elements) !=
          SMESH_SUCCESS) {
        SMESH_ERROR("Failed to read index file %s\n", i_path.c_str());
        ret = SMESH_FAILURE;
      }
      // End of Selection
      elems[ii] = idx;

      if (d == 0) {
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
      free(elems[d]);
    }
    free(elems);
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

  std::vector<Path> x_file = detect_files(
      folder / "x.*", {"raw", "float16", "float32", "float64"});
  std::vector<Path> y_file = detect_files(
      folder / "y.*", {"raw", "float16", "float32", "float64"});
  std::vector<Path> z_file = detect_files(
      folder / "z.*", {"raw", "float16", "float32", "float64"});

  if (x_file.empty()) {
    x_file = detect_files(folder / "x0.*",
                          {"raw", "float16", "float32", "float64"});
  }

  if (y_file.empty()) {
    y_file = detect_files(folder / "x1.*",
                          {"raw", "float16", "float32", "float64"});
  }

  if (z_file.empty()) {
    z_file = detect_files(folder / "x2.*",
                          {"raw", "float16", "float32", "float64"});
  }

  int ndims = x_file.empty() ? 0 : 1; // x only
  ndims += y_file.empty() ? 0 : 1;    // x and y
  ndims += z_file.empty() ? 0 : 1;    // x, y and z

  if (!ndims) {
    SMESH_ERROR("No coordinates found in input folder %s\n", folder.c_str());
    return SMESH_FAILURE;
  }

  geom_t **points = (geom_t **)malloc(sizeof(geom_t *) * ndims);
  for (int d = 0; d < ndims; d++) {
    points[d] = 0;
  }

  std::vector<Path> points_paths = x_file;
  if (!y_file.empty()) {
    points_paths.push_back(y_file[0]);
  }
  if (!z_file.empty()) {
    points_paths.push_back(z_file[0]);
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
      free(points[d]);
    }

    free(points);
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

  if (mesh_block_from_folder(folder, nnodesxelem_out, nelements_out,
                             elems_out) != SMESH_SUCCESS) {
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