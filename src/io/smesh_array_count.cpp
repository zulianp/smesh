#include "smesh_array_count.hpp"

#include "smesh_file_extensions.hpp"
#include "smesh_smz_format.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_glob.hpp"
#include "smesh_types.hpp"

#include <cstdio>

namespace smesh {

namespace {
size_t value_bytes_from_extension(std::string_view ext, size_t raw_value_bytes) {
  if (ext == "raw") return raw_value_bytes;
  if (ext == "int8" || ext == "uint8") return 1;
  if (ext == "int16" || ext == "uint16" || ext == "float16") return 2;
  if (ext == "int32" || ext == "uint32" || ext == "float32") return 4;
  if (ext == "int64" || ext == "uint64" || ext == "float64") return 8;
  return 0;
}
} // namespace

ptrdiff_t array_count_from_file(const Path &path, size_t raw_value_bytes) {
  if (!path.exists()) {
    fprintf(stderr, "array_count_from_file: no file %s\n", path.c_str());
    return -1;
  }
  if (is_compressed_path(path)) {
    compression::SmzInfo info;
    if (compression::smz_read_info(path, &info) != SMESH_SUCCESS) {
      return -1;
    }
    return static_cast<ptrdiff_t>(info.plan.n_values);
  }
  const size_t vb = value_bytes_from_extension(typed_extension(std::string_view(path.to_string())),
                                               raw_value_bytes);
  if (vb == 0) {
    fprintf(stderr, "array_count_from_file: unknown value type for %s\n", path.c_str());
    return -1;
  }
  const size_t bytes = file_size(path);
  if (bytes % vb != 0) {
    fprintf(stderr, "array_count_from_file: size of %s (%zu B) is not a multiple of %zu\n",
            path.c_str(), bytes, vb);
    return -1;
  }
  return static_cast<ptrdiff_t>(bytes / vb);
}

} // namespace smesh
