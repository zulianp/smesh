#ifndef SMESH_FILE_EXTENSIONS_HPP
#define SMESH_FILE_EXTENSIONS_HPP

#include "smesh_base.hpp"
#include "smesh_glob.hpp"
#include "smesh_path.hpp"
#include "smesh_types.hpp"

#include <initializer_list>
#include <string>
#include <string_view>
#include <vector>

namespace smesh {

static inline std::vector<Path>
detect_files(const Path &pattern,
             const std::initializer_list<std::string_view> extensions) {

  std::vector<std::string> files = find_files(pattern.to_string());
  std::vector<Path> paths;
  for (auto file : files) {
    if (std::find(extensions.begin(), extensions.end(), file.substr(file.find_last_of('.') + 1)) !=
        extensions.end()) {
      paths.push_back(Path(file));
    }
  }
  return paths;
}

static inline PrimitiveType detect_real_type(std::string_view file) {
  std::string_view extension = file.substr(file.find_last_of('.') + 1);
  return to_real_type(extension);
}

static inline PrimitiveType detect_integer_type(std::string_view file) {
  std::string_view extension = file.substr(file.find_last_of('.') + 1);
  return to_integer_type(extension);
}

} // namespace smesh

#endif // SMESH_FILE_EXTENSIONS_HPP
