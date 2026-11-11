#ifndef SMESH_FILE_EXTENSIONS_HPP
#define SMESH_FILE_EXTENSIONS_HPP

#include "smesh_base.hpp"
#include "smesh_glob.hpp"
#include "smesh_path.hpp"
#include "smesh_types.hpp"

#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>

namespace smesh {

/// Compressed array files are named `<stem>.<dtype>.smz` (e.g. `x.float32.smz`, `i0.int32.smz`):
/// the dtype extension stays, and `.smz` is appended (see src/compression).
static constexpr std::string_view kCompressedExtension = "smz";

/// The extension that names the value type, for raw and compressed files alike:
/// "float32" for `x.float32` and for `x.float32.smz`.
static inline std::string_view typed_extension(std::string_view file) {
  const size_t slash = file.find_last_of("/\\");
  std::string_view name = slash == std::string_view::npos ? file : file.substr(slash + 1);
  size_t dot = name.find_last_of('.');
  if (dot == std::string_view::npos) {
    return {};
  }
  if (name.substr(dot + 1) == kCompressedExtension) {
    name = name.substr(0, dot);
    dot = name.find_last_of('.');
    if (dot == std::string_view::npos) {
      return {};
    }
  }
  return name.substr(dot + 1);
}

static inline bool is_compressed_path(std::string_view file) {
  const size_t dot = file.find_last_of('.');
  return dot != std::string_view::npos && file.substr(dot + 1) == kCompressedExtension;
}

static inline bool is_compressed_path(const Path &path) { return is_compressed_path(std::string_view(path.to_string())); }

static inline std::string typed_extension(const Path &path) {
  return std::string(typed_extension(std::string_view(path.to_string())));
}

/// The file name without its type extension(s): "x" for `x.float32` and for `x.float32.smz`.
static inline std::string typed_stem(const Path &path) {
  return is_compressed_path(path) ? path.without_extension().file_name() : path.file_name();
}

/// Files matching `pattern` whose type extension is one of `extensions`, raw or compressed.
static inline std::vector<Path>
detect_files(const Path &pattern,
             const std::initializer_list<std::string_view> extensions) {

  std::vector<std::string> files = find_files(pattern.to_string());
  std::vector<Path> paths;
  for (const auto &file : files) {
    if (std::find(extensions.begin(), extensions.end(), typed_extension(std::string_view(file))) !=
        extensions.end()) {
      paths.push_back(Path(file));
    }
  }
  return paths;
}

static inline PrimitiveType detect_real_type(std::string_view file) {
  return to_real_type(typed_extension(file));
}

static inline PrimitiveType detect_integer_type(std::string_view file) {
  return to_integer_type(typed_extension(file));
}

static inline std::string trim(std::string_view str) {
  size_t start = str.find_first_not_of(" \t\n\r\f\v");
  if (start == std::string_view::npos) {
    return "";
  }
  size_t end = str.find_last_not_of(" \t\n\r\f\v");
  return std::string(str.substr(start, end - start + 1));
}

static inline std::string path_basename(const Path &path) {
  const std::string &s = path.to_string();
  const size_t slash = s.find_last_of("/\\");
  return (slash == std::string::npos) ? s : s.substr(slash + 1);
}

static inline bool parse_soa_index_stem(const std::string &stem, int *ii) {
  if (stem.size() < 2 || stem[0] != 'i' || ii == nullptr) {
    return false;
  }
  char *end = nullptr;
  const long v = std::strtol(stem.c_str() + 1, &end, 10);
  if (end == stem.c_str() + 1 || *end != '\0' || v < 0 || v > 1024) {
    return false;
  }
  *ii = static_cast<int>(v);
  return true;
}

/// Prefer float64/int64 over 32-bit leftovers when several typed files share a stem.
static inline int typed_file_rank(std::string_view ext) {
  if (ext == "float64" || ext == "int64") {
    return 0;
  }
  if (ext == "float32" || ext == "int32") {
    return 1;
  }
  if (ext == "float16" || ext == "int16") {
    return 2;
  }
  return 3;
}

/// One file out of several that share a stem (`x.float32`, `x.float64`, `x.float32.smz`, ...):
/// `preferred_basename` if present (raw first, then compressed), otherwise the widest type, raw
/// before compressed.
static inline Path select_one_typed_file(std::vector<Path> files,
                                         const std::string &preferred_basename) {
  if (files.empty()) {
    return Path();
  }
  if (!preferred_basename.empty()) {
    const std::string compressed = preferred_basename + "." + std::string(kCompressedExtension);
    for (const std::string *name : {&preferred_basename, &compressed}) {
      for (const auto &p : files) {
        if (path_basename(p) == *name) {
          return p;
        }
      }
    }
  }
  std::sort(files.begin(), files.end(), [](const Path &a, const Path &b) {
    const int ra = typed_file_rank(typed_extension(std::string_view(a.to_string())));
    const int rb = typed_file_rank(typed_extension(std::string_view(b.to_string())));
    if (ra != rb) {
      return ra < rb;
    }
    const bool ca = is_compressed_path(a), cb = is_compressed_path(b);
    if (ca != cb) {
      return !ca;
    }
    return a.to_string() < b.to_string();
  });
  return files.front();
}

static inline Path select_axis_file(const Path &folder,
                                    const std::vector<Path> &globbed,
                                    const std::string &preferred_basename) {
  if (!preferred_basename.empty()) {
    const Path listed = folder / preferred_basename;
    if (listed.exists()) {
      return listed;
    }
  }
  return select_one_typed_file(globbed, preferred_basename);
}

/// Filenames listed under `points:` in meta.yaml (`x.float64`, ...).
static inline void read_coord_filenames_from_meta(const Path &folder,
                                                  std::string *x,
                                                  std::string *y,
                                                  std::string *z) {
  const Path meta = folder / "meta.yaml";
  if (!meta.exists()) {
    return;
  }
  std::ifstream ifs(meta.c_str());
  if (!ifs) {
    return;
  }
  bool in_points = false;
  std::string line;
  while (std::getline(ifs, line)) {
    const size_t hash = line.find('#');
    if (hash != std::string::npos) {
      line.resize(hash);
    }
    const std::string t = trim(line);
    if (t.empty()) {
      continue;
    }
    if (!in_points) {
      if (t.rfind("points:", 0) == 0) {
        in_points = true;
      }
      continue;
    }
    const bool indented =
        !line.empty() && (line[0] == ' ' || line[0] == '\t' || line[0] == '-');
    if (!indented) {
      break;
    }
    auto take = [&t](const char *key) -> std::string {
      const std::string k(key);
      if (t.size() < k.size() || t.compare(0, k.size(), k) != 0) {
        return {};
      }
      return trim(t.substr(k.size()));
    };
    if (std::string v = take("- x:"); !v.empty() && x) {
      *x = v;
    } else if (std::string v = take("- y:"); !v.empty() && y) {
      *y = v;
    } else if (std::string v = take("- z:"); !v.empty() && z) {
      *z = v;
    }
  }
}

static inline std::vector<Path> select_coordinate_files(const Path &folder) {
  std::string pref_x, pref_y, pref_z;
  read_coord_filenames_from_meta(folder, &pref_x, &pref_y, &pref_z);

  auto glob_axis = [&](const char *pattern) {
    return detect_files(folder / pattern, {"raw", "float16", "float32", "float64"});
  };

  Path x = select_axis_file(folder, glob_axis("x.*"), pref_x);
  Path y = select_axis_file(folder, glob_axis("y.*"), pref_y);
  Path z = select_axis_file(folder, glob_axis("z.*"), pref_z);
  if (x.empty()) {
    x = select_axis_file(folder, glob_axis("x0.*"), pref_x);
  }
  if (y.empty()) {
    y = select_axis_file(folder, glob_axis("x1.*"), pref_y);
  }
  if (z.empty()) {
    z = select_axis_file(folder, glob_axis("x2.*"), pref_z);
  }

  std::vector<Path> out;
  if (!x.empty()) {
    out.push_back(x);
  }
  if (!y.empty()) {
    out.push_back(y);
  }
  if (!z.empty()) {
    out.push_back(z);
  }
  return out;
}

} // namespace smesh

#endif // SMESH_FILE_EXTENSIONS_HPP
