#include "smesh_glob.hpp"

#include "smesh_base.hpp"

#ifndef _WIN32
#include <glob.h>
#include <sys/stat.h>
#else
#include "glob/glob.h"
#include <cassert>
#include <filesystem>
#include <iostream>
#endif

namespace smesh {

std::vector<std::string> find_files(const std::string_view &pattern) {
#ifndef _WIN32
  glob_t gl;
  glob(pattern.data(), GLOB_MARK, NULL, &gl);

  int n_files = gl.gl_pathc;
  std::vector<std::string> ret;
  for (int np = 0; np < n_files; np++) {
    ret.push_back(gl.gl_pathv[np]);
  }

  globfree(&gl);
  return ret;
#else
  std::vector<std::string> ret;
  for (auto path : glob::glob(pattern)) {
    ret.push_back(path.generic_string());
  }
  return ret;
#endif
}

size_t count_files(const std::string_view &pattern) {
  return find_files(pattern).size();
}

int create_directory(const std::string_view &path) {
#ifdef _WIN32
  namespace fs = std::filesystem;
  try {
    return fs::create_directory(path) ? SMESH_SUCCESS : SMESH_FAILURE;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return SMESH_FAILURE;
  }
#else
  struct stat st = {0};
  if (stat(path.data(), &st) == -1) {
    return mkdir(path.data(), 0700) == 0 ? SMESH_SUCCESS : SMESH_FAILURE;
  }

  return SMESH_SUCCESS;
#endif
}

} // namespace smesh
