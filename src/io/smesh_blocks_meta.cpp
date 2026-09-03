#include "smesh_blocks_meta.hpp"

#include <fstream>
#include <string>

#ifdef SMESH_ENABLE_RYAML
#include <ryml.hpp>
#include <ryml_std.hpp>
#endif

namespace smesh {

static std::string trim_copy(const std::string &s) {
  const size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) {
    return "";
  }
  const size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

/// Parse the `meta.yaml` written by `mesh_multiblock_write_yaml` (no RapidYAML).
static bool read_blocks_meta_generated(const Path &meta_file,
                                       std::vector<std::string> &block_names,
                                       std::vector<enum ElemType> &element_types) {
  std::ifstream ifs(meta_file.c_str());
  if (!ifs.good()) {
    return false;
  }

  bool in_blocks = false;
  std::string line;
  while (std::getline(ifs, line)) {
    const size_t hash = line.find('#');
    if (hash != std::string::npos) {
      line.resize(hash);
    }

    const std::string t = trim_copy(line);
    if (t.empty()) {
      continue;
    }

    if (!in_blocks) {
      if (t.rfind("blocks:", 0) == 0) {
        in_blocks = true;
      }
      continue;
    }

    const bool indented =
        !line.empty() && (line[0] == ' ' || line[0] == '\t' || line[0] == '-');
    if (!indented) {
      break;
    }

    if (t.rfind("- name:", 0) == 0) {
      const std::string name = trim_copy(t.substr(7));
      if (name.empty()) {
        continue;
      }
      block_names.push_back(name);
      element_types.push_back(INVALID);
    } else if (t.rfind("element_type:", 0) == 0 && !element_types.empty()) {
      const std::string et = trim_copy(t.substr(13));
      element_types.back() = type_from_string(et.c_str());
    }
  }

  return !block_names.empty();
}

bool read_blocks_meta(const Path &path, std::vector<std::string> &block_names,
                      std::vector<enum ElemType> &element_types) {
  block_names.clear();
  element_types.clear();

  auto meta_file = Path(path) / "meta.yaml";
  if (!meta_file.exists()) {
    return false;
  }

#if defined(SMESH_ENABLE_RYAML)
  std::ifstream ifs(meta_file.c_str(), std::ios::binary);
  if (!ifs.good()) {
    return false;
  }

  std::string yaml((std::istreambuf_iterator<char>(ifs)),
                   std::istreambuf_iterator<char>());
  if (yaml.empty()) {
    return false;
  }

  ryml::Tree tree = ryml::parse_in_arena(ryml::to_csubstr(yaml));
  auto root = tree.rootref();
  if (!root.has_child("blocks")) {
    return false;
  }

  auto blocks = root["blocks"];
  const size_t n = blocks.num_children();
  block_names.reserve(n);
  element_types.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    auto blk = blocks[i];
    if (!blk.has_child("name")) {
      continue;
    }

    auto name_val = blk["name"].val();
    std::string name(name_val.str, name_val.len);

    enum ElemType et = INVALID;
    if (blk.has_child("element_type")) {
      auto et_val = blk["element_type"].val();
      std::string et_str(et_val.str, et_val.len);
      et = type_from_string(et_str.c_str());
    }

    block_names.push_back(name);
    element_types.push_back(et);
  }

  return !block_names.empty();
#else
  return read_blocks_meta_generated(meta_file, block_names, element_types);
#endif
}

} // namespace smesh
