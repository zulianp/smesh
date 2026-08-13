#include "smesh_blocks_meta.hpp"

#include <fstream>
#include <string>

#ifdef SMESH_ENABLE_RYAML
#include <ryml.hpp>
#include <ryml_std.hpp>
#endif

namespace smesh {

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
  (void)path;
  return false;
#endif
}

} // namespace smesh
