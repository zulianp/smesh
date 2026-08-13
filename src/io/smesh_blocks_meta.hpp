#ifndef SMESH_BLOCKS_META_HPP
#define SMESH_BLOCKS_META_HPP

#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"

#include <string>
#include <vector>

namespace smesh {

/// Read `meta.yaml` block list from a serial multi-block mesh folder.
/// Returns false when the folder is legacy single-block (no `blocks:` entry).
bool read_blocks_meta(const Path &path, std::vector<std::string> &block_names,
                      std::vector<enum ElemType> &element_types);

} // namespace smesh

#endif // SMESH_BLOCKS_META_HPP
