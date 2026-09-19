#ifndef SMESH_BLOCKS_META_HPP
#define SMESH_BLOCKS_META_HPP

#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"

#include <string>
#include <vector>

namespace smesh {

/// Read `meta.yaml` block list from a serial multi-block mesh folder.
/// Returns false when the folder is legacy single-block (no `blocks:` entry).
/// Missing `geom_map` keys default to `ISOPARAMETRIC`.
bool read_blocks_meta(const Path &path, std::vector<std::string> &block_names,
                      std::vector<enum ElemType> &element_types,
                      std::vector<enum GeomMap> &geom_maps);

} // namespace smesh

#endif // SMESH_BLOCKS_META_HPP
