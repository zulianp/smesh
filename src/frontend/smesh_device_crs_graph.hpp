#ifndef SMESH_DEVICE_CRS_GRAPH_HPP
#define SMESH_DEVICE_CRS_GRAPH_HPP

#include "smesh_crs_graph.hpp"
#include "smesh_device_buffer.hpp"

namespace smesh {

template <typename count_t, typename idx_t>
std::shared_ptr<CRSGraph<count_t, idx_t>>
to_device(const std::shared_ptr<CRSGraph<count_t, idx_t>> &in) {
  if (in->rowptr()->mem_space() == MEMORY_SPACE_DEVICE) {
    return in;
  }

  return std::make_shared<CRSGraph<count_t, idx_t>>(to_device(in->rowptr()),
                                                    to_device(in->colidx()));
}

} // namespace smesh

#endif // SMESH_DEVICE_CRS_GRAPH_HPP