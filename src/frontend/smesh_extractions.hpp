#ifndef SMESH_EXTRACTIONS_HPP
#define SMESH_EXTRACTIONS_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_types.hpp"

#include <memory>

namespace smesh {

std::shared_ptr<Edgeset> extract_sharp_edges(Mesh &mesh, const geom_t cos_angle_threshold);

std::shared_ptr<Nodeset> extract_sharp_corners(Mesh                       &mesh,
                                               std::shared_ptr<Edgeset>   &sharp_edges,
                                               const bool                  edge_clean_up);

SharedBuffer<element_idx_t> extract_disconnected_faces(Mesh &mesh, const Edgeset &sharp_edges);

}  // namespace smesh

#endif
