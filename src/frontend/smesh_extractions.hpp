#ifndef SMESH_EXTRACTIONS_HPP
#define SMESH_EXTRACTIONS_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_extract_shape_features.hpp"
#include "smesh_buffer.hpp"
#include "smesh_types.hpp"
#include "smesh_forward_declarations.hpp"

#include <utility>

namespace smesh {
    
SharedBuffer<idx_t *> extract_sharp_edges(
    Mesh &mesh,
    const geom_t cos_angle_threshold);

SharedBuffer<idx_t> extract_sharp_corners(
    const ptrdiff_t n_nodes,
    SharedBuffer<idx_t *> &sharp_edges,
    const bool edge_clean_up);

SharedBuffer<element_idx_t> extract_disconnected_faces(
    Mesh &mesh,
    SharedBuffer<idx_t *> &sharp_edges);
    
}   

#endif // SMESH_EXTRACTIONS_HPP
