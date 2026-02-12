#ifndef SMESH_EXTRACTIONS_HPP
#define SMESH_EXTRACTIONS_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_forward_declarations.hpp"

namespace smesh {
    
SharedBuffer<idx_t *> extract_sharp_edges(
    const &mesh,
    const geom_t angle_threshold);

// SharedBuffer<idx_t *> extract_sharp_corners(

}


#endif // SMESH_EXTRACTIONS_HPP
