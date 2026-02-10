#ifndef SMESH_SSHEX8_RESTRICTION_HPP
#define SMESH_SSHEX8_RESTRICTION_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename real_t>
int sshex8_hierarchical_restriction(
    int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const real_t *const SMESH_RESTRICT from,
    real_t *const SMESH_RESTRICT to);

}

#endif // SMESH_SSHEX8_RESTRICTION_HPP
