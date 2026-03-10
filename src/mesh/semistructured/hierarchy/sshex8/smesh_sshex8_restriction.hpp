#ifndef SMESH_SSHEX8_RESTRICTION_HPP
#define SMESH_SSHEX8_RESTRICTION_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

    template <typename idx_t, typename real_t>
    int sshex8_hierarchical_restriction(int                                                     level,
                                        const ptrdiff_t                                         nelements,
                                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                        const u16 *const SMESH_RESTRICT                         element_to_node_incidence_count,
                                        const int                                               vec_size,
                                        const real_t *const SMESH_RESTRICT                      from,
                                        real_t *const SMESH_RESTRICT                            to);

    template <typename idx_t>
    int sshex8_element_node_incidence_count(const int                     level,
                                            const int                     stride,
                                            const ptrdiff_t               nelements,
                                            const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT   elements,
                                            uint16_t *const SMESH_RESTRICT count);

    template <typename idx_t, typename real_t>
    int sshex8_restrict(const ptrdiff_t                     nelements,
                        const int                           from_level,
                        const int                           from_level_stride,
                        const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT         from_elements,
                        const uint16_t *const SMESH_RESTRICT from_element_to_node_incidence_count,
                        const int                           to_level,
                        const int                           to_level_stride,
                        const idx_t *const SMESH_RESTRICT*const SMESH_RESTRICT         to_elements,
                        const int                           vec_size,
                        const real_t *const SMESH_RESTRICT   from,
                        real_t *const SMESH_RESTRICT         to);
}  // namespace smesh

#endif  // SMESH_SSHEX8_RESTRICTION_HPP
