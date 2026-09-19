#ifndef SMESH_SSEDGE_RESTRICTION_HPP
#define SMESH_SSEDGE_RESTRICTION_HPP

#include <stddef.h>
#include "smesh_base.hpp"

namespace smesh {

    template <typename idx_t>
    int ssedge_element_node_incidence_count(const int                                               level,
                                            const int                                               stride,
                                            const ptrdiff_t                                         nelements,
                                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                            uint16_t *const SMESH_RESTRICT                          count);

    template <typename idx_t, typename T>
    int ssedge_restrict(const ptrdiff_t                                         nelements,
                        const int                                               from_level,
                        const int                                               from_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                        const uint16_t *const SMESH_RESTRICT                    from_element_to_node_incidence_count,
                        const int                                               to_level,
                        const int                                               to_level_stride,
                        const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                        const int                                               vec_size,
                        const T *const SMESH_RESTRICT                           from,
                        T *const SMESH_RESTRICT                                 to);

}  // namespace smesh

#endif  // SMESH_SSEDGE_RESTRICTION_HPP
