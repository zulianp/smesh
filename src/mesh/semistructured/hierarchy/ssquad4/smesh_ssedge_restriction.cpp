#include "smesh_ssedge_restriction.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SSEDGE_RESTRICTION(IDX_T, T)                                                              \
    template int ssedge_restrict<IDX_T, T>(const ptrdiff_t                                         nelements,               \
                                           const int                                               from_level,              \
                                           const int                                               from_level_stride,       \
                                           const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,           \
                                           const uint16_t *const SMESH_RESTRICT                    from_element_to_node_incidence_count, \
                                           const int                                               to_level,                \
                                           const int                                               to_level_stride,         \
                                           const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,             \
                                           const int                                               vec_size,                \
                                           const T *const SMESH_RESTRICT                           from,                    \
                                           T *const SMESH_RESTRICT                                 to)

#define SMESH_EXPLICIT_INSTANTIATE_SSEDGE_ELEMENT_NODE_INCIDENCE_COUNT(IDX_T)                                          \
    template int ssedge_element_node_incidence_count<IDX_T>(const int                                               level,    \
                                                            const int                                               stride,   \
                                                            const ptrdiff_t                                         nelements, \
                                                            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements, \
                                                            uint16_t *const SMESH_RESTRICT                          count)

namespace smesh {
    SMESH_EXPLICIT_INSTANTIATE_SSEDGE_RESTRICTION(i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSEDGE_RESTRICTION(i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSEDGE_RESTRICTION(i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSEDGE_RESTRICTION(i64, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSEDGE_ELEMENT_NODE_INCIDENCE_COUNT(i32);
    SMESH_EXPLICIT_INSTANTIATE_SSEDGE_ELEMENT_NODE_INCIDENCE_COUNT(i64);
}  // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_SSEDGE_RESTRICTION
#undef SMESH_EXPLICIT_INSTANTIATE_SSEDGE_ELEMENT_NODE_INCIDENCE_COUNT
