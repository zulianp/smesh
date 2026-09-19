#include "smesh_edgesets.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_EDGESETS(T)                                                     \
    template int extract_edges_from_edgeset<T, T>(                                                 \
            const enum ElemType,                                                                   \
            const T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                   \
            const ptrdiff_t,                                                                       \
            const T *const SMESH_RESTRICT,                                                         \
            const i16 *const SMESH_RESTRICT,                                                       \
            T *const SMESH_RESTRICT *const SMESH_RESTRICT);                                        \
    template int extract_nodeset_from_edgeset<T, T>(                                               \
            const enum ElemType,                                                                   \
            const T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                   \
            const ptrdiff_t,                                                                       \
            const T *const SMESH_RESTRICT,                                                         \
            const i16 *const SMESH_RESTRICT,                                                       \
            ptrdiff_t *,                                                                           \
            T **SMESH_RESTRICT);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_EDGESETS(i32);
SMESH_EXPLICIT_INSTANTIATE_EDGESETS(i64);
}  // namespace smesh
