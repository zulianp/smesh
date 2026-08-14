#include "smesh_sstet4_graph.hpp"
#include "smesh_sstet4_graph.impl.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE(IDX_T)                                                    \
    template int sstet4_generate_elements<IDX_T>(const int,                                                \
                                                 const ptrdiff_t,                                          \
                                                 const ptrdiff_t,                                          \
                                                 const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,  \
                                                 IDX_T *const SMESH_RESTRICT *const       SMESH_RESTRICT,  \
                                                 ptrdiff_t *,                                              \
                                                 ptrdiff_t *);                                             \
    template int sstet4_generate_elements_blocks<IDX_T>(                                                   \
            const int,                                                                                     \
            const ptrdiff_t,                                                                               \
            const ptrdiff_t *const,                                                                        \
            const ptrdiff_t,                                                                               \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
            const element_idx_t *const *const,                                                             \
            const block_idx_t *const *const,                                                               \
            ptrdiff_t *,                                                                                   \
            ptrdiff_t *);                                                                                  \
    template int sstet4_hierarchical_renumbering<IDX_T>(const int,                                         \
                                                        const int,                                         \
                                                        int *const,                                        \
                                                        const ptrdiff_t,                                   \
                                                        const ptrdiff_t,                                   \
                                                        IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                        IDX_T *const                       SMESH_RESTRICT, \
                                                        const bool);                                       \
    template int sstet4_hierarchical_renumbering_blocks<IDX_T>(                                            \
            const int,                                                                                     \
            const int,                                                                                     \
            int *const,                                                                                    \
            const ptrdiff_t,                                                                               \
            const ptrdiff_t *const,                                                                        \
            const ptrdiff_t,                                                                               \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
            IDX_T *const SMESH_RESTRICT,                                                                   \
            const bool)

    SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE(i32);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE(i64);

#undef SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE

}  // namespace smesh

