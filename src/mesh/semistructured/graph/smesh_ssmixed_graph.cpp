#include "smesh_ssmixed_graph.hpp"
#include "smesh_ssmixed_graph.impl.hpp"
#include "smesh_ssmixed_hex_dominant.impl.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_SSMIXED(IDX_T)                                                          \
    template int ssmixed_hex_tet_generate_elements_blocks<IDX_T>(                                          \
            const int,                                                                                     \
            const ptrdiff_t,                                                                               \
            const enum ElemType *const,                                                                    \
            const ptrdiff_t *const,                                                                        \
            const ptrdiff_t,                                                                               \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
            const element_idx_t *const *const,                                                             \
            const block_idx_t *const *const,                                                               \
            ptrdiff_t *,                                                                                   \
            ptrdiff_t *);                                                                                  \
    template int ssmixed_hex_tet_hierarchical_renumbering_blocks<IDX_T>(                                   \
            const int,                                                                                     \
            const int,                                                                                     \
            int *const,                                                                                    \
            const ptrdiff_t,                                                                               \
            const enum ElemType *const,                                                                    \
            const ptrdiff_t *const,                                                                        \
            const ptrdiff_t,                                                                               \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
            IDX_T *const SMESH_RESTRICT,                                                                   \
            const bool);                                                                                   \
    template int ssmixed_hex_dominant_generate_elements_blocks<IDX_T>(                                     \
            const int,                                                                                     \
            const ptrdiff_t,                                                                               \
            const enum ElemType *const,                                                                    \
            const ptrdiff_t *const,                                                                        \
            const ptrdiff_t,                                                                               \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                 \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
            const element_idx_t *const *const,                                                             \
            const block_idx_t *const *const,                                                               \
            ptrdiff_t *,                                                                                   \
            ptrdiff_t *);                                                                                  \
    template int ssmixed_hex_dominant_hierarchical_renumbering_blocks<IDX_T>(                              \
            const int,                                                                                     \
            const int,                                                                                     \
            int *const,                                                                                    \
            const ptrdiff_t,                                                                               \
            const enum ElemType *const,                                                                    \
            const ptrdiff_t *const,                                                                        \
            const ptrdiff_t,                                                                               \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                       \
            IDX_T *const SMESH_RESTRICT,                                                                   \
            const bool)

    SMESH_EXPLICIT_INSTANTIATE_SSMIXED(i32);
    SMESH_EXPLICIT_INSTANTIATE_SSMIXED(i64);

#undef SMESH_EXPLICIT_INSTANTIATE_SSMIXED

}  // namespace smesh
