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

#define SMESH_EXPLICIT_INSTANTIATE_SSTET4_SIDESET(IDX_T, ELEM_IDX_T)                                                             \
    template int sstet4_extract_surface_from_sideset<IDX_T, ELEM_IDX_T>(const int,                                               \
                                                                        const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                                        const ptrdiff_t,                                         \
                                                                        const ELEM_IDX_T *const SMESH_RESTRICT,                  \
                                                                        const i16 *const        SMESH_RESTRICT,                  \
                                                                        IDX_T **const           SMESH_RESTRICT);                           \
    template int sstet4_extract_nodeset_from_sideset<IDX_T, ELEM_IDX_T>(const int,                                               \
                                                                        const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                                        const ptrdiff_t,                                         \
                                                                        const ELEM_IDX_T *const SMESH_RESTRICT,                  \
                                                                        const i16 *const        SMESH_RESTRICT,                  \
                                                                        ptrdiff_t *,                                             \
                                                                        IDX_T **SMESH_RESTRICT)

#define SMESH_EXPLICIT_INSTANTIATE_SSTRI_REMAP(IDX_T)                                   \
    template int sstri_hierarchical_remapping<IDX_T>(const int,                         \
                                                     const int,                         \
                                                     int *const,                        \
                                                     const ptrdiff_t,                   \
                                                     const ptrdiff_t,                   \
                                                     IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                     IDX_T **,                          \
                                                     ptrdiff_t *)

    SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE(i32);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE(i64);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_SIDESET(i32, i32);
    SMESH_EXPLICIT_INSTANTIATE_SSTET4_SIDESET(i64, i64);
    SMESH_EXPLICIT_INSTANTIATE_SSTRI_REMAP(i32);
    SMESH_EXPLICIT_INSTANTIATE_SSTRI_REMAP(i64);

#undef SMESH_EXPLICIT_INSTANTIATE_SSTRI_REMAP
#undef SMESH_EXPLICIT_INSTANTIATE_SSTET4_SIDESET
#undef SMESH_EXPLICIT_INSTANTIATE_SSTET4_SINGLE

}  // namespace smesh

