#include "smesh_sspyramid_graph.hpp"
#include "smesh_sspyramid_graph.impl.hpp"
#include "smesh_types.hpp"

namespace smesh {

#define SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SINGLE(IDX_T)                                                   \
    template int sspyramid_generate_elements<IDX_T>(const int,                                               \
                                                    const ptrdiff_t,                                         \
                                                    const ptrdiff_t,                                         \
                                                    const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                    IDX_T *const SMESH_RESTRICT *const       SMESH_RESTRICT, \
                                                    ptrdiff_t *,                                             \
                                                    ptrdiff_t *);                                            \
    template int sspyramid_generate_elements_blocks<IDX_T>(                                                  \
            const int,                                                                                       \
            const ptrdiff_t,                                                                                 \
            const ptrdiff_t *const,                                                                          \
            const ptrdiff_t,                                                                                 \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                   \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                         \
            const element_idx_t *const *const,                                                               \
            const block_idx_t *const *const,                                                                 \
            ptrdiff_t *,                                                                                     \
            ptrdiff_t *);                                                                                    \
    template int sspyramid_hierarchical_renumbering<IDX_T>(const int,                                        \
                                                           const int,                                        \
                                                           int *const,                                       \
                                                           const ptrdiff_t,                                  \
                                                           const ptrdiff_t,                                  \
                                                           IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                           IDX_T *const                       SMESH_RESTRICT, \
                                                           const bool);                                      \
    template int sspyramid_hierarchical_renumbering_blocks<IDX_T>(                                           \
            const int,                                                                                       \
            const int,                                                                                       \
            int *const,                                                                                      \
            const ptrdiff_t,                                                                                 \
            const ptrdiff_t *const,                                                                          \
            const ptrdiff_t,                                                                                 \
            IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT *const SMESH_RESTRICT,                         \
            IDX_T *const SMESH_RESTRICT,                                                                     \
            const bool)

#define SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SIDESET(IDX_T, ELEM_IDX_T)                                                            \
    template int sspyramid_extract_surface_from_sideset<IDX_T, ELEM_IDX_T>(const int,                                              \
                                                                           const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                                           const ptrdiff_t,                                        \
                                                                           const ELEM_IDX_T *const SMESH_RESTRICT,                 \
                                                                           const i16 *const        SMESH_RESTRICT,                 \
                                                                           IDX_T **const           SMESH_RESTRICT);                          \
    template int sspyramid_extract_nodeset_from_sideset<IDX_T, ELEM_IDX_T>(const int,                                              \
                                                                           const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, \
                                                                           const ptrdiff_t,                                        \
                                                                           const ELEM_IDX_T *const SMESH_RESTRICT,                 \
                                                                           const i16 *const        SMESH_RESTRICT,                 \
                                                                           ptrdiff_t *,                                            \
                                                                           IDX_T **SMESH_RESTRICT)

SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SINGLE(i32);
SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SINGLE(i64);
SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SIDESET(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SIDESET(i64, i64);

#undef SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SIDESET
#undef SMESH_EXPLICIT_INSTANTIATE_SSPYRAMID_SINGLE

}  // namespace smesh
