#include "smesh_sshex8_graph.hpp"
#include "smesh_sshex8_graph.impl.hpp"

namespace smesh {

template int sshex8_generate_elements<i32>(
    const int, const ptrdiff_t, const ptrdiff_t,
    const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    i32 *const SMESH_RESTRICT *const SMESH_RESTRICT, ptrdiff_t *, ptrdiff_t *);

template int sshex8_crs_graph<i32, i32, i32>(
    const int, const ptrdiff_t, const ptrdiff_t,
    const i32 *const SMESH_RETRICT *const SMESH_RESTRICT, i32 **, i32 **);

template int sshex8_hierarchical_renumbering<i32>(
    const int, const int, int *const, const ptrdiff_t, const ptrdiff_t,
    i32 *const SMESH_RETRICT *const SMESH_RESTRICT);

template int sshex8_extract_surface_from_sideset<i32, i32>(
    const int, const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    const ptrdiff_t, const i32 *const SMESH_RESTRICT,
    const i16 *const SMESH_RESTRICT, i32 **const SMESH_RESTRICT);

template int sshex8_extract_nodeset_from_sideset<i32, i32>(
    const int, const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
    const ptrdiff_t, const i32 *const SMESH_RESTRICT,
    const i16 *const SMESH_RESTRICT, ptrdiff_t *, i32 **SMESH_RESTRICT);
    
} // namespace smesh