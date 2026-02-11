#ifndef SMESH_SSHEX8_GRAPH_HPP
#define SMESH_SSHEX8_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

#ifndef SMESH_RETRICT
#define SMESH_RETRICT SMESH_RESTRICT
#endif

template <typename element_idx_t, typename idx_t, typename count_t>
int hex8_build_edge_graph(
    const ptrdiff_t nelements, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    count_t **out_rowptr, idx_t **out_colidx);

template <typename element_idx_t, typename idx_t, typename count_t>
int sshex8_skeleton_crs_graph(
    const int L, const ptrdiff_t nelements, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RETRICT *const SMESH_RESTRICT elements,
    count_t **out_rowptr, idx_t **out_colidx);

template <typename idx_t>
ptrdiff_t nxe_max_node_id(const ptrdiff_t nelements, const int nxe,
                          idx_t **const SMESH_RESTRICT elements);

template <typename idx_t>
int sshex8_generate_elements(
    const int L, const ptrdiff_t m_nelements, const ptrdiff_t m_nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT m_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    ptrdiff_t *n_unique_nodes_out, ptrdiff_t *interior_start_out);

template <typename idx_t>
int sshex8_build_n2e(const int L, const ptrdiff_t nelements,
                     const ptrdiff_t nnodes,
                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                     count_t **out_n2eptr, element_idx_t **out_elindex);

template <typename count_t, typename idx_t, typename element_idx_t>
int sshex8_build_crs_graph_from_n2e(
    const int L, const ptrdiff_t nelements, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT elindex, count_t **out_rowptr,
    idx_t **out_colidx);

template <typename element_idx_t, typename count_t, typename idx_t>
int sshex8_crs_graph(
    const int L, const ptrdiff_t nelements, const ptrdiff_t nnodes,
    const idx_t *const SMESH_RETRICT *const SMESH_RESTRICT elements,
    count_t **out_rowptr, idx_t **out_colidx);

int sshex8_hierarchical_n_levels(const int L);

void sshex8_hierarchical_mesh_levels(const int L, const int nlevels,
                                     int *const levels);

template <typename idx_t>
int sshex8_hierarchical_renumbering(
    const int L, const int nlevels, int *const levels,
    const ptrdiff_t nelements, const ptrdiff_t nnodes,
    idx_t *const SMESH_RETRICT *const SMESH_RESTRICT elements);

template <typename idx_t, typename element_idx_t>
int sshex8_extract_surface_from_sideset(
    const int L, const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const i16 *const SMESH_RESTRICT side_idx,
    idx_t **const SMESH_RESTRICT sides);

template <typename idx_t, typename element_idx_t>
int sshex8_extract_nodeset_from_sideset(
    const int L, const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const i16 *const SMESH_RESTRICT side_idx, ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out);

template <typename idx_t, typename element_idx_t>
int sshex8_extract_quadshell4_surface_from_sideset(
    const int L, const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const i16 *const SMESH_RESTRICT side_idx,
    idx_t **const SMESH_RESTRICT sides);

template <typename idx_t>
int ssquad4_hierarchical_remapping(
    const int L, const int nlevels, int *const levels, const ptrdiff_t nelements,
    const ptrdiff_t nnodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t **SMESH_RESTRICT node_mapping_out, ptrdiff_t *count_out);

} // namespace smesh

#endif // SMESH_SSHEX8_GRAPH_HPP
