#ifndef SMESH_VOLUME_TO_SURFACE_IMPL_HPP
#define SMESH_VOLUME_TO_SURFACE_IMPL_HPP

#include "smesh_adjacency.hpp"
#include "smesh_volume_to_surface.hpp"

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int extract_skin_sideset_from_n2e(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const SMESH_RESTRICT n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx,
    ptrdiff_t *SMESH_RESTRICT n_surf_elements,
    element_idx_t **SMESH_RESTRICT parent_element,
    i16 **SMESH_RESTRICT side_idx) {
    (void)n_nodes;

    LocalSideTable lst;
    lst.fill(element_type);

    const int ns = elem_num_sides(element_type);
    const int nnxe = elem_num_nodes(element_type);
    const int nnxs = elem_num_nodes(side_type(element_type));
    unsigned char *const is_boundary =
        (unsigned char *)malloc((size_t)n_elements * (size_t)ns * sizeof(unsigned char));

    ptrdiff_t surf_count = 0;

#pragma omp parallel for reduction(+ : surf_count)
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
        idx_t side_nodes[LocalSideTable::MAX_NUM_NODES_PER_SIDE];

        for (int s = 0; s < ns; ++s) {
            idx_t pivot_node = elems[lst(s, 0)][e];
            count_t pivot_begin = n2e_ptr[pivot_node];
            count_t pivot_end = n2e_ptr[pivot_node + 1];
            count_t pivot_degree = pivot_end - pivot_begin;

            side_nodes[0] = pivot_node;

            for (int n = 1; n < nnxs; ++n) {
                const idx_t node = elems[lst(s, n)][e];
                side_nodes[n] = node;

                const count_t begin = n2e_ptr[node];
                const count_t end = n2e_ptr[node + 1];
                const count_t degree = end - begin;
                if (degree < pivot_degree) {
                    pivot_node = node;
                    pivot_begin = begin;
                    pivot_end = end;
                    pivot_degree = degree;
                }
            }

            unsigned char boundary = 1;
            for (count_t it = pivot_begin; it < pivot_end; ++it) {
                const element_idx_t e_adj = n2e_idx[it];
                if (e_adj == e) {
                    continue;
                }

                int matches = 0;
                for (int n = 0; n < nnxs; ++n) {
                    const idx_t side_node = side_nodes[n];
                    int found = 0;

                    for (int en = 0; en < nnxe; ++en) {
                        found += (elems[en][e_adj] == side_node);
                    }

                    if (!found) {
                        break;
                    }

                    matches++;
                }

                if (matches == nnxs) {
                    boundary = 0;
                    break;
                }
            }

            is_boundary[e * ns + s] = boundary;
            surf_count += boundary;
        }
    }

    *n_surf_elements = surf_count;
    *parent_element = (element_idx_t *)malloc((size_t)surf_count * sizeof(element_idx_t));
    *side_idx = (i16 *)malloc((size_t)surf_count * sizeof(i16));

    ptrdiff_t side_offset = 0;
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
        for (int s = 0; s < ns; ++s) {
            if (!is_boundary[e * ns + s]) {
                continue;
            }

            (*parent_element)[side_offset] = e;
            (*side_idx)[side_offset] = (i16)s;
            side_offset++;
        }
    }

    free(is_boundary);
    return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_VOLUME_TO_SURFACE_IMPL_HPP
