#ifndef SMESH_RESTRICTION_HPP
#define SMESH_RESTRICTION_HPP

#include <stddef.h>
#include <algorithm>
#include "smesh_base.hpp"
#include "smesh_alloc.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

namespace smesh {
    template <typename idx_t>
    ptrdiff_t max_node_id(const enum ElemType type, const ptrdiff_t nelements, idx_t **const SMESH_RESTRICT elements) {
        const int nxe = elem_num_nodes(type);
        idx_t     ret = 0;
        for (int i = 0; i < nxe; i++) {
            for (ptrdiff_t e = 0; e < nelements; e++) {
                ret = std::max(ret, elements[i][e]);
            }
        }
        return ret;
    }

    // We assume that fine indices have a higher id
    template <typename idx_t>
    int hierarchical_create_coarse_indices(const idx_t                 max_coarse_idx,
                                           const ptrdiff_t             n_indices,
                                           idx_t *const SMESH_RESTRICT fine_indices,
                                           ptrdiff_t                  *n_coarse_indices,
                                           idx_t **SMESH_RESTRICT      coarse_indices) {
        ptrdiff_t count = 0;
#pragma omp parallel for reduction(+ : count)
        for (ptrdiff_t i = 0; i < n_indices; i++) {
            count += fine_indices[i] <= max_coarse_idx;
        }

        *n_coarse_indices = count;
        *coarse_indices   = (idx_t *)SMESH_ALLOC(count * sizeof(idx_t));

        count = 0;
        for (ptrdiff_t i = 0; i < n_indices; i++) {
            if (fine_indices[i] <= max_coarse_idx) {
                (*coarse_indices)[count++] = fine_indices[i];
            }
        }

        return 0;
    }

    template <typename idx_t, typename real_t>
    int hierarchical_collect_coarse_values(const idx_t                        max_coarse_idx,
                                           const ptrdiff_t                    n_indices,
                                           idx_t *const SMESH_RESTRICT        fine_indices,
                                           const real_t *const SMESH_RESTRICT fine_values,
                                           real_t *const SMESH_RESTRICT       coarse_values) {
        ptrdiff_t count = 0;
        for (ptrdiff_t i = 0; i < n_indices; i++) {
            if (fine_indices[i] <= max_coarse_idx) {
                coarse_values[count++] = fine_values[i];
            }
        }

        return 0;
    }

    template <typename count_t, typename idx_t>
    int build_p1_to_p2_edge_map(const ptrdiff_t                     nnodes,
                                const count_t *const SMESH_RESTRICT coarse_rowptr,
                                const idx_t *const SMESH_RESTRICT   coarse_colidx,
                                idx_t *const SMESH_RESTRICT         p2_vertices) {
        ptrdiff_t fine_idx = nnodes;
        for (ptrdiff_t i = 0; i < nnodes; i++) {
            const count_t      start  = coarse_rowptr[i];
            const int          extent = coarse_rowptr[i + 1] - start;
            const idx_t *const cols   = &coarse_colidx[start];
            idx_t *const       verts  = &p2_vertices[start];

            for (int k = 0; k < extent; k++) {
                const idx_t j = cols[k];
                if (i < j) {
                    verts[k] = fine_idx++;
                } else {
                    verts[k] = invalid_idx<idx_t>();
                }
            }
        }

        return SMESH_SUCCESS;
    }

    template <typename count_t, typename idx_t, typename real_t>
    int hierarchical_restriction_with_edge_map(
            // CRS-node-graph of the coarse mesh
            const ptrdiff_t                     nnodes,
            const count_t *const SMESH_RESTRICT coarse_rowptr,
            const idx_t *const SMESH_RESTRICT   coarse_colidx,
            const idx_t *const SMESH_RESTRICT   p2_vertices,
            const int                           vec_size,
            const real_t *const SMESH_RESTRICT  from,
            real_t *const SMESH_RESTRICT        to) {
        // double tick = MPI_Wtime();

#pragma omp parallel for
        for (ptrdiff_t i = 0; i < nnodes; i++) {
            for (int v = 0; v < vec_size; v++) {
#pragma omp atomic update
                to[i * vec_size + v] += from[i * vec_size + v];
            }

            const count_t      start  = coarse_rowptr[i];
            const count_t      end    = coarse_rowptr[i + 1];
            const int          extent = end - start;
            const idx_t *const cols   = &coarse_colidx[start];
            const idx_t *const verts  = &p2_vertices[start];

            for (int k = 0; k < extent; k++) {
                const idx_t j    = cols[k];
                const idx_t edge = verts[k];
                if (i < j) {
                    for (int v = 0; v < vec_size; v++) {
                        const real_t edge_value = 0.5 * from[edge * vec_size + v];

#pragma omp atomic update
                        to[i * vec_size + v] += edge_value;

#pragma omp atomic update
                        to[j * vec_size + v] += edge_value;
                    }
                }
            }
        }

        // double tock = MPI_Wtime();
        // printf("hierarchical_restriction_with_edge_map %g [s]\n", tock - tick);
        return SMESH_SUCCESS;
    }

    template <typename idx_t, typename T>
    int hierarchical_restriction(const enum ElemType                                     from_element,
                                 const enum ElemType                                     to_element,
                                 const ptrdiff_t                                         nelements,
                                 const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                 const u16 *const SMESH_RESTRICT                         e2n_count,
                                 const int                                               vec_size,
                                 const T *const SMESH_RESTRICT                           from,
                                 T *const SMESH_RESTRICT                                 to);

    template <typename count_t, typename idx_t, typename T>
    int hierarchical_restriction(
            // CRS-node-graph of the coarse mesh
            const ptrdiff_t                     nnodes,
            const count_t *const SMESH_RESTRICT coarse_rowptr,
            const idx_t *const SMESH_RESTRICT   coarse_colidx,
            const int                           vec_size,
            const T *const SMESH_RESTRICT       from,
            T *const SMESH_RESTRICT             to);

}  // namespace smesh

#endif  // SMESH_RESTRICTION_HPP
