#ifndef SMESH_PROLONGATION_HPP
#define SMESH_PROLONGATION_HPP

#include <stddef.h>
#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {
    template <typename idx_t, typename T>
    int hierarchical_prolongation(const enum ElemType                                     from_element,
                                  const enum ElemType                                     to_element,
                                  const ptrdiff_t                                         nelements,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                  const int                                               vec_size,
                                  const T *const SMESH_RESTRICT                           from,
                                  T *const SMESH_RESTRICT                                 to);

    template <typename count_t, typename idx_t, typename T>
    int hierarchical_prolongation_with_edge_map(const ptrdiff_t                     nnodes,
                                                const count_t *const SMESH_RESTRICT coarse_rowptr,
                                                const idx_t *const SMESH_RESTRICT   coarse_colidx,
                                                const idx_t *const SMESH_RESTRICT   p2_vertices,
                                                const int                           vec_size,
                                                const T *const SMESH_RESTRICT  from,
                                                T *const SMESH_RESTRICT        to) {
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < nnodes; i++) {
            const ptrdiff_t i_offset = i * vec_size;
            for (int v = 0; v < vec_size; v++) {
                to[i_offset + v] = from[i_offset + v];
            }

            const count_t      start  = coarse_rowptr[i];
            const count_t      end    = coarse_rowptr[i + 1];
            const int          extent = end - start;
            const idx_t *const cols   = &coarse_colidx[start];
            const idx_t *const verts  = &p2_vertices[start];

            for (int k = 0; k < extent; k++) {
                const ptrdiff_t j    = cols[k];
                const idx_t     edge = verts[k];

                if (i < j) {
                    assert(edge >= nnodes);

                    const ptrdiff_t edge_offset = edge * vec_size;
                    const ptrdiff_t j_offset    = j * vec_size;

                    for (int v = 0; v < vec_size; v++) {
                        const T edge_value = 0.5 * (from[i_offset + v] + from[j_offset + v]);

                        to[edge_offset + v] = edge_value;
                    }
                }
            }
        }

        return SMESH_SUCCESS;
    }

}  // namespace smesh

#endif  // SMESH_PROLONGATION_HPP
