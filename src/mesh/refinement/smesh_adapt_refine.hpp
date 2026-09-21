#ifndef SMESH_ADAPT_REFINE_HPP
#define SMESH_ADAPT_REFINE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

inline bool adapt_refine_type_supported(const enum ElemType t) {
    return t == TET4 || t == TRI3 || t == TRISHELL3 || t == QUAD4 || t == QUADSHELL4;
}

/// Conforming local refine. TRI/TET: geometry edges with ℓ > min(h_a,h_b) are
/// marked directly. Grading/quality marks the longest edge of an element that
/// violates the size field or has mean-ratio < `q_min`. Longest marked first;
/// every incident element splits that same edge. QUAD: 4-split with 2:1 closure.
/// `element_mark` (nullable) is over the *input* elements. Output arrays are
/// SMESH_ALLOC'd; free with `mesh_adapt_refine_free`.
template <typename idx_t, typename count_t, typename geom_t>
int mesh_adapt_refine(const enum ElemType                                     element_type,
                      const ptrdiff_t                                         n_elements,
                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                      const int                                               sdim,
                      const ptrdiff_t                                         n_nodes,
                      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                      const geom_t *const SMESH_RESTRICT                      h,
                      const uint8_t *const SMESH_RESTRICT                     element_mark,
                      const geom_t                                            q_min,
                      const uint8_t *const SMESH_RESTRICT                     geom_node,
                      const int                                               max_levels,
                      ptrdiff_t                                              *n_elements_out,
                      idx_t                                                ***elements_out,
                      ptrdiff_t                                              *n_nodes_out,
                      geom_t                                               ***points_out,
                      ptrdiff_t                                             **parent_elem_out,
                      count_t                                               **parent_ptr_out,
                      idx_t                                                 **child_id_out,
                      idx_t                                                 **node_a_out,
                      idx_t                                                 **node_b_out);

void mesh_adapt_refine_free(const int         nxe,
                            const int         sdim,
                            idx_t           **elements,
                            geom_t          **points,
                            ptrdiff_t        *parent_elem,
                            void             *parent_ptr,
                            idx_t            *child_id,
                            idx_t            *node_a,
                            idx_t            *node_b);

}  // namespace smesh

#endif
