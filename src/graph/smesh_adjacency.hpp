#ifndef SMESH_ADJACENCY_HPP
#define SMESH_ADJACENCY_HPP

#include <cstddef>
#include <cstdint>

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

namespace smesh {

    struct LocalSideTable {
        static constexpr int MAX_NUM_SIDES          = 6;
        static constexpr int MAX_NUM_NODES_PER_SIDE = 9;
        int                  nnxs{-1};
        int                  table[MAX_NUM_SIDES * MAX_NUM_NODES_PER_SIDE];

        inline void fill(enum ElemType element_type) {
            enum ElemType st = side_type(element_type);
            this->nnxs       = elem_num_nodes(st);
            if (element_type == TET10 || element_type == TET4) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;
                (*this)(0, 2) = 4 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;
                (*this)(1, 2) = 4 - 1;

                (*this)(2, 0) = 1 - 1;
                (*this)(2, 1) = 4 - 1;
                (*this)(2, 2) = 3 - 1;

                (*this)(3, 0) = 1 - 1;
                (*this)(3, 1) = 3 - 1;
                (*this)(3, 2) = 2 - 1;

                if (element_type == TET10) {
                    (*this)(0, 3) = 5 - 1;
                    (*this)(0, 4) = 9 - 1;
                    (*this)(0, 5) = 8 - 1;

                    (*this)(1, 3) = 6 - 1;
                    (*this)(1, 4) = 10 - 1;
                    (*this)(1, 5) = 9 - 1;

                    (*this)(2, 3) = 8 - 1;
                    (*this)(2, 4) = 10 - 1;
                    (*this)(2, 5) = 7 - 1;

                    (*this)(3, 3) = 7 - 1;
                    (*this)(3, 4) = 6 - 1;
                    (*this)(3, 5) = 5 - 1;
                }

            } else if (element_type == TRI3) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 1 - 1;
            } else if (element_type == TRI6) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;
                (*this)(0, 2) = 4 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;
                (*this)(1, 2) = 5 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 1 - 1;
                (*this)(2, 2) = 6 - 1;
            } else if (element_type == QUAD4) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 4 - 1;

                (*this)(3, 0) = 4 - 1;
                (*this)(3, 1) = 1 - 1;
            } else if (element_type == QUAD9 || element_type == QUADSHELL9) {
                (*this)(0, 0) = 0;
                (*this)(0, 1) = 1;
                (*this)(0, 2) = 4;

                (*this)(1, 0) = 1;
                (*this)(1, 1) = 2;
                (*this)(1, 2) = 5;

                (*this)(2, 0) = 2;
                (*this)(2, 1) = 3;
                (*this)(2, 2) = 6;

                (*this)(3, 0) = 3;
                (*this)(3, 1) = 0;
                (*this)(3, 2) = 7;
            } else if (element_type == HEX8) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;
                (*this)(0, 2) = 6 - 1;
                (*this)(0, 3) = 5 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;
                (*this)(1, 2) = 7 - 1;
                (*this)(1, 3) = 6 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 4 - 1;
                (*this)(2, 2) = 8 - 1;
                (*this)(2, 3) = 7 - 1;

                (*this)(3, 0) = 4 - 1;
                (*this)(3, 1) = 1 - 1;
                (*this)(3, 2) = 5 - 1;
                (*this)(3, 3) = 8 - 1;

                (*this)(4, 0) = 4 - 1;
                (*this)(4, 1) = 3 - 1;
                (*this)(4, 2) = 2 - 1;
                (*this)(4, 3) = 1 - 1;

                (*this)(5, 0) = 5 - 1;
                (*this)(5, 1) = 6 - 1;
                (*this)(5, 2) = 7 - 1;
                (*this)(5, 3) = 8 - 1;
            } else if (element_type == HEX27) {
                // Corners 0..7, edge nodes 8..19, face centers 20..25,
                // and volume center 26. Face order follows HEX8 above.
                const int faces[6][9] = {
                        {0, 1, 5, 4, 8, 17, 12, 16, 20},
                        {1, 2, 6, 5, 9, 18, 13, 17, 21},
                        {2, 3, 7, 6, 10, 19, 14, 18, 22},
                        {3, 0, 4, 7, 11, 16, 15, 19, 23},
                        {3, 2, 1, 0, 10, 9, 8, 11, 24},
                        {4, 5, 6, 7, 12, 13, 14, 15, 25},
                };
                for (int side = 0; side < 6; ++side) {
                    for (int node = 0; node < 9; ++node) {
                        (*this)(side, node) = faces[side][node];
                    }
                }
            } else if (element_type == PROTEUS_HEX8) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;
                (*this)(0, 2) = 6 - 1;
                (*this)(0, 3) = 5 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 4 - 1;
                (*this)(1, 2) = 8 - 1;
                (*this)(1, 3) = 6 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 7 - 1;
                (*this)(2, 2) = 8 - 1;
                (*this)(2, 3) = 4 - 1;

                (*this)(3, 0) = 1 - 1;
                (*this)(3, 1) = 5 - 1;
                (*this)(3, 2) = 7 - 1;
                (*this)(3, 3) = 3 - 1;

                (*this)(4, 0) = 1 - 1;
                (*this)(4, 1) = 3 - 1;
                (*this)(4, 2) = 4 - 1;
                (*this)(4, 3) = 2 - 1;

                (*this)(5, 0) = 5 - 1;
                (*this)(5, 1) = 6 - 1;
                (*this)(5, 2) = 8 - 1;
                (*this)(5, 3) = 7 - 1;
            } else {
                SMESH_ERROR("fill_local_side_table: Unsupported element type: %s\n", type_to_string(element_type));
            }
        }

        inline int &operator()(const int side, const int node) { return table[side * nnxs + node]; }

        inline int operator()(const int side, const int node) const { return table[side * nnxs + node]; }
    };

    template <typename idx_t, typename count_t, typename element_idx_t>
    void create_element_adj_table_from_dual_graph(const ptrdiff_t                                         n_elements,
                                                  enum ElemType                                           element_type,
                                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                                  const count_t *const                                    adj_ptr,
                                                  const element_idx_t *const                              adj_idx,
                                                  element_idx_t *const SMESH_RESTRICT                     table);

    template <typename idx_t, typename count_t, typename element_idx_t>
    void create_element_adj_table_from_dual_graph_soa(const ptrdiff_t                                           n_elements,
                                                      enum ElemType                                             element_type,
                                                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT   elems,
                                                      const count_t *const                                      adj_ptr,
                                                      const element_idx_t *const                                adj_idx,
                                                      element_idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT table);

    template <typename idx_t, typename count_t = idx_t, typename element_idx_t>
    void create_element_adj_table(const ptrdiff_t                                         n_elements,
                                  const ptrdiff_t                                         n_nodes,
                                  enum ElemType                                           element_type,
                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                  element_idx_t **SMESH_RESTRICT                          table_out);

    template <typename idx_t, typename count_t = idx_t, typename element_idx_t>
    void extract_surface_connectivity_with_adj_table(const ptrdiff_t                                         n_elements,
                                                     const ptrdiff_t                                         n_nodes,
                                                     enum ElemType                                           element_type,
                                                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                                                     ptrdiff_t                                              *n_surf_elements,
                                                     idx_t **SMESH_RESTRICT                                  surf_elems,
                                                     element_idx_t **SMESH_RESTRICT                          parent_element);

    template <typename element_idx_t>
    int extract_sideset_from_adj_table(const enum ElemType                       element_type,
                                       const ptrdiff_t                           n_elements,
                                       const element_idx_t *const SMESH_RESTRICT table,
                                       ptrdiff_t *SMESH_RESTRICT                 n_surf_elements,
                                       element_idx_t **SMESH_RESTRICT            parent_element,
                                       i16 **SMESH_RESTRICT                      side_idx);

    template <typename idx_t, typename count_t = idx_t, typename element_idx_t>
    int extract_skin_sideset(const ptrdiff_t                                         n_elements,
                             const ptrdiff_t                                         n_nodes,
                             const enum ElemType                                     element_type,
                             const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
                             ptrdiff_t *SMESH_RESTRICT                               n_surf_elements,
                             element_idx_t **SMESH_RESTRICT                          parent_element,
                             i16 **SMESH_RESTRICT                                    side_idx);


}  // namespace smesh

#endif  // SMESH_ADJACENCY_HPP
