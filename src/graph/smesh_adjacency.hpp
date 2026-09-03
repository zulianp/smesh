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
        int                  nnxs_side[MAX_NUM_SIDES];
        int                  table[MAX_NUM_SIDES * MAX_NUM_NODES_PER_SIDE];

        inline int &operator()(const int side, const int node) {
            return table[side * MAX_NUM_NODES_PER_SIDE + node];
        }

        inline int operator()(const int side, const int node) const {
            return table[side * MAX_NUM_NODES_PER_SIDE + node];
        }

        inline void fill_side_arities(enum ElemType element_type) {
            const int ns = elem_num_sides(element_type);
            if (elem_sides_homogeneous(element_type)) {
                this->nnxs = elem_num_nodes(side_type(element_type));
                for (int s = 0; s < ns; ++s) {
                    nnxs_side[s] = this->nnxs;
                }
            } else {
                this->nnxs = -1;
                for (int s = 0; s < ns; ++s) {
                    nnxs_side[s] = elem_num_nodes(side_type(element_type, s));
                }
            }
        }

        static void report_unsupported(const char *operation, enum ElemType element_type) {
            fprintf(stderr,
                    "%s: LocalSideTable has no layout for %s "
                    "(supported: HEX/TET/WEDGE/PYRAMID/QUAD/TRI/EDGE and shell aliases)\n",
                    operation,
                    type_to_string(element_type));
        }

        static bool supported(enum ElemType element_type) {
            LocalSideTable table;
            return table.fill(element_type) == SMESH_SUCCESS;
        }

        inline int fill(enum ElemType element_type) {
            enum ElemType arities = element_type;
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
                    // Exodus TET10 mids: 4 on 0-1, 5 on 1-2, 6 on 2-0, 7 on 0-3, 8 on 1-3, 9 on 2-3.
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

            } else if (element_type == TRI3 || element_type == TRISHELL3) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 1 - 1;
                arities = TRI3;
            } else if (element_type == TRI6 || element_type == TRISHELL6) {
                // Exodus TRI6: mids 3 on 0-1, 4 on 1-2, 5 on 2-0.
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;
                (*this)(0, 2) = 4 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;
                (*this)(1, 2) = 5 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 1 - 1;
                (*this)(2, 2) = 6 - 1;
                arities = TRI6;
            } else if (element_type == QUAD4 || element_type == QUADSHELL4) {
                (*this)(0, 0) = 1 - 1;
                (*this)(0, 1) = 2 - 1;

                (*this)(1, 0) = 2 - 1;
                (*this)(1, 1) = 3 - 1;

                (*this)(2, 0) = 3 - 1;
                (*this)(2, 1) = 4 - 1;

                (*this)(3, 0) = 4 - 1;
                (*this)(3, 1) = 1 - 1;
                arities = QUAD4;
            } else if (element_type == QUAD9 || element_type == QUADSHELL9) {
                // Exodus QUAD9 edges: mids 4,5,6,7 on 0-1, 1-2, 2-3, 3-0. Face center 8 is not a side node.
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
                // Exodus/PATRAN face centers: 20=-y, 21=+x, 22=+y, 23=-x, 24=-z, 25=+z.
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
            } else if (element_type == WEDGE6) {
                // Exodus/VTK: bottom 0,1,2, top 3,4,5. Sides 0-2 QUAD4, 3-4 TRI3.
                (*this)(0, 0) = 0;
                (*this)(0, 1) = 1;
                (*this)(0, 2) = 4;
                (*this)(0, 3) = 3;

                (*this)(1, 0) = 1;
                (*this)(1, 1) = 2;
                (*this)(1, 2) = 5;
                (*this)(1, 3) = 4;

                (*this)(2, 0) = 2;
                (*this)(2, 1) = 0;
                (*this)(2, 2) = 3;
                (*this)(2, 3) = 5;

                (*this)(3, 0) = 0;
                (*this)(3, 1) = 1;
                (*this)(3, 2) = 2;

                (*this)(4, 0) = 3;
                (*this)(4, 1) = 4;
                (*this)(4, 2) = 5;
            } else if (element_type == PYRAMID5) {
                // Base 0,1,2,3 (HEX-compatible), apex 4. Sides 0-3 TRI3, side 4 QUAD4.
                (*this)(0, 0) = 0;
                (*this)(0, 1) = 1;
                (*this)(0, 2) = 4;

                (*this)(1, 0) = 1;
                (*this)(1, 1) = 2;
                (*this)(1, 2) = 4;

                (*this)(2, 0) = 2;
                (*this)(2, 1) = 3;
                (*this)(2, 2) = 4;

                (*this)(3, 0) = 3;
                (*this)(3, 1) = 0;
                (*this)(3, 2) = 4;

                (*this)(4, 0) = 0;
                (*this)(4, 1) = 1;
                (*this)(4, 2) = 2;
                (*this)(4, 3) = 3;
            } else if (element_type == EDGE2 || element_type == EDGESHELL2 || element_type == EDGE3 ||
                       element_type == EDGESHELL3) {
                (*this)(0, 0) = 0;
                (*this)(1, 0) = 1;
            } else {
                return SMESH_FAILURE;
            }

            fill_side_arities(arities);
            return SMESH_SUCCESS;
        }
    };

    struct LocalEdgeTable {
        static constexpr int MAX_NUM_EDGES          = 12;
        static constexpr int MAX_NUM_NODES_PER_EDGE = 3;
        int                  nnxe{-1};
        int                  nnxe_edge[MAX_NUM_EDGES];
        int                  table[MAX_NUM_EDGES * MAX_NUM_NODES_PER_EDGE];

        inline int &operator()(const int edge, const int node) {
            return table[edge * MAX_NUM_NODES_PER_EDGE + node];
        }

        inline int operator()(const int edge, const int node) const {
            return table[edge * MAX_NUM_NODES_PER_EDGE + node];
        }

        inline void fill_edge_arities(enum ElemType element_type, const int nne) {
            const int ne = elem_num_edges(element_type);
            this->nnxe   = nne;
            for (int e = 0; e < ne; ++e) {
                nnxe_edge[e] = nne;
            }
        }

        static void report_unsupported(const char *operation, enum ElemType element_type) {
            fprintf(stderr,
                    "%s: LocalEdgeTable has no layout for %s "
                    "(supported: HEX/TET/WEDGE/PYRAMID/QUAD/TRI/EDGE and shell aliases)\n",
                    operation,
                    type_to_string(element_type));
        }

        static bool supported(enum ElemType element_type) {
            LocalEdgeTable table;
            return table.fill(element_type) == SMESH_SUCCESS;
        }

        inline int fill(enum ElemType element_type) {
            if (is_semistructured_type(element_type)) {
                const enum ElemType fam = ss_source_family(element_type);
                if (fam != element_type) {
                    return fill(fam);
                }
            }

            if (element_type == TRI3 || element_type == TRISHELL3 || element_type == TRI6 ||
                element_type == TRISHELL6 || element_type == QUAD4 || element_type == QUADSHELL4 ||
                element_type == QUAD9 || element_type == QUADSHELL9) {
                LocalSideTable lst;
                if (lst.fill(element_type) != SMESH_SUCCESS) {
                    return SMESH_FAILURE;
                }
                const int ne = elem_num_edges(element_type);
                this->nnxe   = lst.nnxs;
                for (int e = 0; e < ne; ++e) {
                    nnxe_edge[e] = lst.nnxs_side[e];
                    for (int n = 0; n < nnxe_edge[e]; ++n) {
                        (*this)(e, n) = lst(e, n);
                    }
                }
                return SMESH_SUCCESS;
            }

            if (element_type == EDGE2 || element_type == EDGESHELL2 || element_type == BEAM2) {
                (*this)(0, 0) = 0;
                (*this)(0, 1) = 1;
                fill_edge_arities(element_type, 2);
                return SMESH_SUCCESS;
            }
            if (element_type == EDGE3 || element_type == EDGESHELL3) {
                (*this)(0, 0) = 0;
                (*this)(0, 1) = 1;
                (*this)(0, 2) = 2;
                fill_edge_arities(element_type, 3);
                return SMESH_SUCCESS;
            }

            if (element_type == HEX8 || element_type == PROTEUS_HEX8) {
                const int edges[12][2] = {{0, 1},
                                          {1, 2},
                                          {2, 3},
                                          {3, 0},
                                          {4, 5},
                                          {5, 6},
                                          {6, 7},
                                          {7, 4},
                                          {0, 4},
                                          {1, 5},
                                          {2, 6},
                                          {3, 7}};
                for (int e = 0; e < 12; ++e) {
                    (*this)(e, 0) = edges[e][0];
                    (*this)(e, 1) = edges[e][1];
                }
                fill_edge_arities(HEX8, 2);
                return SMESH_SUCCESS;
            }

            if (element_type == HEX27) {
                const int edges[12][3] = {{0, 1, 8},
                                          {1, 2, 9},
                                          {2, 3, 10},
                                          {3, 0, 11},
                                          {4, 5, 12},
                                          {5, 6, 13},
                                          {6, 7, 14},
                                          {7, 4, 15},
                                          {0, 4, 16},
                                          {1, 5, 17},
                                          {2, 6, 18},
                                          {3, 7, 19}};
                for (int e = 0; e < 12; ++e) {
                    (*this)(e, 0) = edges[e][0];
                    (*this)(e, 1) = edges[e][1];
                    (*this)(e, 2) = edges[e][2];
                }
                fill_edge_arities(HEX27, 3);
                return SMESH_SUCCESS;
            }

            if (element_type == TET4 || element_type == TET10) {
                const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
                for (int e = 0; e < 6; ++e) {
                    (*this)(e, 0) = edges[e][0];
                    (*this)(e, 1) = edges[e][1];
                }
                if (element_type == TET10) {
                    const int mids[6] = {4, 5, 6, 7, 8, 9};
                    for (int e = 0; e < 6; ++e) {
                        (*this)(e, 2) = mids[e];
                    }
                    fill_edge_arities(TET10, 3);
                } else {
                    fill_edge_arities(TET4, 2);
                }
                return SMESH_SUCCESS;
            }

            if (element_type == WEDGE6) {
                const int edges[9][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
                for (int e = 0; e < 9; ++e) {
                    (*this)(e, 0) = edges[e][0];
                    (*this)(e, 1) = edges[e][1];
                }
                fill_edge_arities(WEDGE6, 2);
                return SMESH_SUCCESS;
            }

            if (element_type == PYRAMID5) {
                const int edges[8][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
                for (int e = 0; e < 8; ++e) {
                    (*this)(e, 0) = edges[e][0];
                    (*this)(e, 1) = edges[e][1];
                }
                fill_edge_arities(PYRAMID5, 2);
                return SMESH_SUCCESS;
            }

            return SMESH_FAILURE;
        }
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
