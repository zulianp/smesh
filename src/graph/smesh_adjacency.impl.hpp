#ifndef SMESH_ADJACENCY_IMPL_HPP
#define SMESH_ADJACENCY_IMPL_HPP

#include "smesh_graph.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_sort.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <iterator>

namespace smesh {


template <typename idx_t, typename count_t, typename element_idx_t>
void create_element_adj_table_from_dual_graph(
    const ptrdiff_t n_elements, enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const adj_ptr,
    const element_idx_t *const adj_idx,
    element_idx_t *const SMESH_RESTRICT table) {
  enum ElemType element_type_for_algo = element_type;

  if (element_type == TET10) {
    // This is enough for many operations
    element_type_for_algo = TET4;
  } else if (element_type == TRI6) {
    element_type_for_algo = TRI3;
  }

  LocalSideTable lst;
  lst.fill(element_type_for_algo);

  enum ElemType st = side_type(element_type_for_algo);
  const int nn = elem_num_nodes(st);
  const int ns = elem_num_sides(element_type_for_algo);

#pragma omp parallel for
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    idx_t nodes1[LocalSideTable::MAX_NUM_NODES_PER_SIDE];
    idx_t nodes2[LocalSideTable::MAX_NUM_NODES_PER_SIDE];
    int assigned[LocalSideTable::MAX_NUM_SIDES];

    const count_t begin = adj_ptr[e];
    const count_t end = adj_ptr[e + 1];
    const count_t range = end - begin;

    memset(assigned, 0, range * sizeof(int));

    for (int s1 = 0; s1 < ns; s1++) {
      table[e * ns + s1] = invalid_idx<element_idx_t>();

      for (int j = 0; j < nn; j++) {
        nodes1[j] = elems[lst(s1, j)][e];
      }

      std::sort(nodes1, nodes1 + nn);

      for (count_t k = 0; k < range; k++) {
        if (assigned[k])
          continue;
        const element_idx_t e_adj = adj_idx[begin + k];

        for (int s2 = 0; s2 < ns; s2++) {
          for (int j = 0; j < nn; j++) {
            nodes2[j] = elems[lst(s2, j)][e_adj];
          }

          std::sort(nodes2, nodes2 + nn);

          int diffs = 0;
          for (int j = 0; j < nn; j++) {
            diffs += nodes1[j] != nodes2[j];
          }

          if (!diffs) {
            // Array of structures
            table[e * ns + s1] = e_adj;
            assigned[k] = 1;
            break;
          }
        }
      }
    }
  }
}

template <typename idx_t, typename count_t, typename element_idx_t>
void create_element_adj_table_from_dual_graph_soa(
    const ptrdiff_t n_elements, enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const count_t *const adj_ptr, const element_idx_t *const adj_idx,
    element_idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT table) {
  enum ElemType element_type_for_algo = element_type;

  if (element_type == TET10) {
    // This is enough for many operations
    element_type_for_algo = TET4;
  } else if (element_type == TRI6) {
    element_type_for_algo = TRI3;
  }

  LocalSideTable lst;
  lst.fill(element_type_for_algo);

  idx_t nodes1[LocalSideTable::MAX_NUM_NODES_PER_SIDE];
  idx_t nodes2[LocalSideTable::MAX_NUM_NODES_PER_SIDE];
  int assigned[LocalSideTable::MAX_NUM_SIDES];

  enum ElemType st = side_type(element_type_for_algo);
  const int nn = elem_num_nodes(st);
  const int ns = elem_num_sides(element_type_for_algo);

  for (ptrdiff_t e = 0; e < n_elements; e++) {
    const count_t begin = adj_ptr[e];
    const count_t end = adj_ptr[e + 1];
    const count_t range = end - begin;

    memset(assigned, 0, range * sizeof(int));

    for (int s1 = 0; s1 < ns; s1++) {
      table[s1][e] = invalid_idx<element_idx_t>();

      for (int j = 0; j < nn; j++) {
        nodes1[j] = elems[lst(s1, j)][e];
      }

      std::sort(nodes1, nodes1 + nn);

      for (count_t k = 0; k < range; k++) {
        if (assigned[k])
          continue;
        const element_idx_t e_adj = adj_idx[begin + k];

        for (int s2 = 0; s2 < ns; s2++) {
          for (int j = 0; j < nn; j++) {
            nodes2[j] = elems[lst(s2, j)][e_adj];
          }

          std::sort(nodes2, nodes2 + nn);

          int diffs = 0;
          for (int j = 0; j < nn; j++) {
            diffs += nodes1[j] != nodes2[j];
          }

          if (!diffs) {
            // Array of structures
            table[s1][e] = e_adj;
            assigned[k] = 1;
            break;
          }
        }
      }
    }
  }
}

template <typename idx_t, typename count_t, typename element_idx_t>
void create_element_adj_table(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    element_idx_t **SMESH_RESTRICT table_out) {
  enum ElemType element_type_for_algo = element_type;

  if (element_type == TET10) {
    // This is enough for many operations
    element_type_for_algo = TET4;
  } else if (element_type == TRI6) {
    element_type_for_algo = TRI3;
  }

  count_t *adj_ptr = 0;
  element_idx_t *adj_idx = 0;
  create_dual_graph(n_elements, n_nodes, element_type_for_algo, elems, &adj_ptr,
                    &adj_idx);

  const int ns = elem_num_sides(element_type);
  element_idx_t *table =
      (element_idx_t *)malloc(n_elements * ns * sizeof(element_idx_t));
  create_element_adj_table_from_dual_graph(n_elements, element_type, elems,
                                           adj_ptr, adj_idx, table);

  free(adj_ptr);
  free(adj_idx);

  *table_out = table;
}

template <typename idx_t, typename count_t, typename element_idx_t>
void extract_surface_connectivity_with_adj_table(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    ptrdiff_t *n_surf_elements,
    idx_t **SMESH_RESTRICT surf_elems,
    element_idx_t **SMESH_RESTRICT parent_element)

{
  double tick = time_seconds();

  const int ns = elem_num_sides(element_type);
  element_idx_t *table = 0;
  create_element_adj_table<idx_t, count_t, element_idx_t>(
      n_elements, n_nodes, element_type, elems, &table);

  LocalSideTable lst;
  lst.fill(element_type);
  enum ElemType st = side_type(element_type);
  const int nn = elem_num_nodes(st);

  *n_surf_elements = 0;
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    for (int s = 0; s < ns; s++) {
      // Array of structures
      const element_idx_t e_adj = table[e * ns + s];
      if (e_adj == invalid_idx<element_idx_t>()) {
        (*n_surf_elements)++;
      }
    }
  }

  *parent_element =
      (element_idx_t *)malloc((*n_surf_elements) * sizeof(element_idx_t));
  for (int s = 0; s < nn; s++) {
    surf_elems[s] = (idx_t *)malloc((*n_surf_elements) * sizeof(idx_t));
  }

  ptrdiff_t side_offset = 0;
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    for (int s = 0; s < ns; s++) {
      // Array of structures
      const element_idx_t e_adj = table[e * ns + s];
      if (e_adj == invalid_idx<element_idx_t>()) {
        for (int n = 0; n < nn; n++) {
          idx_t node = elems[lst(s, n)][e];
          surf_elems[n][side_offset] = node;
        }

        (*parent_element)[side_offset] = e;
        side_offset++;
      }
    }
  }

  free(table);

  double tock = time_seconds();
  printf(
      "adj_table.c: extract_surface_connectivity_with_adj_table\t%g seconds\n",
      tock - tick);
}

template <typename element_idx_t>
int extract_sideset_from_adj_table(
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const element_idx_t *const SMESH_RESTRICT table,
    ptrdiff_t *SMESH_RESTRICT n_surf_elements,
    element_idx_t **SMESH_RESTRICT parent_element,
    int16_t **SMESH_RESTRICT side_idx) {
  double tick = time_seconds();

  LocalSideTable lst;
  lst.fill(element_type);
  const int ns = elem_num_sides(element_type);

  *n_surf_elements = 0;
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    for (int s = 0; s < ns; s++) {
      const element_idx_t e_adj = table[e * ns + s];
      if (e_adj == invalid_idx<element_idx_t>()) {
        (*n_surf_elements)++;
      }
    }
  }

  *parent_element =
      (element_idx_t *)malloc((*n_surf_elements) * sizeof(element_idx_t));
  *side_idx = (i16 *)malloc((*n_surf_elements) * sizeof(i16));

  ptrdiff_t side_offset = 0;
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    for (int s = 0; s < ns; s++) {
      const element_idx_t e_adj = table[e * ns + s];
      if (e_adj == invalid_idx<element_idx_t>()) {
        (*parent_element)[side_offset] = e;
        (*side_idx)[side_offset] = s;
        side_offset++;
      }
    }
  }

  double tock = time_seconds();
  printf("adj_table.c: extract_sideset_from_adj_table\t%g seconds\n",
         tock - tick);

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int extract_skin_sideset(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    ptrdiff_t *SMESH_RESTRICT n_surf_elements,
    element_idx_t **SMESH_RESTRICT parent_element,
    int16_t **SMESH_RESTRICT side_idx)

{
  element_idx_t *table = 0;
  create_element_adj_table<idx_t, count_t, element_idx_t>(
      n_elements, n_nodes, element_type, elems, &table);
  int err =
      extract_sideset_from_adj_table(element_type, n_elements, table,
                                     n_surf_elements, parent_element, side_idx);
  free(table);
  return err;
}

} // namespace smesh

#endif // SMESH_ADJACENCY_IMPL_HPP
