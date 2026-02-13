#ifndef SMESH_SIDESETS_IMPL_HPP
#define SMESH_SIDESETS_IMPL_HPP

#include "smesh_adjacency.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_mask.hpp"
#include "smesh_sidesets.hpp"
#include "smesh_sort.hpp"

namespace smesh {

template <typename element_idx_t>
int sideset_to_e2s_count(const ptrdiff_t n_sides, const ptrdiff_t n_elements,
                         const element_idx_t *const SMESH_RESTRICT
                             parent_element,
                         element_idx_t *const SMESH_RESTRICT e2s_ptr) {
  memset(e2s_ptr, 0, (n_elements + 1) * sizeof(element_idx_t));
  for (ptrdiff_t i = 0; i < n_sides; i++) {
    const element_idx_t e = parent_element[i];
    e2s_ptr[e + 1]++;
  }

  for (ptrdiff_t i = 0; i < n_elements; i++) {
    e2s_ptr[i + 1] += e2s_ptr[i];
  }

  return SMESH_SUCCESS;
}

template <typename element_idx_t>
int sideset_to_e2s_fill(const ptrdiff_t n_sides, const ptrdiff_t n_elements,
                        const element_idx_t *const SMESH_RESTRICT
                            parent_element,
                        const element_idx_t *const SMESH_RESTRICT e2s_ptr,
                        element_idx_t *const SMESH_RESTRICT e2s_idx) {
  element_idx_t *bk =
      (element_idx_t *)calloc(n_elements, sizeof(element_idx_t));

  for (ptrdiff_t i = 0; i < n_sides; i++) {
    element_idx_t e = parent_element[i];
    e2s_idx[e2s_ptr[e] + bk[e]++] = i;
  }

  free(bk);
  return SMESH_SUCCESS;
}

template <typename element_idx_t, typename count_t, typename idx_t,
          typename mask_t, typename selector_t>
int sideset_select_propagate(
    // Sideset
    const ptrdiff_t n_sides,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const i16 *const SMESH_RESTRICT side_idx,
    const count_t *const SMESH_RESTRICT n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx,
    // Mesh connectivity
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    // Selection
    const element_idx_t sideset_seed, mask_t *const SMESH_RESTRICT selected,
    selector_t &&selector) {

  mask_t *visited = (mask_t *)calloc(mask_count(n_sides), sizeof(mask_t));

  memset(selected, 0, mask_count(n_sides) * sizeof(mask_t));
  element_idx_t *e2s_ptr =
      (element_idx_t *)malloc((n_elements + 1) * sizeof(element_idx_t));
  sideset_to_e2s_count(n_sides, n_elements, parent_element, e2s_ptr);

  element_idx_t *e2s_idx =
      (element_idx_t *)malloc(e2s_ptr[n_elements] * sizeof(element_idx_t));
  sideset_to_e2s_fill(n_sides, n_elements, parent_element, e2s_ptr, e2s_idx);

  element_idx_t *queue =
      (element_idx_t *)malloc(n_sides * sizeof(element_idx_t));
  ptrdiff_t queue_size = 0;
  queue[queue_size++] = sideset_seed;

  LocalSideTable lst;
  lst.fill(element_type);

  const int nnxs = elem_num_nodes(side_type(element_type));
  const enum ElemType st = side_type(element_type);
  const int side_dim = elem_manifold_dim(st);
  const int min_shared = (side_dim == 1) ? 1 : 2;

  for (ptrdiff_t cursor = 0; cursor < queue_size; cursor++) {
    const ptrdiff_t side = queue[cursor];
    mask_set(side, selected);

    const element_idx_t e = parent_element[side];
    const i16 s = side_idx[side];
    for (int vn = 0; vn < nnxs; ++vn) {
      const idx_t vtx = elements[lst(s, vn)][e];
      const count_t beg = n2e_ptr[vtx];
      const count_t end = n2e_ptr[vtx + 1];

      for (count_t it = beg; it < end; ++it) {
        const element_idx_t esp = n2e_idx[it];

        for (ptrdiff_t k = e2s_ptr[esp]; k < e2s_ptr[esp + 1]; ++k) {
          const element_idx_t nside = e2s_idx[k];
          if (nside == side || nside == invalid_idx<element_idx_t>() ||
              mask_get(nside, visited) || mask_get(nside, selected))
            continue;

          const i16 ns = side_idx[nside];
          mask_set(nside, visited);

          // count shared nodes (edge adjacency requires >=2)
          int shared = 0;
          for (int a = 0; a < nnxs; ++a) {
            for (int b = 0; b < nnxs; ++b) {
              shared += (elements[lst(s, a)][e] == elements[lst(ns, b)][esp]);
            }
          }
          if (shared < min_shared)
            continue;

          
          if (selector(side, nside)) {
            queue[queue_size++] = nside;
          }
        }
      }
    }
  }

  free(queue);
  free(e2s_ptr);
  free(e2s_idx);
  free(visited);
  return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int extract_nodeset_from_sideset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT side_idx, ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out) {
  const enum ElemType st = side_type(element_type);
  const int nn = elem_num_nodes(st);

  const ptrdiff_t n = nn * n_surf_elements;
  idx_t *nodes = (idx_t *)malloc(n * sizeof(idx_t));
  LocalSideTable lst;
  lst.fill(element_type);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n_surf_elements; i++) {
    const ptrdiff_t e = parent_element[i];
    const int s = side_idx[i];

    for (int k = 0; k < nn; k++) {
      idx_t node = elems[lst(s, k)][e];
      nodes[i * nn + k] = node;
    }
  }

  *n_nodes_out = (ptrdiff_t)sort_and_unique(nodes, (size_t)n);
  *nodes_out = (idx_t *)realloc(nodes, (size_t)(*n_nodes_out) * sizeof(idx_t));

  return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int extract_nodeset_from_sidesets(
    ptrdiff_t n_sidesets, const enum ElemType element_type[],
    idx_t **const SMESH_RESTRICT elems[], const ptrdiff_t n_surf_elements[],
    const element_idx_t *const SMESH_RESTRICT parent_element[],
    const int16_t *const SMESH_RESTRICT side_idx[], ptrdiff_t *n_nodes_out,
    idx_t **SMESH_RESTRICT nodes_out) {
  ptrdiff_t n_nodes = 0;
  for (ptrdiff_t ss = 0; ss < n_sidesets; ss++) {
    const enum ElemType st = side_type(element_type[ss]);
    const int nn = elem_num_nodes(st);
    n_nodes += n_surf_elements[ss] * nn;
  }

  idx_t *nodes = (idx_t *)malloc((size_t)n_nodes * sizeof(idx_t));
  ptrdiff_t node_offset = 0;

  for (ptrdiff_t ss = 0; ss < n_sidesets; ss++) {
    const enum ElemType st = side_type(element_type[ss]);
    const int nn = elem_num_nodes(st);

    const ptrdiff_t n = nn * n_surf_elements[ss];
    LocalSideTable lst;
    lst.fill(element_type[ss]);

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n_surf_elements[ss]; i++) {
      const ptrdiff_t e = parent_element[ss][i];
      const int s = side_idx[ss][i];

      for (int k = 0; k < nn; k++) {
        idx_t node = elems[ss][lst(s, k)][e];
        nodes[node_offset + i * nn + k] = node;
      }
    }

    node_offset += n;
  }

  *n_nodes_out = (ptrdiff_t)sort_and_unique(nodes, (size_t)n_nodes);
  *nodes_out = (idx_t *)realloc(nodes, (size_t)(*n_nodes_out) * sizeof(idx_t));

  return SMESH_SUCCESS;
}

template <typename idx_t, typename element_idx_t>
int extract_surface_from_sideset(
    const enum ElemType element_type,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
    const ptrdiff_t n_surf_elements,
    const element_idx_t *const SMESH_RESTRICT parent_element,
    const int16_t *const SMESH_RESTRICT side_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT sides) {
  LocalSideTable lst;
  lst.fill(element_type);
  const int nn = elem_num_nodes(side_type(element_type));

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n_surf_elements; i++) {
    const ptrdiff_t e = parent_element[i];
    const int s = side_idx[i];

    for (int n = 0; n < nn; n++) {
      idx_t node = elems[lst(s, n)][e];
      sides[n][i] = node;
    }
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_SIDESETS_IMPL_HPP
