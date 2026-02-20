#ifndef SMESH_REFINE_IMPL_HPP
#define SMESH_REFINE_IMPL_HPP

#include "smesh_elem_type.hpp"
#include "smesh_refine.hpp"
#include "smesh_search.hpp"

namespace smesh {

static const int tet4_refine_pattern[8][4] = {
    // Corner tests
    {0, 4, 6, 7},
    {4, 1, 5, 8},
    {6, 5, 2, 9},
    {7, 8, 9, 3},
    // Octahedron tets
    {4, 5, 6, 8},
    {7, 4, 6, 8},
    {6, 5, 9, 8},
    {7, 6, 9, 8}};

static const int tri3_refine_pattern[4][3] = {
    // Corner triangles
    {0, 3, 5},
    {3, 1, 4},
    {5, 4, 2},
    // Center triangle
    {3, 4, 5}};

template <typename idx_t, typename count_t,
          typename geom_t>
int mesh_refine(
    const enum ElemType element_type, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT coarse_elements,
    const int spatial_dim, const ptrdiff_t n_coarse_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT coarse_points,
    const count_t *const SMESH_RESTRICT n2n_ptr,
    const idx_t *const SMESH_RESTRICT n2n_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT refined_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT refined_points) {
  if (element_type != TET4 && element_type != TRI3) {
    SMESH_ERROR("mesh_refine: unsupported element_type %d\n", element_type);
    return SMESH_FAILURE;
  }

  const count_t nnz = n2n_ptr[n_coarse_nodes];
  idx_t *edge_idx = (idx_t *)malloc(nnz * sizeof(idx_t));
  memset(edge_idx, 0, nnz * sizeof(idx_t));

  idx_t next_id = n_coarse_nodes;
  for (ptrdiff_t i = 0; i < n_coarse_nodes; i++) {
    const count_t begin = n2n_ptr[i];
    const count_t end = n2n_ptr[i + 1];

    for (count_t k = begin; k < end; k++) {
      const idx_t j = n2n_idx[k];

      if (i < j) {
        edge_idx[k] = next_id++;
      }
    }
  }

  const int n_nodes_per_element = elem_num_nodes(element_type);
  for (int d = 0; d < n_nodes_per_element; d++) {
    // Copy p1 portion
    memcpy(refined_elements[d], coarse_elements[d], n_elements * sizeof(idx_t));
  }

  for (int d = 0; d < spatial_dim; d++) {
    // Copy p1 portion
    memcpy(refined_points[d], coarse_points[d],
           n_coarse_nodes * sizeof(geom_t));
  }

  if (element_type == TET4) {
    // TODO fill p2 node indices in elements
    for (ptrdiff_t e = 0; e < n_elements; e++) {
      idx_t macro_element[10];
      for (int k = 0; k < 4; k++) {
        macro_element[k] = coarse_elements[k][e];
      }

      // Ordering of edges compliant to exodusII spec
      idx_t row[6];
      row[0] = std::min(macro_element[0], macro_element[1]);
      row[1] = std::min(macro_element[1], macro_element[2]);
      row[2] = std::min(macro_element[0], macro_element[2]);
      row[3] = std::min(macro_element[0], macro_element[3]);
      row[4] = std::min(macro_element[1], macro_element[3]);
      row[5] = std::min(macro_element[2], macro_element[3]);

      idx_t key[6];
      key[0] = std::max(macro_element[0], macro_element[1]);
      key[1] = std::max(macro_element[1], macro_element[2]);
      key[2] = std::max(macro_element[0], macro_element[2]);
      key[3] = std::max(macro_element[0], macro_element[3]);
      key[4] = std::max(macro_element[1], macro_element[3]);
      key[5] = std::max(macro_element[2], macro_element[3]);

      for (int l = 0; l < 6; l++) {
        const idx_t r = row[l];
        const count_t row_begin = n2n_ptr[r];
        const count_t len_row = n2n_ptr[r + 1] - row_begin;
        const idx_t *cols = &n2n_idx[row_begin];
        const idx_t k = binary_search(key[l], cols, len_row);
        macro_element[l + 4] = edge_idx[row_begin + k];
      }

      // distribute macro_element to fine_mesh
      ptrdiff_t element_offset = e * 8;
      for (int k = 0; k < 4; k++) {
        for (int sub_e = 0; sub_e < 8; sub_e++) {
          const idx_t ik = macro_element[tet4_refine_pattern[sub_e][k]];
          refined_elements[k][element_offset + sub_e] = ik;
        }
      }
    }

  } else if (element_type == TRI3) {
    // TODO fill p2 node indices in elements
    for (ptrdiff_t e = 0; e < n_elements; e++) {
      idx_t macro_element[6];
      for (int k = 0; k < 3; k++) {
        macro_element[k] = coarse_elements[k][e];
      }

      // Ordering of edges compliant to exodusII spec
      idx_t row[3];
      row[0] = std::min(macro_element[0], macro_element[1]);
      row[1] = std::min(macro_element[1], macro_element[2]);
      row[2] = std::min(macro_element[0], macro_element[2]);

      idx_t key[3];
      key[0] = std::max(macro_element[0], macro_element[1]);
      key[1] = std::max(macro_element[1], macro_element[2]);
      key[2] = std::max(macro_element[0], macro_element[2]);

      for (int l = 0; l < 3; l++) {
        const idx_t r = row[l];
        const count_t row_begin = n2n_ptr[r];
        const count_t len_row = n2n_ptr[r + 1] - row_begin;
        const idx_t *cols = &n2n_idx[row_begin];
        const idx_t k = binary_search(key[l], cols, len_row);
        macro_element[l + 3] = edge_idx[row_begin + k];
      }

      // distribute macro_element to fine_mesh
      ptrdiff_t element_offset = e * 4;
      for (int k = 0; k < 3; k++) {
        for (int sub_e = 0; sub_e < 4; sub_e++) {
          const idx_t ik = macro_element[tri3_refine_pattern[sub_e][k]];
          refined_elements[k][element_offset + sub_e] = ik;
        }
      }
    }
  }

     // Generate p2 coordinates
     for (ptrdiff_t i = 0; i < n_coarse_nodes; i++) {
      const count_t begin = n2n_ptr[i];
      const count_t end   = n2n_ptr[i + 1];

      for (count_t k = begin; k < end; k++) {
          const idx_t j = n2n_idx[k];

          if (i < j) {
              const idx_t nidx = edge_idx[k];

              for (int d = 0; d < spatial_dim; d++) {
                  geom_t xi = refined_points[d][i];
                  geom_t xj = refined_points[d][j];

                  // Midpoint
                  refined_points[d][nidx] = (xi + xj) / 2;
              }
          }
      }
  }

  free(edge_idx);
  return SMESH_SUCCESS;
}
} // namespace smesh

#endif // SMESH_REFINE_IMPL_HPP