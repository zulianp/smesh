#ifndef SMESH_PROMOTIONS_IMPL_HPP
#define SMESH_PROMOTIONS_IMPL_HPP

#include "smesh_promotions.hpp"
#include "smesh_types.hpp"
#include "smesh_search.hpp"

namespace smesh {


template <typename idx_t, typename count_t, typename element_idx_t>
void mesh_tet4_to_tet15(
    const ptrdiff_t n_elements, const ptrdiff_t n_nodes,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements,
    const count_t *const SMESH_RESTRICT n2n_upper_triangular_ptr,
    const idx_t *const SMESH_RESTRICT n2n_upper_triangular_idx,
    const element_idx_t *const SMESH_RESTRICT e2e_table,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_elements,
    ptrdiff_t *const SMESH_RESTRICT tet15_n_nodes) {
  const int nsxe = 4;
  const ptrdiff_t id_offset = n_nodes + n2n_upper_triangular_ptr[n_nodes];
  ptrdiff_t n_unique_sides = 0;

  idx_t *side_ids = (idx_t *)malloc(n_elements * nsxe * sizeof(idx_t));
  {
#ifndef NDEBUG
    for (ptrdiff_t e = 0; e < n_elements * 4; e++) {
      side_ids[e] = invalid_idx<element_idx_t>();
    }
#endif

    for (ptrdiff_t e = 0; e < n_elements; e++) {
      for (int s = 0; s < nsxe; s++) {
        const ptrdiff_t offset = e * nsxe + s;
        const element_idx_t neigh = e2e_table[offset];
        if (neigh == invalid_idx<element_idx_t>() || e < neigh) {
          // Create new id
          side_ids[offset] = id_offset + n_unique_sides;
          n_unique_sides++;
        } else if (neigh != invalid_idx<element_idx_t>() && e > neigh) {
          // Search id in neighbor
          bool found = false;
          for (int sneigh = 0; sneigh < nsxe; sneigh++) {
            const ptrdiff_t offset_neigh = neigh * nsxe + sneigh;
            if (e2e_table[offset_neigh] == e) {
              SMESH_ASSERT(side_ids[offset_neigh] !=
                           invalid_idx<element_idx_t>());
              side_ids[offset] = side_ids[offset_neigh];
              found = true;
              break;
            }
          }

          SMESH_ASSERT(found);
          SMESH_UNUSED(found);
        }
      }
    }
  }

  const ptrdiff_t n_edges = n2n_upper_triangular_ptr[n_nodes];
  {
    // Element connectivity
    for (ptrdiff_t e = 0; e < n_elements; e++) {
      // Nodes
      idx_t ii[15] = {tet4_elements[0][e], tet4_elements[1][e],
                      tet4_elements[2][e], tet4_elements[3][e]};

      // Edges
      {
        idx_t row[6];
        row[0] = std::min(tet4_elements[0][e], tet4_elements[1][e]);
        row[1] = std::min(tet4_elements[1][e], tet4_elements[2][e]);
        row[2] = std::min(tet4_elements[0][e], tet4_elements[2][e]);
        row[3] = std::min(tet4_elements[0][e], tet4_elements[3][e]);
        row[4] = std::min(tet4_elements[1][e], tet4_elements[3][e]);
        row[5] = std::min(tet4_elements[2][e], tet4_elements[3][e]);

        idx_t key[6];
        key[0] = std::max(tet4_elements[0][e], tet4_elements[1][e]);
        key[1] = std::max(tet4_elements[1][e], tet4_elements[2][e]);
        key[2] = std::max(tet4_elements[0][e], tet4_elements[2][e]);
        key[3] = std::max(tet4_elements[0][e], tet4_elements[3][e]);
        key[4] = std::max(tet4_elements[1][e], tet4_elements[3][e]);
        key[5] = std::max(tet4_elements[2][e], tet4_elements[3][e]);

        for (int l = 0; l < 6; l++) {
          const idx_t r = row[l];
          const count_t row_begin = n2n_upper_triangular_ptr[r];
          const count_t len_row = n2n_upper_triangular_ptr[r + 1] - row_begin;
          const idx_t *cols = &n2n_upper_triangular_idx[row_begin];
          const count_t k =
              binary_search<idx_t, count_t>(key[l], cols, len_row);
          SMESH_ASSERT(k < len_row);
          ii[4 + l] = row_begin + k + n_nodes;
        }
      }

      // Faces
      {
        ii[10] = side_ids[e * nsxe + 0];
        ii[11] = side_ids[e * nsxe + 1];
        ii[12] = side_ids[e * nsxe + 2];
        ii[13] = side_ids[e * nsxe + 3];

        SMESH_ASSERT(ii[10] >= n_nodes + n_edges);
        SMESH_ASSERT(ii[11] >= n_nodes + n_edges);
        SMESH_ASSERT(ii[12] >= n_nodes + n_edges);
        SMESH_ASSERT(ii[13] >= n_nodes + n_edges);
      }

      // Volume
      ii[14] = n_nodes + n_edges + n_unique_sides + e;

      for (int node = 0; node < 15; node++) {
        tet15_elements[node][e] = ii[node];
      }
    }
  }

  *tet15_n_nodes = n_nodes + n_edges + n_unique_sides + n_elements;
  free(side_ids);
}

template <typename idx_t, typename count_t, typename geom_t>
void mesh_tet4_to_tet15_points(
    const ptrdiff_t n_elements, const ptrdiff_t tet4_n_nodes,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_pts,
    const count_t *const SMESH_RESTRICT n2n_upper_triangular_ptr,
    const idx_t *const SMESH_RESTRICT n2n_upper_triangular_idx,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_pts) {

  const ptrdiff_t n_edges = n2n_upper_triangular_ptr[tet4_n_nodes];
  const idx_t id_offset = idx_t(tet4_n_nodes + n_edges);
  SMESH_UNUSED(id_offset);

  // Edge points
  for (ptrdiff_t i = 0; i < tet4_n_nodes; i++) {
    const count_t row_begin = n2n_upper_triangular_ptr[i];
    const count_t len_row = n2n_upper_triangular_ptr[i + 1] - row_begin;
    const idx_t *cols = &n2n_upper_triangular_idx[row_begin];

    for (int d = 0; d < 3; d++) {
      tet15_pts[d][i] = tet4_pts[d][i];
    }

    for (count_t k = 0; k < len_row; k++) {
      const idx_t j = cols[k];
      const idx_t edge = tet4_n_nodes + row_begin + k;

      for (int d = 0; d < 3; d++) {
        tet15_pts[d][edge] = (tet4_pts[d][i] + tet4_pts[d][j]) / 2;
      }
    }
  }

  // Rest of the nodes
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    // Collect indices
    // Nodes
    idx_t ii[15] = {tet15_elements[0][e], tet15_elements[1][e],
                    tet15_elements[2][e], tet15_elements[3][e]};

    // Faces
    {
      ii[10] = tet15_elements[10][e];
      ii[11] = tet15_elements[11][e];
      ii[12] = tet15_elements[12][e];
      ii[13] = tet15_elements[13][e];

      SMESH_ASSERT(ii[10] >= id_offset);
      SMESH_ASSERT(ii[11] >= id_offset);
      SMESH_ASSERT(ii[12] >= id_offset);
      SMESH_ASSERT(ii[13] >= id_offset);
    }

    // Volume
    ii[14] = tet15_elements[14][e];

    // Compute points
    for (int d = 0; d < 3; d++) {
      tet15_pts[d][ii[10]] =
          (tet4_pts[d][ii[0]] + tet4_pts[d][ii[1]] + tet4_pts[d][ii[3]]) / 3;
      tet15_pts[d][ii[11]] =
          (tet4_pts[d][ii[1]] + tet4_pts[d][ii[2]] + tet4_pts[d][ii[3]]) / 3;
      tet15_pts[d][ii[12]] =
          (tet4_pts[d][ii[0]] + tet4_pts[d][ii[3]] + tet4_pts[d][ii[2]]) / 3;
      tet15_pts[d][ii[13]] =
          (tet4_pts[d][ii[0]] + tet4_pts[d][ii[2]] + tet4_pts[d][ii[1]]) / 3;
    }

    for (int d = 0; d < 3; d++) {
      tet15_pts[d][ii[14]] = (tet4_pts[d][ii[0]] + tet4_pts[d][ii[1]] +
                              tet4_pts[d][ii[2]] + tet4_pts[d][ii[3]]) /
                             4;
    }
  }
}

} // namespace smesh

#endif // SMESH_PROMOTIONS_IMPL_HPP