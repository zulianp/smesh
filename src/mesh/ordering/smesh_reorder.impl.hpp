#ifndef SMESH_REORDER_IMPL_HPP
#define SMESH_REORDER_IMPL_HPP

#include "smesh_reorder.hpp"

namespace smesh {

template <typename idx_t, typename element_idx_t>
int mesh_block_reorder(
    const int nxe, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT in_elements,
    const element_idx_t *const SMESH_RESTRICT e2e_gather,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT out_elements) {
  if (in_elements == out_elements) {
    // In place use temporary buffer
    idx_t *temp = (idx_t *)malloc(n_elements * sizeof(idx_t));
    for (int d = 0; d < nxe; ++d) {
      memcpy(temp, out_elements[d], n_elements * sizeof(idx_t));
      for (ptrdiff_t i = 0; i < n_elements; i++) {
        out_elements[d][i] = temp[e2e_gather[i]];
      }
    }

    free(temp);
  } else {
    for (int d = 0; d < nxe; ++d) {
      for (ptrdiff_t i = 0; i < n_elements; ++i) {
        out_elements[d][i] = in_elements[d][e2e_gather[i]];
      }
    }
  }
  return SMESH_SUCCESS;
}

template <typename idx_t>
int mesh_block_renumber_element_nodes(
    const int nxe, const ptrdiff_t n_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *const SMESH_RESTRICT next_node_idx,
    idx_t *const SMESH_RESTRICT *n2n_scatter) {
  for (ptrdiff_t e = 0; e < n_elements; ++e) {
    for (int d = 0; d < nxe; ++d) {
      idx_t ii = elements[d][e];

      if (n2n_scatter[ii] == invalid_idx<idx_t>()) {
        n2n_scatter[ii] = (*next_node_idx)++;
      }

      elements[d][e] = n2n_scatter[ii];
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int reorder_scatter(const ptrdiff_t n, const idx_t *const SMESH_RESTRICT *scatter,
                    T *const SMESH_RESTRICT in_array,
                    T *const SMESH_RESTRICT out_array) {
  if (in_array == out_array) {
    SMESH_ERROR("In place reordering of arrays is not supported");
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < n; ++i) {
    out_array[scatter[i]] = in_array[i];
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename T>
int reorder_gather(const ptrdiff_t n, 
                   const idx_t *const SMESH_RESTRICT *gather,
                   T *const SMESH_RESTRICT in_array,
                   T *const SMESH_RESTRICT out_array) {
  if (in_array == out_array) {
    SMESH_ERROR("In place reordering of arrays is not supported");
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < n; ++i) {
    out_array[i] = in_array[gather[i]];
  }

  return SMESH_SUCCESS;
}
} // namespace smesh

#endif // SMESH_REORDER_IMPL_HPP
