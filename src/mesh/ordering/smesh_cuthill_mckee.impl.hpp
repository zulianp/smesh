#ifndef SMESH_CUTHILL_MCKEE_IMPL_HPP
#define SMESH_CUTHILL_MCKEE_IMPL_HPP

#include "smesh_cuthill_mckee.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <cmath>

namespace smesh {

template <typename count_t, typename idx_t>
int eccentricity(const ptrdiff_t n_nodes,
                 const count_t *const SMESH_RESTRICT n2n_rowptr,
                 const idx_t *const SMESH_RESTRICT n2n_idx,
                 idx_t *const SMESH_RESTRICT out) {
  if (n_nodes <= 0) {
    return SMESH_SUCCESS;
  }
  if (!n2n_rowptr || !n2n_idx || !out) {
    SMESH_ERROR("Invalid input");
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    out[i] = idx_t(0);
  }

  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    for (count_t j = n2n_rowptr[i]; j < n2n_rowptr[i + 1]; j++) {
      const ptrdiff_t neighbor = (ptrdiff_t)n2n_idx[j];
      if (neighbor != i) {
        const ptrdiff_t dist = std::abs(neighbor - i);
        if (dist > (ptrdiff_t)out[i]) {
          out[i] = (idx_t)dist;
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

//TODO: check for bugs
template <typename count_t, typename idx_t>
int cuthill_mckee(const ptrdiff_t n_nodes,
                  const count_t *const SMESH_RESTRICT n2n_rowptr,
                  const idx_t *const SMESH_RESTRICT n2n_idx,
                  idx_t *const SMESH_RESTRICT reordering) {
  if (n_nodes <= 0) {
    return SMESH_SUCCESS;
  }
  if (!n2n_rowptr || !n2n_idx || !reordering) {
    SMESH_ERROR("Invalid input");
    return SMESH_FAILURE;
  }

  idx_t *queue = (idx_t *)malloc(n_nodes * sizeof(idx_t));

  auto degree = [&n2n_rowptr](ptrdiff_t i) {
    return n2n_rowptr[i + 1] - n2n_rowptr[i];
  };

  // Initialize arrays
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    reordering[i] = invalid_idx<idx_t>();
  }

  // Find starting node with minimum degree
  ptrdiff_t start_node = 0;
  idx_t min_degree = degree(0);
  idx_t max_degree = min_degree;
  for (ptrdiff_t i = 1; i < n_nodes; i++) {
    idx_t d = degree(i);
    if (d < min_degree) {
      min_degree = d;
      start_node = i;
    } else if (d > max_degree) {
      max_degree = d;
    }
  }

  queue[0] = start_node;
  ptrdiff_t queue_size = 1;
  ptrdiff_t reorder_idx = 0;
  reordering[start_node] = reorder_idx++;

  ptrdiff_t cursor = 0;
  while (cursor < queue_size) {
    const ptrdiff_t current_queue_size = queue_size - cursor;
    const idx_t *const current_queue = &queue[cursor];

    ptrdiff_t end_cursor = cursor;
    for (ptrdiff_t i = 0; i < current_queue_size; i++) {
      const idx_t current = current_queue[i];
      if (current == invalid_idx<idx_t>()) {
        continue;
      }

      const ptrdiff_t start_cursor = end_cursor;
      for (ptrdiff_t j = n2n_rowptr[current]; j < n2n_rowptr[current + 1];
           j++) {
        const ptrdiff_t neighbor = n2n_idx[j];
        if (neighbor != current &&
            reordering[neighbor] == invalid_idx<idx_t>()) {
          queue[end_cursor++] = neighbor;
        }
      }

      std::sort(
          queue + start_cursor, queue + end_cursor,
          [&](ptrdiff_t a, ptrdiff_t b) { return degree(a) < degree(b); });

      for (ptrdiff_t i = start_cursor; i < end_cursor; i++) {
        const idx_t current = queue[i];
        if (current == invalid_idx<idx_t>()) {
          continue;
        }

        if (reordering[current] == invalid_idx<idx_t>()) {
          reordering[current] = reorder_idx++;
        } else {
          queue[i] = invalid_idx<idx_t>();
        }
      }
    }
  }

  // Handle any remaining unvisited nodes (disconnected components)
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    if (reordering[i] == invalid_idx<idx_t>()) {
      reordering[i] = reorder_idx++;
    }
  }

  free(queue);
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_CUTHILL_MCKEE_IMPL_HPP
