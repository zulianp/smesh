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

  const idx_t invalid = invalid_idx<idx_t>();
  const idx_t queued = invalid - idx_t(1);

  idx_t *queue = (idx_t *)malloc((size_t)n_nodes * sizeof(idx_t));
  count_t *degrees = (count_t *)malloc((size_t)n_nodes * sizeof(count_t));
  if (!queue || !degrees) {
    free(queue);
    free(degrees);
    SMESH_ERROR("Out of memory");
    return SMESH_FAILURE;
  }

  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    reordering[i] = invalid;
    degrees[i] = n2n_rowptr[i + 1] - n2n_rowptr[i];
  }

  // Find starting node with minimum degree
  ptrdiff_t start_node = 0;
  count_t min_degree = degrees[0];
  for (ptrdiff_t i = 1; i < n_nodes; i++) {
    const count_t d = degrees[i];
    if (d < min_degree) {
      min_degree = d;
      start_node = i;
    }
  }

  ptrdiff_t reorder_idx = 0;

  auto run_component = [&](const ptrdiff_t seed) {
    ptrdiff_t head = 0;
    ptrdiff_t tail = 0;

    queue[tail++] = (idx_t)seed;
    reordering[seed] = (idx_t)reorder_idx++;

    while (head < tail) {
      const ptrdiff_t current = (ptrdiff_t)queue[head++];
      SMESH_ASSERT(current >= 0 && current < n_nodes);

      const ptrdiff_t seg_begin = tail;
      for (count_t j = n2n_rowptr[current]; j < n2n_rowptr[current + 1]; j++) {
        const ptrdiff_t neighbor = (ptrdiff_t)n2n_idx[j];
        SMESH_ASSERT(neighbor >= 0 && neighbor < n_nodes);
        if (neighbor == current) {
          continue;
        }
        if (reordering[neighbor] == invalid) {
          reordering[neighbor] = queued; // mark discovered (prevents duplicates)
          queue[tail++] = (idx_t)neighbor;
          SMESH_ASSERT(tail <= n_nodes);
        }
      }

      if (tail - seg_begin > 1) {
        std::sort(queue + seg_begin, queue + tail, [&](idx_t a, idx_t b) {
          const count_t da = degrees[(ptrdiff_t)a];
          const count_t db = degrees[(ptrdiff_t)b];
          return (da < db) || (da == db && a < b);
        });
      }

      for (ptrdiff_t k = seg_begin; k < tail; k++) {
        const ptrdiff_t v = (ptrdiff_t)queue[k];
        reordering[v] = (idx_t)reorder_idx++;
      }
    }
  };

  run_component(start_node);

  // Handle any remaining unvisited nodes (disconnected components)
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    if (reordering[i] == invalid) {
      run_component(i);
    }
  }

  free(degrees);
  free(queue);
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_CUTHILL_MCKEE_IMPL_HPP
