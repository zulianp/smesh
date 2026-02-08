#ifndef SMESH_SEARCH_HPP
#define SMESH_SEARCH_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {
template <typename idx_t>
static SMESH_INLINE idx_t linear_search(const idx_t target,
                                        const idx_t *const arr,
                                        const int size) {
  idx_t i;
  for (i = 0; i < size - 4; i += 4) {
    if (arr[i] == target)
      return i;
    if (arr[i + 1] == target)
      return i + 1;
    if (arr[i + 2] == target)
      return i + 2;
    if (arr[i + 3] == target)
      return i + 3;
  }
  for (; i < size; i++) {
    if (arr[i] == target)
      return i;
  }
  return invalid_idx<idx_t>();
}

template <typename idx_t, typename count_t>
static inline count_t binary_search(const idx_t key, const idx_t *const cols,
                                    const count_t len_row) {
  count_t lo = 0;
  count_t hi = len_row;
  while (lo < hi) {
    const count_t mid = lo + (hi - lo) / 2;
    const idx_t v = cols[mid];
    if (v < key) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

template <int N_TARGETS, typename idx_t>
static SMESH_INLINE void find(const idx_t *SMESH_RESTRICT targets,
                              const idx_t *const SMESH_RESTRICT row,
                              const int lenrow, idx_t *SMESH_RESTRICT ks) {
#pragma unroll(N_TARGETS)
  for (int d = 0; d < N_TARGETS; ++d) {
    ks[d] = 0;
  }

  for (int i = 0; i < lenrow; ++i) {
#pragma unroll(N_TARGETS)
    for (int d = 0; d < N_TARGETS; ++d) {
      ks[d] += row[i] < targets[d];
    }
  }
}

} // namespace smesh

#endif // SMESH_SEARCH_HPP