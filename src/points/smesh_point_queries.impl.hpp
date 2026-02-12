#ifndef SMESH_POINT_QUERIES_IMPL_HPP
#define SMESH_POINT_QUERIES_IMPL_HPP

#include "smesh_point_queries.hpp"

#include <cstring>

namespace smesh {

template <typename T>
int sqdist(const int spatial_dim, const ptrdiff_t n_nodes,
           const T *const SMESH_RESTRICT *const SMESH_RESTRICT points,
           const T *const SMESH_RESTRICT point, T *const SMESH_RESTRICT *d) {

  memset(d, 0, n_nodes * sizeof(T));
  for (int d = 0; d < spatial_dim; ++d) {
#pragma omp parallel for
    for (ptrdiff_t node = 0; node < n_nodes; ++node) {
      const T diff = points[d][node] - point[d];
      d[node] += diff * diff;
    }
  }
  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_POINT_QUERIES_IMPL_HPP
