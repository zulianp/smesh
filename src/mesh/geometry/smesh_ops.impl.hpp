#ifndef SMESH_OPS_IMPL_HPP
#define SMESH_OPS_IMPL_HPP

#include "smesh_ops.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
int barycenters(
    const int nxe, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int spatial_dim,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT barycenters) {

#pragma omp parallel for
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    geom_t barycenter[3];
    for (int d = 0; d < spatial_dim; d++) {
      barycenter[d] = 0;
    }

    for (int d = 0; d < spatial_dim; d++) {
      for (int i = 0; i < nxe; i++) {
        barycenter[d] += points[d][elements[i][e]];
      }

      barycenter[d] /= nxe;
    }

    for (int d = 0; d < 3; d++) {
      barycenters[d][e] = barycenter[d];
    }
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_OPS_IMPL_HPP