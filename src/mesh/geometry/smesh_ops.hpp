#ifndef SMESH_OPS_HPP
#define SMESH_OPS_HPP

#include "smesh_base.hpp"

namespace smesh {


template <typename idx_t, typename geom_t>
int barycenters(
    const int nxe, const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int spatial_dim,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT barycenters);

} // namespace smesh

#endif // SMESH_OPS_HPP
