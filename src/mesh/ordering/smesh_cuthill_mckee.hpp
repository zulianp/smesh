#ifndef SMESH_CUTHILL_MCKEE_HPP
#define SMESH_CUTHILL_MCKEE_HPP

#include "smesh_base.hpp"

namespace smesh {

template <typename count_t, typename idx_t>
int eccentricity(const ptrdiff_t n_nodes,
                 const count_t *const SMESH_RESTRICT n2n_rowptr,
                 const idx_t *const SMESH_RESTRICT n2n_idx,
                 idx_t *const SMESH_RESTRICT out);

template <typename count_t, typename idx_t>
int cuthill_mckee(const ptrdiff_t n_nodes,
                  const count_t *const SMESH_RESTRICT n2n_rowptr,
                  const idx_t *const SMESH_RESTRICT n2n_idx,
                  idx_t *const SMESH_RESTRICT reordering);
} // namespace smesh

#endif // SMESH_CUTHILL_MCKEE_HPP