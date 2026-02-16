#ifndef SMESH_DECOMPOSE_HPP
#define SMESH_DECOMPOSE_HPP

#include "smesh_base.hpp"

#include <mpi.h>

namespace smesh {

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(MPI_Comm comm, const ptrdiff_t n_local_elements,
               const ptrdiff_t n_global_elements,
               const ptrdiff_t n_local_nodes,
               const ptrdiff_t n_global_nodes, const int nnodesxelem,
               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
               count_t **out_n2eptr, element_idx_t **out_n2e_idx);

} // namespace smesh

#endif // SMESH_DECOMPOSE_HPP
