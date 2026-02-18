#ifndef SMESH_DECOMPOSE_HPP
#define SMESH_DECOMPOSE_HPP

#include "smesh_base.hpp"

#include <mpi.h>

namespace smesh {

inline ptrdiff_t rank_split(const ptrdiff_t n, const int comm_size,
                            const int comm_rank) {
  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t nlocal = uniform_split;
  ptrdiff_t remainder = n - nlocal * comm_size;

  if (remainder > comm_rank) {
    nlocal += 1;
  }

  return nlocal;
}

inline ptrdiff_t rank_start(const ptrdiff_t n, const int comm_size,
                            const int comm_rank) {
  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t remainder = n - uniform_split * comm_size;

  ptrdiff_t rank = comm_rank;
  ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);
  return rank_start;
}

inline int rank_owner(const ptrdiff_t n, const ptrdiff_t gidx,
                      const int comm_size) {
  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t remainder = n - uniform_split * comm_size;

  ptrdiff_t rank = gidx / uniform_split;
  ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);

  if (gidx >= rank_start) {
#ifndef NDEBUG
    ptrdiff_t rank_end =
        rank_start + uniform_split + (ptrdiff_t)(rank < remainder);
    assert(gidx < rank_end);
#endif
    return rank;
  } else {
    rank -= 1;
#ifndef NDEBUG
    ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);
    ptrdiff_t rank_end =
        rank_start + uniform_split + (ptrdiff_t)(rank < remainder);

    assert(gidx >= rank_start);
    assert(gidx < rank_end);
#endif
    return rank;
  }
}

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2e(MPI_Comm comm, const ptrdiff_t n_local_elements,
               const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
               const ptrdiff_t n_global_nodes, const int nnodesxelem,
               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
               count_t **out_n2eptr, element_idx_t **out_n2e_idx);

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2n_from_n2e(
    MPI_Comm comm, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
    const ptrdiff_t n_global_nodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx, count_t **out_n2n_ptr,
    element_idx_t **out_n2n_idx);

template <typename idx_t, typename count_t, typename element_idx_t>
int create_n2n_from_n2e(
    MPI_Comm comm, const ptrdiff_t n_local_elements,
    const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
    const ptrdiff_t n_global_nodes, const int nnodesxelem,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const count_t *const SMESH_RESTRICT n2eptr,
    const element_idx_t *const SMESH_RESTRICT n2e_idx, count_t **out_n2n_ptr,
    element_idx_t **out_n2n_idx);

template <typename idx_t, typename count_t, typename element_idx_t>
int redistribute_n2e(MPI_Comm comm, const int comm_size, const int comm_rank,
                     const ptrdiff_t n_local2global,
                     const ptrdiff_t n_global_nodes,
                     const ptrdiff_t n_global_elements,
                     const count_t *const SMESH_RESTRICT n2eptr,
                     const element_idx_t *const SMESH_RESTRICT n2e_idx,
                     ptrdiff_t *const SMESH_RESTRICT out_local2global_size,
                     idx_t **const SMESH_RESTRICT out_local2global,
                     count_t **const SMESH_RESTRICT out_local_n2e_ptr,
                     element_idx_t **const SMESH_RESTRICT out_local_n2e_idx);

template <typename idx_t, typename count_t, typename element_idx_t>
int localize_element_indices(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    idx_t *const *const SMESH_RESTRICT elems, const ptrdiff_t local2global_size,
    const count_t *const SMESH_RESTRICT local_n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const idx_t *const SMESH_RESTRICT local2global,
    idx_t **const SMESH_RESTRICT local_elements);

template <typename idx_t, typename count_t, typename element_idx_t>
int rearrange_local_nodes(const int comm_size, const int comm_rank,
                          const ptrdiff_t n_global_elements,
                          const ptrdiff_t n_local_elements,
                          const int nnodesxelem,
                          const ptrdiff_t local2global_size,
                          count_t *const SMESH_RESTRICT local_n2e_ptr,
                          element_idx_t *const SMESH_RESTRICT local_n2e_idx,
                          idx_t *const SMESH_RESTRICT local2global,
                          idx_t **const SMESH_RESTRICT local_elements,
                          ptrdiff_t *const SMESH_RESTRICT out_n_owned,
                          ptrdiff_t *const SMESH_RESTRICT out_n_shared,
                          ptrdiff_t *const SMESH_RESTRICT out_n_ghosts);

template <typename idx_t, typename count_t, typename element_idx_t>
int rearrange_local_elements(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    const ptrdiff_t local2global_size,
    count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    idx_t **const SMESH_RESTRICT local_elements, const ptrdiff_t n_owned_nodes,
    ptrdiff_t *const SMESH_RESTRICT n_owned_not_shared,
    element_idx_t *const SMESH_RESTRICT element_local_to_global);
} // namespace smesh

#endif // SMESH_DECOMPOSE_HPP
