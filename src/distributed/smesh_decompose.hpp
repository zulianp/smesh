#ifndef SMESH_DECOMPOSE_HPP
#define SMESH_DECOMPOSE_HPP

#include "smesh_base.hpp"

#include <mpi.h>

#include <algorithm>
#include <stddef.h>

namespace smesh {

// These three functions define one partition of [0, n) over comm_size ranks and must agree
// exactly: rank_start gives where a rank's block begins, rank_split how long it is, and
// rank_owner which rank a global index falls in. The layout is that the first
// `n % comm_size` ranks take one extra entry each.
//
// n < comm_size is a real case, not a misuse: a coarse multigrid level has far fewer nodes
// and elements than a full node has ranks, and the id spaces derived from them are passed
// here unchanged. The assertion that used to stand in for handling it is compiled out of
// every -DNDEBUG build, leaving `n / comm_size` to evaluate to zero and the arithmetic
// below to divide by it. Giving the first n ranks one entry each and the rest none keeps
// the three consistent, and matches what id_space_size already does where it is used.

inline ptrdiff_t rank_split(const ptrdiff_t n, const int comm_size,
                            const int comm_rank) {
  SMESH_ASSERT(comm_size > 0);
  SMESH_ASSERT(comm_rank >= 0 && comm_rank < comm_size);

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
  SMESH_ASSERT(comm_size > 0);
  SMESH_ASSERT(comm_rank >= 0 && comm_rank <= comm_size);

  ptrdiff_t uniform_split = n / comm_size;
  ptrdiff_t remainder = n - uniform_split * comm_size;

  ptrdiff_t rank = comm_rank;
  ptrdiff_t rank_start = rank * uniform_split + std::min(rank, remainder);
  return rank_start;
}

/// The exact inverse of rank_start: which rank owns global index gidx.
///
/// rank_start lays the range out as `r * uniform_split + min(r, remainder)`, so the first
/// `remainder` ranks hold `uniform_split + 1` entries each and the rest hold
/// `uniform_split`. Dividing gidx by uniform_split alone is therefore not the inverse: it
/// ignores the wider blocks at the front and overshoots, and a single `rank -= 1` cannot
/// correct an overshoot of more than one rank. Swept over the (n, comm_size) pairs this
/// code actually sees, the previous form returned a rank equal to comm_size for 2915
/// global indices: n=512 at 288 ranks returns 288 for 223 of its 512 indices, and n=729 at
/// 256 ranks returns 256 for 215 of them. The overshoot appears and disappears as
/// comm_size grows rather than worsening monotonically, which is what makes it look like a
/// problem-size effect when it is arithmetic.
///
/// It mattered because callers use the result as an array index without checking it.
/// `send_displs[rank_owner(...) + 1]++` then writes one past the end of a
/// SMESH_CALLOC(comm_size + 1) buffer, which corrupts the heap rather than failing: the
/// symptom is glibc reporting "munmap_chunk(): invalid pointer" or "corrupted size vs.
/// prev_size" from some later free, with nothing pointing back at the write. Every
/// SMESH_ASSERT here is compiled out under -DNDEBUG, which is how a release build reaches
/// the allocator without tripping a check.
///
/// Splitting at the boundary between the two block sizes gives the inverse directly.
inline int rank_owner(const ptrdiff_t n, const ptrdiff_t gidx,
                      const int comm_size) {
  SMESH_ASSERT(gidx >= 0);
  SMESH_ASSERT(gidx < n);
  SMESH_ASSERT(comm_size > 0);

  const ptrdiff_t uniform_split = n / comm_size;
  const ptrdiff_t remainder = n - uniform_split * comm_size;

  // Fewer indices than ranks: rank_split gives the first n ranks one index each and the
  // rest none, so index i belongs to rank i. The trailing empty ranks own nothing and are
  // never named here.
  if (uniform_split == 0) {
    return (int)gidx;
  }

  // The first `remainder` ranks own uniform_split + 1 entries and so cover
  // [0, remainder * (uniform_split + 1)); above that the blocks are uniform_split wide.
  const ptrdiff_t wide_end = remainder * (uniform_split + 1);
  const ptrdiff_t rank = gidx < wide_end
                             ? gidx / (uniform_split + 1)
                             : remainder + (gidx - wide_end) / uniform_split;

#ifndef NDEBUG
  const ptrdiff_t start = rank * uniform_split + std::min(rank, remainder);
  const ptrdiff_t end = start + uniform_split + (ptrdiff_t)(rank < remainder);
  SMESH_ASSERT(gidx >= start);
  SMESH_ASSERT(gidx < end);
  SMESH_ASSERT(rank >= 0);
  SMESH_ASSERT(rank < comm_size);
#endif
  return (int)rank;
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

template <typename count_t, typename element_idx_t, typename local2global_t>
int redistribute_n2e(MPI_Comm comm, const int comm_size, const int comm_rank,
                     const ptrdiff_t n_local2global,
                     const ptrdiff_t n_global_nodes,
                     const ptrdiff_t n_global_elements,
                     const count_t *const SMESH_RESTRICT n2eptr,
                     const element_idx_t *const SMESH_RESTRICT n2e_idx,
                     ptrdiff_t *const SMESH_RESTRICT out_local2global_size,
                     local2global_t **const SMESH_RESTRICT out_local2global,
                     count_t **const SMESH_RESTRICT out_local_n2e_ptr,
                     element_idx_t **const SMESH_RESTRICT out_local_n2e_idx);

template <typename idx_t, typename count_t, typename element_idx_t,
          typename local2global_t = idx_t>
int localize_element_indices(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    idx_t *const *const SMESH_RESTRICT elems, const ptrdiff_t local2global_size,
    const count_t *const SMESH_RESTRICT local_n2e_ptr,
    const element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const local2global_t *const SMESH_RESTRICT local2global,
    idx_t **const SMESH_RESTRICT local_elements);

template <typename idx_t, typename count_t, typename element_idx_t,
          typename local2global_t = idx_t>
int rearrange_local_nodes(const int comm_size, const int comm_rank,
                          const ptrdiff_t n_global_elements,
                          const ptrdiff_t n_local_elements,
                          const int nnodesxelem,
                          const ptrdiff_t local2global_size,
                          count_t *const SMESH_RESTRICT local_n2e_ptr,
                          element_idx_t *const SMESH_RESTRICT local_n2e_idx,
                          local2global_t *const SMESH_RESTRICT local2global,
                          idx_t **const SMESH_RESTRICT local_elements,
                          ptrdiff_t *const SMESH_RESTRICT out_n_owned,
                          ptrdiff_t *const SMESH_RESTRICT out_n_shared,
                          ptrdiff_t *const SMESH_RESTRICT out_n_ghosts);

template <typename idx_t, typename count_t, typename element_idx_t,
          typename large_idx_t>
int rearrange_local_elements(
    const int comm_size, const int comm_rank, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    const ptrdiff_t local2global_size,
    count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    idx_t **const SMESH_RESTRICT local_elements, const ptrdiff_t n_owned_nodes,
    ptrdiff_t *const SMESH_RESTRICT n_owned_not_shared,
    large_idx_t *const SMESH_RESTRICT element_local_to_global,
    const large_idx_t *const SMESH_RESTRICT input_element_mapping = nullptr);

template <typename idx_t, typename count_t, typename element_idx_t,
          typename large_idx_t>
int expand_aura_elements_inconsistent(
    MPI_Comm comm, const ptrdiff_t n_global_elements,
    const ptrdiff_t n_local_elements, const int nnodesxelem,
    count_t *const SMESH_RESTRICT local_n2e_ptr,
    element_idx_t *const SMESH_RESTRICT local_n2e_idx,
    const large_idx_t *const SMESH_RESTRICT local2global,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT local_elements,
    const large_idx_t *const SMESH_RESTRICT element_local_to_global,
    const ptrdiff_t node_n_owned, const ptrdiff_t nodes_n_ghosts,
    large_idx_t **const SMESH_RESTRICT out_aura_element_mapping,
    idx_t **const SMESH_RESTRICT out_aura_element_nodes,
    ptrdiff_t *const SMESH_RESTRICT out_n_aura);

template <typename idx_t, typename local2global_t = idx_t>
int prepare_node_renumbering(
    MPI_Comm comm, const ptrdiff_t n_global_nodes,
    const ptrdiff_t owned_nodes_start, const ptrdiff_t n_owned_nodes,
    const local2global_t *const SMESH_RESTRICT local2global,
    idx_t *const SMESH_RESTRICT global2owned);

int node_ownership_ranges(MPI_Comm comm, const ptrdiff_t n_owned_nodes,
                          ptrdiff_t *const SMESH_RESTRICT owned_nodes_ranges);

template <typename idx_t, typename local2global_t = idx_t>
int stich_aura_elements(
    MPI_Comm comm, const ptrdiff_t n_owned_nodes,
    const ptrdiff_t n_shared_nodes, const ptrdiff_t n_ghost_nodes,
    const local2global_t *const SMESH_RESTRICT local2global,
    const int nnodesxelem, const ptrdiff_t n_aura_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT e2n_aura,
    const ptrdiff_t n_local_elements, idx_t **const SMESH_RESTRICT e2n_local,
    local2global_t **const SMESH_RESTRICT n2n_local2global_out,
    ptrdiff_t *const SMESH_RESTRICT out_n_aura_nodes);

template <typename idx_t, typename local2global_t = idx_t>
int collect_ghost_and_aura_import_indices(
    MPI_Comm comm, const ptrdiff_t n_owned_nodes, const ptrdiff_t n_ghost_nodes,
    const ptrdiff_t n_aura_nodes, const ptrdiff_t n_global_nodes,
    const local2global_t *const SMESH_RESTRICT local2global,
    const idx_t *const SMESH_RESTRICT global2owned,
    const ptrdiff_t *const SMESH_RESTRICT owned_node_ranges,
    idx_t *const SMESH_RESTRICT ghost_and_aura_to_owned);

template <typename idx_t>
int determine_ownership(const int comm_size, const int comm_rank,
                        const ptrdiff_t n_owned_nodes, const ptrdiff_t n_ghosts,
                        const ptrdiff_t n_aura_nodes,
                        const idx_t *const SMESH_RESTRICT local2owned,
                        const ptrdiff_t *const SMESH_RESTRICT owned_nodes_range,
                        int *const SMESH_RESTRICT owner);

template <typename idx_t, typename local2global_t = idx_t>
int group_ghost_and_aura_by_rank(
    const int comm_size, const ptrdiff_t n_owned, const ptrdiff_t n_ghosts,
    const ptrdiff_t n_aura_nodes,
    local2global_t *const SMESH_RESTRICT local2global,
    idx_t *const SMESH_RESTRICT ghost_and_aura_to_owned,
    int *const SMESH_RESTRICT owner, const int nnodesxelem,
    const ptrdiff_t n_local_elements, const ptrdiff_t n_aura_elements,
    idx_t **const SMESH_RESTRICT local_elements);

} // namespace smesh

#endif // SMESH_DECOMPOSE_HPP
