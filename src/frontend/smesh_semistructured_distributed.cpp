#include "smesh_semistructured.hpp"

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_line_quadrature.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sstet4_mesh.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#ifdef SMESH_ENABLE_MPI
#include "smesh_alltoallv.impl.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#endif

namespace smesh {

#ifdef SMESH_ENABLE_MPI

namespace {

constexpr i64 k_alltoall_chunk = static_cast<i64>(1) << 20;
constexpr large_idx_t k_key_pad = static_cast<large_idx_t>(-1);

static void sort4(large_idx_t *v, const int n) {
    for (int i = 1; i < n; ++i) {
        const large_idx_t x = v[i];
        int               j = i;
        while (j > 0 && v[j - 1] > x) {
            v[j] = v[j - 1];
            --j;
        }
        v[j] = x;
    }
}

static ptrdiff_t id_space_size(const ptrdiff_t n, const int comm_size) {
    return n >= comm_size ? n : (ptrdiff_t)comm_size;
}

static int node_bucket(const int rank, const int owner, const int shared, const int used_owned, const int used_aura) {
    SMESH_UNUSED(used_aura);
    if (owner == rank) {
        return shared ? 1 : 0;
    }
    return used_owned ? 2 : 3;
}

/// Send tuples to rank_owner(n_space, key[0]). Unique per CRS row (short linear scan).
static int unique_tuples(MPI_Comm               comm,
                         const ptrdiff_t        n_space,
                         const ptrdiff_t        n_req,
                         const large_idx_t     *const keys,
                         const large_idx_t     *const aux,
                         large_idx_t           *const out_gid,
                         int                   *const out_owner,
                         int                   *const out_shared,
                         ptrdiff_t             *const n_global_unique) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const ptrdiff_t n_space_pad = id_space_size(n_space, size);
    const ptrdiff_t nodes_start = rank_start(n_space_pad, size, rank);
    const ptrdiff_t n_owned_space = rank_split(n_space_pad, size, rank);

    i64 *send_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
    i64 *send_count  = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
    i64 *recv_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
    i64 *recv_count  = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));

    for (ptrdiff_t i = 0; i < n_req; ++i) {
        const int p = rank_owner(n_space_pad, (ptrdiff_t)keys[i * 4], size);
        send_displs[p + 1]++;
    }

    SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count, 1, mpi_type<i64>(), comm));

    recv_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        recv_displs[r + 1] = recv_displs[r] + recv_count[r];
    }
    send_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] += send_displs[r];
    }

    const ptrdiff_t n_send = (ptrdiff_t)send_displs[size];
    const ptrdiff_t n_recv = (ptrdiff_t)recv_displs[size];
    constexpr int   k_nfields = 7;

    large_idx_t *send_pack = (large_idx_t *)SMESH_ALLOC((size_t)n_send * k_nfields * sizeof(large_idx_t));
    large_idx_t *recv_pack = (large_idx_t *)SMESH_ALLOC((size_t)n_recv * k_nfields * sizeof(large_idx_t));

    for (ptrdiff_t i = 0; i < n_req; ++i) {
        const int p    = rank_owner(n_space_pad, (ptrdiff_t)keys[i * 4], size);
        const i64 slot = send_displs[p] + send_count[p]++;
        large_idx_t *row = send_pack + (ptrdiff_t)slot * k_nfields;
        row[0] = keys[i * 4 + 0];
        row[1] = keys[i * 4 + 1];
        row[2] = keys[i * 4 + 2];
        row[3] = keys[i * 4 + 3];
        row[4] = aux[i];
        row[5] = (large_idx_t)rank;
        row[6] = (large_idx_t)i;
    }

    if (all_to_allv_64v(send_pack,
                        send_count,
                        send_displs,
                        recv_pack,
                        recv_count,
                        recv_displs,
                        k_nfields,
                        comm,
                        k_alltoall_chunk) != SMESH_SUCCESS) {
        SMESH_FREE(send_displs);
        SMESH_FREE(send_count);
        SMESH_FREE(recv_displs);
        SMESH_FREE(recv_count);
        SMESH_FREE(send_pack);
        SMESH_FREE(recv_pack);
        return SMESH_FAILURE;
    }

    ptrdiff_t *rowptr = (ptrdiff_t *)SMESH_CALLOC((size_t)n_owned_space + 1, sizeof(ptrdiff_t));
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const ptrdiff_t off = (ptrdiff_t)recv_pack[i * k_nfields] - nodes_start;
        rowptr[off + 1]++;
    }
    for (ptrdiff_t a = 0; a < n_owned_space; ++a) {
        rowptr[a + 1] += rowptr[a];
    }

    ptrdiff_t *order = (ptrdiff_t *)SMESH_ALLOC((size_t)n_recv * sizeof(ptrdiff_t));
    ptrdiff_t *book  = (ptrdiff_t *)SMESH_CALLOC((size_t)n_owned_space, sizeof(ptrdiff_t));
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const ptrdiff_t off = (ptrdiff_t)recv_pack[i * k_nfields] - nodes_start;
        order[rowptr[off] + book[off]++] = i;
    }

    ptrdiff_t   *uniq_of     = (ptrdiff_t *)SMESH_ALLOC((size_t)n_recv * sizeof(ptrdiff_t));
    large_idx_t *uniq_k1     = (large_idx_t *)SMESH_ALLOC((size_t)n_recv * 3 * sizeof(large_idx_t));
    large_idx_t *uniq_aux    = (large_idx_t *)SMESH_ALLOC((size_t)n_recv * sizeof(large_idx_t));
    int         *uniq_owner  = (int *)SMESH_ALLOC((size_t)n_recv * sizeof(int));
    int         *uniq_first  = (int *)SMESH_ALLOC((size_t)n_recv * sizeof(int));
    int         *uniq_shared = (int *)SMESH_CALLOC((size_t)n_recv, sizeof(int));
    ptrdiff_t    n_unique_local = 0;

    for (ptrdiff_t a = 0; a < n_owned_space; ++a) {
        const ptrdiff_t begin           = rowptr[a];
        const ptrdiff_t end             = rowptr[a + 1];
        const ptrdiff_t row_uniq_start  = n_unique_local;
        for (ptrdiff_t t = begin; t < end; ++t) {
            const ptrdiff_t   i   = order[t];
            const large_idx_t *k  = recv_pack + i * k_nfields;
            ptrdiff_t         found = -1;
            for (ptrdiff_t u = row_uniq_start; u < n_unique_local; ++u) {
                if (k[1] == uniq_k1[u * 3] && k[2] == uniq_k1[u * 3 + 1] && k[3] == uniq_k1[u * 3 + 2]) {
                    found = u;
                    break;
                }
            }
            const int         src = (int)k[5];
            const large_idx_t ax  = k[4];
            if (found < 0) {
                const ptrdiff_t u = n_unique_local++;
                uniq_k1[u * 3]     = k[1];
                uniq_k1[u * 3 + 1] = k[2];
                uniq_k1[u * 3 + 2] = k[3];
                uniq_aux[u]        = ax;
                uniq_owner[u]      = src;
                uniq_first[u]      = src;
                uniq_of[i]         = u;
            } else {
                if (ax < uniq_aux[found]) {
                    uniq_aux[found]   = ax;
                    uniq_owner[found] = src;
                }
                if (src != uniq_first[found]) {
                    uniq_shared[found] = 1;
                }
                uniq_of[i] = found;
            }
        }
    }

    long long n_unique_ll = (long long)n_unique_local;
    long long base_ll     = 0;
    SMESH_MPI_CATCH(MPI_Exscan(&n_unique_ll, &base_ll, 1, MPI_LONG_LONG, MPI_SUM, comm));
    if (rank == 0) {
        base_ll = 0;
    }
    long long n_unique_global_ll = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&n_unique_ll, &n_unique_global_ll, 1, MPI_LONG_LONG, MPI_SUM, comm));
    *n_global_unique = (ptrdiff_t)n_unique_global_ll;

    large_idx_t *uniq_gid = (large_idx_t *)SMESH_ALLOC((size_t)n_unique_local * sizeof(large_idx_t));
    for (ptrdiff_t u = 0; u < n_unique_local; ++u) {
        uniq_gid[u] = (large_idx_t)base_ll + (large_idx_t)u;
    }

    memset(send_count, 0, (size_t)size * sizeof(i64));
    memset(send_displs, 0, ((size_t)size + 1) * sizeof(i64));
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const int src = (int)recv_pack[i * k_nfields + 5];
        send_displs[src + 1]++;
    }
    SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count, 1, mpi_type<i64>(), comm));
    recv_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        recv_displs[r + 1] = recv_displs[r] + recv_count[r];
    }
    send_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] += send_displs[r];
    }

    constexpr int k_nreply = 4;
    const ptrdiff_t n_reply_send = (ptrdiff_t)send_displs[size];
    const ptrdiff_t n_reply_recv = (ptrdiff_t)recv_displs[size];
    large_idx_t *reply_send = (large_idx_t *)SMESH_ALLOC((size_t)n_reply_send * k_nreply * sizeof(large_idx_t));
    large_idx_t *reply_recv = (large_idx_t *)SMESH_ALLOC((size_t)n_reply_recv * k_nreply * sizeof(large_idx_t));
    memset(send_count, 0, (size_t)size * sizeof(i64));
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const int         src  = (int)recv_pack[i * k_nfields + 5];
        const i64         slot = send_displs[src] + send_count[src]++;
        const ptrdiff_t   u    = uniq_of[i];
        large_idx_t      *row  = reply_send + (ptrdiff_t)slot * k_nreply;
        row[0] = recv_pack[i * k_nfields + 6];
        row[1] = uniq_gid[u];
        row[2] = (large_idx_t)uniq_owner[u];
        row[3] = (large_idx_t)uniq_shared[u];
    }
    if (all_to_allv_64v(reply_send,
                        send_count,
                        send_displs,
                        reply_recv,
                        recv_count,
                        recv_displs,
                        k_nreply,
                        comm,
                        k_alltoall_chunk) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    for (ptrdiff_t i = 0; i < n_reply_recv; ++i) {
        const large_idx_t *row = reply_recv + i * k_nreply;
        const ptrdiff_t    orig = (ptrdiff_t)row[0];
        out_gid[orig]    = row[1];
        out_owner[orig]  = (int)row[2];
        out_shared[orig] = (int)row[3];
    }

    SMESH_FREE(send_displs);
    SMESH_FREE(send_count);
    SMESH_FREE(recv_displs);
    SMESH_FREE(recv_count);
    SMESH_FREE(send_pack);
    SMESH_FREE(recv_pack);
    SMESH_FREE(rowptr);
    SMESH_FREE(order);
    SMESH_FREE(book);
    SMESH_FREE(uniq_of);
    SMESH_FREE(uniq_k1);
    SMESH_FREE(uniq_aux);
    SMESH_FREE(uniq_owner);
    SMESH_FREE(uniq_first);
    SMESH_FREE(uniq_shared);
    SMESH_FREE(uniq_gid);
    SMESH_FREE(reply_send);
    SMESH_FREE(reply_recv);
    return SMESH_SUCCESS;
}

static int local_unique_by_index(const ptrdiff_t        n_inc,
                                 const large_idx_t     *const keys,
                                 const large_idx_t     *const aux,
                                 const idx_t           *const loc,
                                 const ptrdiff_t        n_index,
                                 ptrdiff_t             *const n_uniq,
                                 large_idx_t          **const out_uniq_keys,
                                 large_idx_t          **const out_uniq_aux,
                                 ptrdiff_t            **const out_inc_to_uniq) {
    if (n_inc == 0) {
        *n_uniq          = 0;
        *out_uniq_keys   = nullptr;
        *out_uniq_aux    = nullptr;
        *out_inc_to_uniq = nullptr;
        return SMESH_SUCCESS;
    }

    ptrdiff_t *rowptr = (ptrdiff_t *)SMESH_CALLOC((size_t)n_index + 1, sizeof(ptrdiff_t));
    for (ptrdiff_t i = 0; i < n_inc; ++i) {
        rowptr[(ptrdiff_t)loc[i] + 1]++;
    }
    for (ptrdiff_t a = 0; a < n_index; ++a) {
        rowptr[a + 1] += rowptr[a];
    }

    ptrdiff_t *order = (ptrdiff_t *)SMESH_ALLOC((size_t)n_inc * sizeof(ptrdiff_t));
    ptrdiff_t *book  = (ptrdiff_t *)SMESH_CALLOC((size_t)n_index, sizeof(ptrdiff_t));
    for (ptrdiff_t i = 0; i < n_inc; ++i) {
        const ptrdiff_t a = (ptrdiff_t)loc[i];
        order[rowptr[a] + book[a]++] = i;
    }

    ptrdiff_t   *inc_to_uniq = (ptrdiff_t *)SMESH_ALLOC((size_t)n_inc * sizeof(ptrdiff_t));
    large_idx_t *uniq_keys   = (large_idx_t *)SMESH_ALLOC((size_t)n_inc * 4 * sizeof(large_idx_t));
    large_idx_t *uniq_aux    = (large_idx_t *)SMESH_ALLOC((size_t)n_inc * sizeof(large_idx_t));
    ptrdiff_t    n_u         = 0;

    for (ptrdiff_t a = 0; a < n_index; ++a) {
        const ptrdiff_t begin          = rowptr[a];
        const ptrdiff_t end            = rowptr[a + 1];
        const ptrdiff_t row_uniq_start = n_u;
        for (ptrdiff_t t = begin; t < end; ++t) {
            const ptrdiff_t    i = order[t];
            const large_idx_t *k = keys + i * 4;
            ptrdiff_t          found = -1;
            for (ptrdiff_t u = row_uniq_start; u < n_u; ++u) {
                if (k[0] == uniq_keys[u * 4] && k[1] == uniq_keys[u * 4 + 1] && k[2] == uniq_keys[u * 4 + 2] &&
                    k[3] == uniq_keys[u * 4 + 3]) {
                    found = u;
                    break;
                }
            }
            if (found < 0) {
                uniq_keys[n_u * 4 + 0] = k[0];
                uniq_keys[n_u * 4 + 1] = k[1];
                uniq_keys[n_u * 4 + 2] = k[2];
                uniq_keys[n_u * 4 + 3] = k[3];
                uniq_aux[n_u]          = aux[i];
                inc_to_uniq[i]         = n_u;
                n_u += 1;
            } else {
                if (aux[i] < uniq_aux[found]) {
                    uniq_aux[found] = aux[i];
                }
                inc_to_uniq[i] = found;
            }
        }
    }

    SMESH_FREE(rowptr);
    SMESH_FREE(order);
    SMESH_FREE(book);
    *n_uniq          = n_u;
    *out_uniq_keys   = uniq_keys;
    *out_uniq_aux    = uniq_aux;
    *out_inc_to_uniq = inc_to_uniq;
    return SMESH_SUCCESS;
}

/// Dense ids: gid = id. Owner/shared via rank_owner(n_space, id) SoA slots.
static int unique_by_id(MPI_Comm           comm,
                        const ptrdiff_t    n_space,
                        const ptrdiff_t    n_req,
                        const large_idx_t *const ids,
                        const large_idx_t *const aux,
                        large_idx_t       *const out_gid,
                        int               *const out_owner,
                        int               *const out_shared) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const ptrdiff_t n_space_pad = id_space_size(n_space, size);
    const ptrdiff_t ids_start   = rank_start(n_space_pad, size, rank);
    const ptrdiff_t n_owned_ids = rank_split(n_space_pad, size, rank);

    i64 *send_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
    i64 *send_count  = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));
    i64 *recv_displs = (i64 *)SMESH_CALLOC((size_t)size + 1, sizeof(i64));
    i64 *recv_count  = (i64 *)SMESH_CALLOC((size_t)size, sizeof(i64));

    for (ptrdiff_t i = 0; i < n_req; ++i) {
        const int p = rank_owner(n_space_pad, (ptrdiff_t)ids[i], size);
        send_displs[p + 1]++;
    }
    SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count, 1, mpi_type<i64>(), comm));
    recv_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        recv_displs[r + 1] = recv_displs[r] + recv_count[r];
    }
    send_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] += send_displs[r];
    }

    const ptrdiff_t n_send = (ptrdiff_t)send_displs[size];
    const ptrdiff_t n_recv = (ptrdiff_t)recv_displs[size];
    constexpr int   k_nfields = 4;
    large_idx_t *send_pack = (large_idx_t *)SMESH_ALLOC((size_t)n_send * k_nfields * sizeof(large_idx_t));
    large_idx_t *recv_pack = (large_idx_t *)SMESH_ALLOC((size_t)n_recv * k_nfields * sizeof(large_idx_t));

    for (ptrdiff_t i = 0; i < n_req; ++i) {
        const int p    = rank_owner(n_space_pad, (ptrdiff_t)ids[i], size);
        const i64 slot = send_displs[p] + send_count[p]++;
        large_idx_t *row = send_pack + (ptrdiff_t)slot * k_nfields;
        row[0] = ids[i];
        row[1] = aux[i];
        row[2] = (large_idx_t)rank;
        row[3] = (large_idx_t)i;
    }
    if (all_to_allv_64v(send_pack,
                        send_count,
                        send_displs,
                        recv_pack,
                        recv_count,
                        recv_displs,
                        k_nfields,
                        comm,
                        k_alltoall_chunk) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }

    int         *seen       = (int *)SMESH_CALLOC((size_t)n_owned_ids, sizeof(int));
    large_idx_t *slot_aux   = (large_idx_t *)SMESH_ALLOC((size_t)n_owned_ids * sizeof(large_idx_t));
    int         *slot_owner = (int *)SMESH_ALLOC((size_t)n_owned_ids * sizeof(int));
    int         *slot_first = (int *)SMESH_ALLOC((size_t)n_owned_ids * sizeof(int));
    int         *slot_shared = (int *)SMESH_CALLOC((size_t)n_owned_ids, sizeof(int));

    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const large_idx_t *row = recv_pack + i * k_nfields;
        const ptrdiff_t    off = (ptrdiff_t)row[0] - ids_start;
        const int          src = (int)row[2];
        if (!seen[off]) {
            seen[off]        = 1;
            slot_aux[off]    = row[1];
            slot_owner[off]  = src;
            slot_first[off]  = src;
            slot_shared[off] = 0;
        } else {
            if (row[1] < slot_aux[off]) {
                slot_aux[off]   = row[1];
                slot_owner[off] = src;
            }
            if (src != slot_first[off]) {
                slot_shared[off] = 1;
            }
        }
    }

    memset(send_count, 0, (size_t)size * sizeof(i64));
    memset(send_displs, 0, ((size_t)size + 1) * sizeof(i64));
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const int src = (int)recv_pack[i * k_nfields + 2];
        send_displs[src + 1]++;
    }
    SMESH_MPI_CATCH(MPI_Alltoall(&send_displs[1], 1, mpi_type<i64>(), recv_count, 1, mpi_type<i64>(), comm));
    recv_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        recv_displs[r + 1] = recv_displs[r] + recv_count[r];
    }
    send_displs[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_displs[r + 1] += send_displs[r];
    }

    constexpr int k_nreply = 3;
    const ptrdiff_t n_reply_send = (ptrdiff_t)send_displs[size];
    const ptrdiff_t n_reply_recv = (ptrdiff_t)recv_displs[size];
    large_idx_t *reply_send = (large_idx_t *)SMESH_ALLOC((size_t)n_reply_send * k_nreply * sizeof(large_idx_t));
    large_idx_t *reply_recv = (large_idx_t *)SMESH_ALLOC((size_t)n_reply_recv * k_nreply * sizeof(large_idx_t));
    memset(send_count, 0, (size_t)size * sizeof(i64));
    for (ptrdiff_t i = 0; i < n_recv; ++i) {
        const large_idx_t *in  = recv_pack + i * k_nfields;
        const int          src = (int)in[2];
        const ptrdiff_t    off = (ptrdiff_t)in[0] - ids_start;
        const i64          slot = send_displs[src] + send_count[src]++;
        large_idx_t       *row  = reply_send + (ptrdiff_t)slot * k_nreply;
        row[0] = in[3];
        row[1] = (large_idx_t)slot_owner[off];
        row[2] = (large_idx_t)slot_shared[off];
    }
    if (all_to_allv_64v(reply_send,
                        send_count,
                        send_displs,
                        reply_recv,
                        recv_count,
                        recv_displs,
                        k_nreply,
                        comm,
                        k_alltoall_chunk) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    for (ptrdiff_t i = 0; i < n_req; ++i) {
        out_gid[i] = ids[i];
    }
    for (ptrdiff_t i = 0; i < n_reply_recv; ++i) {
        const large_idx_t *row  = reply_recv + i * k_nreply;
        const ptrdiff_t    orig = (ptrdiff_t)row[0];
        out_owner[orig]  = (int)row[1];
        out_shared[orig] = (int)row[2];
    }

    SMESH_FREE(send_displs);
    SMESH_FREE(send_count);
    SMESH_FREE(recv_displs);
    SMESH_FREE(recv_count);
    SMESH_FREE(send_pack);
    SMESH_FREE(recv_pack);
    SMESH_FREE(seen);
    SMESH_FREE(slot_aux);
    SMESH_FREE(slot_owner);
    SMESH_FREE(slot_first);
    SMESH_FREE(slot_shared);
    SMESH_FREE(reply_send);
    SMESH_FREE(reply_recv);
    return SMESH_SUCCESS;
}

static void count_entity_nodes(const ptrdiff_t n_uniq,
                               const int       n_per,
                               const int      *const owner,
                               const int      *const shared,
                               const int      *const uo,
                               const int      *const ua,
                               const int       rank,
                               ptrdiff_t       n_bkt[4]) {
    for (ptrdiff_t u = 0; u < n_uniq; ++u) {
        if (!uo[u] && !ua[u]) {
            continue;
        }
        n_bkt[node_bucket(rank, owner[u], shared[u], uo[u], ua[u])] += n_per;
    }
}

static void pack_entity_nodes(const ptrdiff_t    n_uniq,
                              const int          n_per,
                              const large_idx_t  base,
                              const large_idx_t *const gid,
                              const int         *const owner,
                              const int         *const shared,
                              const int         *const uo,
                              const int         *const ua,
                              const int          rank,
                              ptrdiff_t          cur[4],
                              large_idx_t       *const nmap,
                              int               *const nown,
                              idx_t             *const ss_local) {
    for (ptrdiff_t u = 0; u < n_uniq; ++u) {
        if (!uo[u] && !ua[u]) {
            continue;
        }
        const int b = node_bucket(rank, owner[u], shared[u], uo[u], ua[u]);
        for (int t = 0; t < n_per; ++t) {
            const ptrdiff_t w = cur[b]++;
            nmap[w]           = base + gid[u] * (large_idx_t)n_per + (large_idx_t)t;
            nown[w]           = owner[u];
            ss_local[u * n_per + t] = (idx_t)w;
        }
    }
}


static large_idx_t element_gid(const Mesh::Block &block, const large_idx_t concat0, const ptrdiff_t e) {
    const ptrdiff_t   n_owned = block.n_elements_owned();
    const large_idx_t local   = (e < n_owned) ? block.element_mapping()->data()[e]
                                              : block.aura_element_mapping()->data()[e - n_owned];
    return concat0 + local;
}

static void hex_fill_lattice(const int L, int **coords) {
    const int nxe = sshex8_nxe(L);
    for (int d = 0; d < 3; ++d) {
        coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(int));
    }
    for (int zi = 0; zi <= L; ++zi) {
        for (int yi = 0; yi <= L; ++yi) {
            for (int xi = 0; xi <= L; ++xi) {
                const int lidx    = sshex8_lidx(L, xi, yi, zi);
                coords[0][lidx] = xi;
                coords[1][lidx] = yi;
                coords[2][lidx] = zi;
            }
        }
    }
}

static void tet_fill_lattice(const int L, int **coords) {
    const int nxe = sstet4_nxe(L);
    for (int d = 0; d < 3; ++d) {
        coords[d] = (int *)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(int));
    }
    for (int z = 0; z <= L; ++z) {
        for (int y = 0; y <= L - z; ++y) {
            for (int x = 0; x <= L - z - y; ++x) {
                const int lidx    = sstet4_lidx(L, x, y, z);
                coords[0][lidx] = x;
                coords[1][lidx] = y;
                coords[2][lidx] = z;
            }
        }
    }
}

static void write_hex_edge(const int          L,
                           const int         *lagr_to_proteus,
                           int              **coords,
                           const int          d1,
                           const int          d2,
                           const idx_t         edge_start,
                           const ptrdiff_t    e,
                           idx_t **const      ss) {
    const int lid1 = lagr_to_proteus[d1];
    const int lid2 = lagr_to_proteus[d2];
    int       start[3], len[3], dir[3];
    for (int d = 0; d < 3; ++d) {
        start[d] = coords[d][lid1];
        dir[d]   = 1;
        len[d]   = 1;
        int x    = coords[d][lid2] - coords[d][lid1];
        if (x > 0) {
            x -= 1;
            len[d]   = x;
            start[d] = 1;
        } else if (x < 0) {
            x += 1;
            len[d]   = x;
            dir[d]   = -1;
            start[d] = L - 1;
        }
    }
    int en = 0;
    for (int zi = 0; zi != len[2]; zi += dir[2]) {
        for (int yi = 0; yi != len[1]; yi += dir[1]) {
            for (int xi = 0; xi != len[0]; xi += dir[0]) {
                const int lidx = sshex8_lidx(L, start[0] + xi, start[1] + yi, start[2] + zi);
                ss[lidx][e]    = edge_start + en;
                en += 1;
            }
        }
    }
}

static void write_hex_face(const int               L,
                           const large_idx_t       corners[8],
                           const LocalSideTable   &lst,
                           const int              *lagr_to_proteus,
                           int                   **coords,
                           const idx_t              face_offset,
                           const ptrdiff_t         e,
                           const int               f,
                           idx_t **const           ss) {
    int         argmin = 0;
    large_idx_t valmin = corners[lst(f, 0)];
    for (int i = 1; i < 4; ++i) {
        const large_idx_t temp = corners[lst(f, i)];
        if (temp < valmin) {
            argmin = i;
            valmin = temp;
        }
    }
    int lst_o = argmin;
    int lst_u = (lst_o + 1) % 4;
    int lst_v = (lst_o + 3) % 4;
    if (corners[lst(f, lst_u)] > corners[lst(f, lst_v)]) {
        const int tmp = lst_v;
        lst_v         = lst_u;
        lst_u         = tmp;
    }
    const int lidx_o = lagr_to_proteus[lst(f, lst_o)];
    const int lidx_u = lagr_to_proteus[lst(f, lst_u)];
    const int lidx_v = lagr_to_proteus[lst(f, lst_v)];

    int o_start[3];
    int u_len[3], u_dir[3];
    int v_len[3], v_dir[3];
    for (int d = 0; d < 3; ++d) {
        o_start[d] = coords[d][lidx_o];
    }
    for (int d = 0; d < 3; ++d) {
        int x    = coords[d][lidx_u] - coords[d][lidx_o];
        u_dir[d] = 1;
        u_len[d] = 1;
        if (x > 0) {
            x -= 1;
            u_len[d]   = x;
            o_start[d] = 1;
        } else if (x < 0) {
            x += 1;
            u_len[d]   = x;
            u_dir[d]   = -1;
            o_start[d] = L - 1;
        }
    }
    for (int d = 0; d < 3; ++d) {
        int x    = coords[d][lidx_v] - coords[d][lidx_o];
        v_dir[d] = 1;
        v_len[d] = 1;
        if (x > 0) {
            x -= 1;
            v_len[d]   = x;
            o_start[d] = 1;
        } else if (x < 0) {
            x += 1;
            v_len[d]   = x;
            v_dir[d]   = -1;
            o_start[d] = L - 1;
        }
    }
    int local_offset = 0;
    for (int vzi = 0; vzi != v_len[2]; vzi += v_dir[2]) {
        for (int vyi = 0; vyi != v_len[1]; vyi += v_dir[1]) {
            for (int vxi = 0; vxi != v_len[0]; vxi += v_dir[0]) {
                for (int uzi = 0; uzi != u_len[2]; uzi += u_dir[2]) {
                    for (int uyi = 0; uyi != u_len[1]; uyi += u_dir[1]) {
                        for (int uxi = 0; uxi != u_len[0]; uxi += u_dir[0]) {
                            const int pidx = sshex8_lidx(L, uxi + vxi + o_start[0], uyi + vyi + o_start[1], uzi + vzi + o_start[2]);
                            ss[pidx][e]    = face_offset + local_offset++;
                        }
                    }
                }
            }
        }
    }
}

static void write_tet_edge(const int          L,
                           const int         *lagr_to_proteus,
                           int              **coords,
                           const int          d1,
                           const int          d2,
                           const idx_t         edge_start,
                           const ptrdiff_t    e,
                           idx_t **const      ss) {
    const int lid1 = lagr_to_proteus[d1];
    const int lid2 = lagr_to_proteus[d2];
    int       P1[3], P2[3];
    for (int d = 0; d < 3; ++d) {
        P1[d] = coords[d][lid1];
        P2[d] = coords[d][lid2];
    }
    for (int t = 1; t < L; ++t) {
        const int xi   = (P1[0] * (L - t) + P2[0] * t) / L;
        const int yi   = (P1[1] * (L - t) + P2[1] * t) / L;
        const int zi   = (P1[2] * (L - t) + P2[2] * t) / L;
        const int lidx = sstet4_lidx(L, xi, yi, zi);
        ss[lidx][e]    = edge_start + (t - 1);
    }
}

static void write_tet_face(const int               L,
                           const large_idx_t       corners[4],
                           const LocalSideTable   &lst,
                           const int              *lagr_to_proteus,
                           int                   **coords,
                           const idx_t              face_offset,
                           const ptrdiff_t         e,
                           const int               f,
                           idx_t **const           ss) {
    int         argmin = 0;
    large_idx_t valmin = corners[lst(f, 0)];
    for (int i = 1; i < 3; ++i) {
        const large_idx_t temp = corners[lst(f, i)];
        if (temp < valmin) {
            argmin = i;
            valmin = temp;
        }
    }
    int lst_o = argmin;
    int lst_u = (argmin + 1) % 3;
    int lst_v = (argmin + 2) % 3;
    if (corners[lst(f, lst_u)] > corners[lst(f, lst_v)]) {
        const int tmp = lst_v;
        lst_v         = lst_u;
        lst_u         = tmp;
    }
    const int lidx_o = lagr_to_proteus[lst(f, lst_o)];
    const int lidx_u = lagr_to_proteus[lst(f, lst_u)];
    const int lidx_v = lagr_to_proteus[lst(f, lst_v)];
    int       Po[3], Pu[3], Pv[3];
    for (int d = 0; d < 3; ++d) {
        Po[d] = coords[d][lidx_o];
        Pu[d] = coords[d][lidx_u];
        Pv[d] = coords[d][lidx_v];
    }
    int local_offset = 0;
    for (int t = 1; t <= L - 2; ++t) {
        for (int s = 1; s <= L - 1 - t; ++s) {
            const int w    = L - s - t;
            const int xi   = (Po[0] * w + Pu[0] * s + Pv[0] * t) / L;
            const int yi   = (Po[1] * w + Pu[1] * s + Pv[1] * t) / L;
            const int zi   = (Po[2] * w + Pu[2] * s + Pv[2] * t) / L;
            const int pidx = sstet4_lidx(L, xi, yi, zi);
            ss[pidx][e]    = face_offset + local_offset++;
        }
    }
}

}  // namespace

std::shared_ptr<Mesh> to_semistructured_distributed(const int                    level,
                                                    const std::shared_ptr<Mesh> &mesh,
                                                    const bool                   use_GLL) {
    SMESH_TRACE_SCOPE("to_semistructured_distributed");

    enum ElemType family = INVALID;
    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        const enum ElemType f = ss_source_family(mesh->element_type(static_cast<block_idx_t>(b)));
        if (family == INVALID) {
            family = f;
        } else if (f != family) {
            fprintf(stderr, "to_semistructured: mixed-family semistructured conversion is not implemented\n");
            return nullptr;
        }
    }
    if (family != HEX8 && family != TET4) {
        fprintf(stderr,
                "to_semistructured: SS family %s is not implemented\n",
                family == INVALID ? "INVALID" : type_to_string(family));
        return nullptr;
    }
    if (family == TET4 && use_GLL) {
        fprintf(stderr, "to_semistructured: GLL nodes are not implemented for TET SS\n");
        return nullptr;
    }
    if (family == HEX8 && level < 2) {
        fprintf(stderr, "to_semistructured: HEX SS requires level >= 2\n");
        return nullptr;
    }
    if (family == TET4 && level < 1) {
        fprintf(stderr, "to_semistructured: TET SS requires level >= 1\n");
        return nullptr;
    }

    auto              dist = mesh->distributed();
    auto              comm = mesh->comm();
    const int         rank = comm->rank();
    const int         size = comm->size();
    const ptrdiff_t   n_coarse_global = dist->n_nodes_global();
    const ptrdiff_t   n_elem_global   = dist->n_elements_global();
    const large_idx_t *coarse_nmap    = dist->node_mapping()->data();
    const int         *coarse_owner   = dist->node_owner()->data();
    const ptrdiff_t    n_coarse_local = dist->n_nodes_local();
    const ptrdiff_t    n_coarse_owned = dist->n_nodes_owned();
    const ptrdiff_t    n_coarse_ons   = dist->n_nodes_owned_not_shared();

    const int nxe     = (family == HEX8) ? sshex8_nxe(level) : sstet4_nxe(level);
    const int n_macro = (family == HEX8) ? 8 : 4;
    const int nsides  = (family == HEX8) ? 6 : 4;
    const int nnxs    = (family == HEX8) ? 4 : 3;
    const int nxedge  = (family == HEX8) ? (level - 1) : sstet4_nxedge(level);
    const int nxf     = (family == HEX8) ? ((level - 1) * (level - 1)) : sstet4_nxface(level);
    const int nxvol   = (family == HEX8) ? ((level - 1) * (level - 1) * (level - 1)) : sstet4_nxvol(level);

    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        if (mesh->block(b)->n_nodes_per_element() != n_macro) {
            fprintf(stderr,
                    "to_semistructured: block '%s' does not have %d nodes per element\n",
                    mesh->block(b)->name().c_str(),
                    n_macro);
            return nullptr;
        }
    }

    const ptrdiff_t n_blocks = (ptrdiff_t)mesh->n_blocks();
    ptrdiff_t *n_owned_b  = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    ptrdiff_t *n_global_b = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        n_owned_b[b] = mesh->block((size_t)b)->n_elements_owned();
    }
    SMESH_MPI_CATCH(MPI_Allreduce(n_owned_b,
                                  n_global_b,
                                  (int)n_blocks,
                                  mpi_type<ptrdiff_t>(),
                                  MPI_SUM,
                                  comm->get()));
    large_idx_t *concat0 = (large_idx_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(large_idx_t));
    {
        large_idx_t acc = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            concat0[b] = acc;
            acc += (large_idx_t)n_global_b[b];
        }
        SMESH_ASSERT((ptrdiff_t)acc == n_elem_global);
    }

    ptrdiff_t             *n_e        = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    const idx_t *const   **coarse_soa = (const idx_t *const **)SMESH_ALLOC((size_t)n_blocks * sizeof(const idx_t *const *));
    ptrdiff_t              n_e_tot    = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto block    = mesh->block((size_t)b);
        n_e[b]        = block->n_elements();
        coarse_soa[b] = block->elements()->data();
        n_e_tot += n_e[b];
    }

    static const int hex_lagr_conn[8][3] = {{1, 3, 4}, {0, 2, 5}, {1, 3, 6}, {0, 2, 7}, {0, 5, 7}, {1, 4, 6}, {2, 5, 7}, {3, 4, 6}};
    static const int tet_lagr_conn[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

    LocalSideTable lst;
    lst.fill(family == HEX8 ? HEX8 : TET4);

    int hex_corners[8];
    int tet_corners[4];
    if (family == HEX8) {
        hex_corners[0] = sshex8_lidx(level, 0, 0, 0);
        hex_corners[1] = sshex8_lidx(level, level, 0, 0);
        hex_corners[2] = sshex8_lidx(level, level, level, 0);
        hex_corners[3] = sshex8_lidx(level, 0, level, 0);
        hex_corners[4] = sshex8_lidx(level, 0, 0, level);
        hex_corners[5] = sshex8_lidx(level, level, 0, level);
        hex_corners[6] = sshex8_lidx(level, level, level, level);
        hex_corners[7] = sshex8_lidx(level, 0, level, level);
    } else {
        tet_corners[0] = sstet4_lidx(level, 0, 0, 0);
        tet_corners[1] = sstet4_lidx(level, level, 0, 0);
        tet_corners[2] = sstet4_lidx(level, 0, level, 0);
        tet_corners[3] = sstet4_lidx(level, 0, 0, level);
    }
    const int *lagr = (family == HEX8) ? hex_corners : tet_corners;

    int *c_uo = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));
    int *c_ua = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));

    ptrdiff_t    n_edge_inc = 0;
    ptrdiff_t    n_face_inc = 0;
    if (nxedge > 0 || nxf > 0) {
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
                large_idx_t gc[8];
                for (int d = 0; d < n_macro; ++d) {
                    gc[d] = coarse_nmap[coarse_soa[b][d][e]];
                }
                if (nxedge > 0) {
                    for (int d1 = 0; d1 < n_macro; ++d1) {
                        const int *conn = (family == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                        for (int k = 0; k < 3; ++k) {
                            if (gc[d1] <= gc[conn[k]]) {
                                n_edge_inc++;
                            }
                        }
                    }
                }
                n_face_inc += (nxf > 0) ? nsides : 0;
            }
        }
    }

    large_idx_t *edge_keys = nullptr;
    large_idx_t *edge_aux  = nullptr;
    idx_t       *edge_loc  = nullptr;
    if (nxedge > 0 && n_edge_inc > 0) {
        edge_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_edge_inc * 4 * sizeof(large_idx_t));
        edge_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_edge_inc * sizeof(large_idx_t));
        edge_loc  = (idx_t *)SMESH_ALLOC((size_t)n_edge_inc * sizeof(idx_t));
    }
    large_idx_t *face_keys = nullptr;
    large_idx_t *face_aux  = nullptr;
    idx_t       *face_loc  = nullptr;
    if (nxf > 0 && n_face_inc > 0) {
        face_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_face_inc * 4 * sizeof(large_idx_t));
        face_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_face_inc * sizeof(large_idx_t));
        face_loc  = (idx_t *)SMESH_ALLOC((size_t)n_face_inc * sizeof(idx_t));
    }
    large_idx_t *vol_ids = nullptr;
    large_idx_t *vol_aux = nullptr;
    if (nxvol > 0 && n_e_tot > 0) {
        vol_ids = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tot * sizeof(large_idx_t));
        vol_aux = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tot * sizeof(large_idx_t));
    }

    ptrdiff_t ie = 0;
    ptrdiff_t iff = 0;
    ptrdiff_t iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block   = mesh->block((size_t)b);
        const ptrdiff_t n_owned = block->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[8];
            idx_t       lc[8];
            for (int d = 0; d < n_macro; ++d) {
                lc[d] = coarse_soa[b][d][e];
                gc[d] = coarse_nmap[lc[d]];
                if (from_owned) {
                    c_uo[lc[d]] = 1;
                } else {
                    c_ua[lc[d]] = 1;
                }
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    const int *conn = (family == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                    for (int k = 0; k < 3; ++k) {
                        const int d2 = conn[k];
                        if (gc[d1] > gc[d2]) {
                            continue;
                        }
                        edge_keys[ie * 4 + 0] = gc[d1];
                        edge_keys[ie * 4 + 1] = gc[d2];
                        edge_keys[ie * 4 + 2] = k_key_pad;
                        edge_keys[ie * 4 + 3] = k_key_pad;
                        edge_aux[ie]          = (large_idx_t)rank;
                        edge_loc[ie]          = lc[d1];
                        ie++;
                    }
                }
            }
            if (nxf > 0) {
                for (int f = 0; f < nsides; ++f) {
                    large_idx_t fk[4] = {k_key_pad, k_key_pad, k_key_pad, k_key_pad};
                    idx_t       loc_min = lc[lst(f, 0)];
                    large_idx_t gmin    = gc[lst(f, 0)];
                    for (int i = 0; i < nnxs; ++i) {
                        fk[i] = gc[lst(f, i)];
                        if (gc[lst(f, i)] < gmin) {
                            gmin    = gc[lst(f, i)];
                            loc_min = lc[lst(f, i)];
                        }
                    }
                    sort4(fk, nnxs);
                    face_keys[iff * 4 + 0] = fk[0];
                    face_keys[iff * 4 + 1] = fk[1];
                    face_keys[iff * 4 + 2] = fk[2];
                    face_keys[iff * 4 + 3] = fk[3];
                    face_aux[iff]          = element_gid(*block, concat0[b], e);
                    face_loc[iff]          = loc_min;
                    iff++;
                }
            }
            if (nxvol > 0) {
                vol_ids[iv] = element_gid(*block, concat0[b], e);
                vol_aux[iv] = from_owned ? 0 : 1;
                iv++;
            }
        }
    }
    SMESH_ASSERT(ie == n_edge_inc);
    SMESH_ASSERT(iff == n_face_inc);
    SMESH_ASSERT(iv == ((nxvol > 0) ? n_e_tot : 0));

    ptrdiff_t    n_edge_uniq     = 0;
    large_idx_t *edge_uniq_keys  = nullptr;
    large_idx_t *edge_uniq_aux   = nullptr;
    ptrdiff_t   *edge_inc_to_uniq = nullptr;
    large_idx_t *edge_gid        = nullptr;
    int         *edge_owner      = nullptr;
    int         *edge_shared     = nullptr;
    ptrdiff_t    n_edges_global  = 0;
    if (nxedge > 0) {
        if (local_unique_by_index(n_edge_inc,
                                  edge_keys,
                                  edge_aux,
                                  edge_loc,
                                  n_coarse_local,
                                  &n_edge_uniq,
                                  &edge_uniq_keys,
                                  &edge_uniq_aux,
                                  &edge_inc_to_uniq) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(edge_keys);
        SMESH_FREE(edge_aux);
        SMESH_FREE(edge_loc);
        edge_keys = nullptr;
        edge_gid    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1) * sizeof(large_idx_t));
        edge_owner  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1) * sizeof(int));
        edge_shared = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1) * sizeof(int));
        if (unique_tuples(comm->get(),
                          n_coarse_global,
                          n_edge_uniq,
                          edge_uniq_keys,
                          edge_uniq_aux,
                          edge_gid,
                          edge_owner,
                          edge_shared,
                          &n_edges_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(edge_uniq_keys);
        SMESH_FREE(edge_uniq_aux);
        edge_uniq_keys = nullptr;
        edge_uniq_aux  = nullptr;
    }

    ptrdiff_t    n_face_uniq      = 0;
    large_idx_t *face_uniq_keys   = nullptr;
    large_idx_t *face_uniq_aux    = nullptr;
    ptrdiff_t   *face_inc_to_uniq = nullptr;
    large_idx_t *face_gid         = nullptr;
    int         *face_owner       = nullptr;
    int         *face_shared      = nullptr;
    ptrdiff_t    n_faces_global   = 0;
    if (nxf > 0) {
        if (local_unique_by_index(n_face_inc,
                                  face_keys,
                                  face_aux,
                                  face_loc,
                                  n_coarse_local,
                                  &n_face_uniq,
                                  &face_uniq_keys,
                                  &face_uniq_aux,
                                  &face_inc_to_uniq) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(face_keys);
        SMESH_FREE(face_aux);
        SMESH_FREE(face_loc);
        face_keys = nullptr;
        face_gid    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_face_uniq, 1) * sizeof(large_idx_t));
        face_owner  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_face_uniq, 1) * sizeof(int));
        face_shared = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_face_uniq, 1) * sizeof(int));
        if (unique_tuples(comm->get(),
                          n_coarse_global,
                          n_face_uniq,
                          face_uniq_keys,
                          face_uniq_aux,
                          face_gid,
                          face_owner,
                          face_shared,
                          &n_faces_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(face_uniq_keys);
        SMESH_FREE(face_uniq_aux);
        face_uniq_keys = nullptr;
        face_uniq_aux  = nullptr;
    }

    large_idx_t *vol_gid    = nullptr;
    int         *vol_owner  = nullptr;
    int         *vol_shared = nullptr;
    if (nxvol > 0) {
        vol_gid    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1) * sizeof(large_idx_t));
        vol_owner  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1) * sizeof(int));
        vol_shared = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1) * sizeof(int));
        if (unique_by_id(comm->get(), n_elem_global, n_e_tot, vol_ids, vol_aux, vol_gid, vol_owner, vol_shared) !=
            SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(vol_ids);
        SMESH_FREE(vol_aux);
        vol_ids = nullptr;
        vol_aux = nullptr;
    }

    const large_idx_t edge_base = (large_idx_t)n_coarse_global;
    const large_idx_t face_base = edge_base + (large_idx_t)n_edges_global * (large_idx_t)nxedge;
    const large_idx_t vol_base  = face_base + (large_idx_t)n_faces_global * (large_idx_t)nxf;
    const ptrdiff_t   n_ss_global =
            (ptrdiff_t)vol_base + n_elem_global * (ptrdiff_t)nxvol;
    if (n_ss_global < size) {
        fprintf(stderr, "to_semistructured: SS node count smaller than communicator size\n");
        return nullptr;
    }

    int *edge_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *edge_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *face_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_face_uniq, 1), sizeof(int));
    int *face_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_face_uniq, 1), sizeof(int));
    int *vol_uo  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1), sizeof(int));
    int *vol_ua  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1), sizeof(int));

    ie  = 0;
    iff = 0;
    iv  = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const ptrdiff_t n_owned = mesh->block((size_t)b)->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    const int *conn = (family == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                    for (int k = 0; k < 3; ++k) {
                        if (gc[d1] > gc[conn[k]]) {
                            continue;
                        }
                        const ptrdiff_t u = edge_inc_to_uniq[ie++];
                        if (from_owned) {
                            edge_uo[u] = 1;
                        } else {
                            edge_ua[u] = 1;
                        }
                    }
                }
            }
            if (nxf > 0) {
                for (int f = 0; f < nsides; ++f) {
                    const ptrdiff_t u = face_inc_to_uniq[iff++];
                    if (from_owned) {
                        face_uo[u] = 1;
                    } else {
                        face_ua[u] = 1;
                    }
                }
            }
            if (nxvol > 0) {
                if (from_owned) {
                    vol_uo[iv] = 1;
                } else {
                    vol_ua[iv] = 1;
                }
                iv++;
            }
        }
    }

    ptrdiff_t n_bkt[4] = {0, 0, 0, 0};
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        if (!c_uo[i] && !c_ua[i]) {
            continue;
        }
        const int sh = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
        n_bkt[node_bucket(rank, coarse_owner[i], sh, c_uo[i], c_ua[i])]++;
    }
    if (nxedge > 0) {
        count_entity_nodes(n_edge_uniq, nxedge, edge_owner, edge_shared, edge_uo, edge_ua, rank, n_bkt);
    }
    if (nxf > 0) {
        count_entity_nodes(n_face_uniq, nxf, face_owner, face_shared, face_uo, face_ua, rank, n_bkt);
    }
    if (nxvol > 0) {
        count_entity_nodes(n_e_tot, nxvol, vol_owner, vol_shared, vol_uo, vol_ua, rank, n_bkt);
    }

    ptrdiff_t off[5];
    off[0] = 0;
    for (int k = 0; k < 4; ++k) {
        off[k + 1] = off[k] + n_bkt[k];
    }
    const ptrdiff_t n_owned  = off[2];
    const ptrdiff_t n_shared = n_bkt[1];
    const ptrdiff_t n_ghosts = n_bkt[2];
    const ptrdiff_t n_aura   = n_bkt[3];
    const ptrdiff_t n_local  = off[4];
    if (n_local > (ptrdiff_t)std::numeric_limits<idx_t>::max()) {
        fprintf(stderr, "to_semistructured: local SS node count exceeds idx_t\n");
        return nullptr;
    }

    auto        node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto        node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap        = node_mapping->data();
    int         *nown        = node_owner->data();
    auto         points      = create_host_buffer<geom_t>((size_t)mesh->spatial_dimension(), (size_t)n_local);
    auto         coarse_p    = mesh->points()->data();
    auto         p           = points->data();
    const int    sdim        = mesh->spatial_dimension();

    idx_t *corner_ss = (idx_t *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(idx_t));
    idx_t *edge_ss   = (nxedge > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_edge_uniq * (size_t)nxedge * sizeof(idx_t)) : nullptr;
    idx_t *face_ss   = (nxf > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_face_uniq * (size_t)nxf * sizeof(idx_t)) : nullptr;
    idx_t *vol_ss    = (nxvol > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_e_tot * (size_t)nxvol * sizeof(idx_t)) : nullptr;

    ptrdiff_t cur[4] = {off[0], off[1], off[2], off[3]};
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        if (!c_uo[i] && !c_ua[i]) {
            continue;
        }
        const int       sh = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
        const int       b  = node_bucket(rank, coarse_owner[i], sh, c_uo[i], c_ua[i]);
        const ptrdiff_t w  = cur[b]++;
        nmap[w]            = coarse_nmap[i];
        nown[w]            = coarse_owner[i];
        corner_ss[i]       = (idx_t)w;
        for (int d = 0; d < sdim; ++d) {
            p[d][w] = coarse_p[d][i];
        }
    }
    if (nxedge > 0) {
        pack_entity_nodes(n_edge_uniq, nxedge, edge_base, edge_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    }
    if (nxf > 0) {
        pack_entity_nodes(n_face_uniq, nxf, face_base, face_gid, face_owner, face_shared, face_uo, face_ua, rank, cur, nmap, nown, face_ss);
    }
    if (nxvol > 0) {
        pack_entity_nodes(n_e_tot, nxvol, vol_base, vol_gid, vol_owner, vol_shared, vol_uo, vol_ua, rank, cur, nmap, nown, vol_ss);
    }

    int *coords[3] = {nullptr, nullptr, nullptr};
    if (family == HEX8) {
        hex_fill_lattice(level, coords);
    } else {
        tet_fill_lattice(level, coords);
    }

    const enum ElemType ss_type = semistructured_type(family, level);
    std::vector<std::shared_ptr<Mesh::Block>> ss_blocks((size_t)n_blocks);
    ie  = 0;
    iff = 0;
    iv  = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            coarse_block = mesh->block((size_t)b);
        auto     ss_elems = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_e[b]);
        idx_t **out       = ss_elems->data();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                const idx_t local = coarse_soa[b][d][e];
                gc[d]             = coarse_nmap[local];
                out[lagr[d]][e]   = corner_ss[local];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    const int *conn = (family == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                    for (int k = 0; k < 3; ++k) {
                        const int d2 = conn[k];
                        if (gc[d1] > gc[d2]) {
                            continue;
                        }
                        const idx_t estart = edge_ss[edge_inc_to_uniq[ie++] * nxedge];
                        if (family == HEX8) {
                            write_hex_edge(level, hex_corners, coords, d1, d2, estart, e, out);
                        } else {
                            write_tet_edge(level, tet_corners, coords, d1, d2, estart, e, out);
                        }
                    }
                }
            }
            if (nxf > 0) {
                for (int f = 0; f < nsides; ++f) {
                    const idx_t foff = face_ss[face_inc_to_uniq[iff++] * nxf];
                    if (family == HEX8) {
                        write_hex_face(level, gc, lst, hex_corners, coords, foff, e, f, out);
                    } else {
                        write_tet_face(level, gc, lst, tet_corners, coords, foff, e, f, out);
                    }
                }
            }
            if (nxvol > 0) {
                const idx_t voff = vol_ss[iv * nxvol];
                if (family == HEX8) {
                    const int Lm1 = level - 1;
                    for (int zi = 1; zi < level; ++zi) {
                        for (int yi = 1; yi < level; ++yi) {
                            for (int xi = 1; xi < level; ++xi) {
                                const int lidx = sshex8_lidx(level, xi, yi, zi);
                                const int en   = (zi - 1) * Lm1 * Lm1 + (yi - 1) * Lm1 + (xi - 1);
                                out[lidx][e]   = voff + (idx_t)en;
                            }
                        }
                    }
                } else {
                    for (int z = 1; z <= level - 3; ++z) {
                        for (int y = 1; y <= level - 2 - z; ++y) {
                            for (int x = 1; x <= level - 1 - z - y; ++x) {
                                const int lidx = sstet4_lidx(level, x, y, z);
                                const int en   = sstet4_lidx(level - 4, x - 1, y - 1, z - 1);
                                out[lidx][e]   = voff + (idx_t)en;
                            }
                        }
                    }
                }
                iv++;
            }
        }
        auto ss_block = std::make_shared<Mesh::Block>();
        ss_block->set_name(coarse_block->name());
        ss_block->set_element_type(ss_type);
        ss_block->set_elements(ss_elems);
        ss_block->set_distributed_elements(coarse_block->n_elements_owned(),
                                           coarse_block->n_elements_shared(),
                                           coarse_block->n_elements_ghosts(),
                                           coarse_block->element_mapping(),
                                           coarse_block->aura_element_mapping());
        ss_blocks[(size_t)b] = ss_block;
    }
    for (int d = 0; d < 3; ++d) {
        SMESH_FREE(coords[d]);
    }

    auto ghosts_and_aura = create_host_buffer<idx_t>((size_t)(n_ghosts + n_aura));
    auto node_offsets    = create_host_buffer<ptrdiff_t>((size_t)size + 1);
    node_ownership_ranges(comm->get(), n_owned, node_offsets->data());
    SMESH_ASSERT(node_offsets->data()[size] == n_ss_global);

    if (n_ghosts + n_aura > 0) {
        const ptrdiff_t n_id = rank_split(n_ss_global, size, rank);
        idx_t *global2owned = (idx_t *)SMESH_CALLOC((size_t)n_id, sizeof(idx_t));
        if (prepare_node_renumbering(comm->get(),
                                     n_ss_global,
                                     node_offsets->data()[rank],
                                     n_owned,
                                     nmap,
                                     global2owned) != SMESH_SUCCESS) {
            return nullptr;
        }
        if (collect_ghost_and_aura_import_indices(comm->get(),
                                                  n_owned,
                                                  n_ghosts,
                                                  n_aura,
                                                  n_ss_global,
                                                  nmap,
                                                  global2owned,
                                                  node_offsets->data(),
                                                  ghosts_and_aura->data()) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(global2owned);
    }

    const double *qx = nullptr;
    if (use_GLL) {
        switch (level) {
            case 2:
                qx = line_GL_q3_x;
                break;
            case 4:
                qx = line_GL_q5_x;
                break;
            case 8:
                qx = line_GL_q9_x;
                break;
            case 16:
                qx = line_GL_q17_x;
                break;
            default:
                SMESH_ERROR("Unsupported GLL order %d!", level);
                return nullptr;
        }
    }

    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        idx_t **el = ss_blocks[(size_t)b]->elements()->data();
        if (family == HEX8) {
            if (use_GLL) {
                sshex8_fill_points_1D_map(level, n_e[b], el, p, qx, p);
            } else {
                sshex8_fill_points(level, n_e[b], el, p, p);
            }
        } else {
            sstet4_fill_points(level, n_e[b], el, p, p);
        }
    }

    SMESH_FREE(n_owned_b);
    SMESH_FREE(n_global_b);
    SMESH_FREE(concat0);
    SMESH_FREE(n_e);
    SMESH_FREE(coarse_soa);
    SMESH_FREE(c_uo);
    SMESH_FREE(c_ua);
    SMESH_FREE(edge_inc_to_uniq);
    SMESH_FREE(edge_gid);
    SMESH_FREE(edge_owner);
    SMESH_FREE(edge_shared);
    SMESH_FREE(face_inc_to_uniq);
    SMESH_FREE(face_gid);
    SMESH_FREE(face_owner);
    SMESH_FREE(face_shared);
    SMESH_FREE(vol_gid);
    SMESH_FREE(vol_owner);
    SMESH_FREE(vol_shared);
    SMESH_FREE(edge_uo);
    SMESH_FREE(edge_ua);
    SMESH_FREE(face_uo);
    SMESH_FREE(face_ua);
    SMESH_FREE(vol_uo);
    SMESH_FREE(vol_ua);
    SMESH_FREE(corner_ss);
    SMESH_FREE(edge_ss);
    SMESH_FREE(face_ss);
    SMESH_FREE(vol_ss);

    auto ret     = std::make_shared<Mesh>(comm, ss_blocks, points);
    auto ss_dist = std::make_shared<Distributed>();
    ss_dist->set_nodes(n_ss_global,
                       n_owned,
                       n_shared,
                       n_ghosts,
                       n_aura,
                       node_mapping,
                       node_owner,
                       node_offsets,
                       ghosts_and_aura);
    ss_dist->set_elements(dist->n_elements_global(),
                          dist->n_elements_owned(),
                          dist->n_elements_shared(),
                          dist->n_elements_ghosts(),
                          dist->element_mapping(),
                          dist->aura_element_mapping());
    ret->set_distributed(ss_dist);
    return ret;
}

#else

std::shared_ptr<Mesh> to_semistructured_distributed(const int /*level*/,
                                                    const std::shared_ptr<Mesh> & /*mesh*/,
                                                    const bool /*use_GLL*/) {
    return nullptr;
}

#endif

}  // namespace smesh

