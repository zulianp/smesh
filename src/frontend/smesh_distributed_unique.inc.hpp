// Extracted A8 unique-numbering / pack helpers.
// Included from smesh_semistructured_distributed.cpp (anonymous namespace)
// and smesh_mesh_transforms_distributed.cpp. Do not edit HEX/TET unique-map loops.

constexpr i64 k_alltoall_chunk = static_cast<i64>(1) << 20;
constexpr large_idx_t k_key_pad = static_cast<large_idx_t>(-1);

/// unique_tuples / unique_by_id pick the min aux incidence as owner. Aura-only
/// hits on a low rank must not win over owned-element hits on another rank.
static large_idx_t owned_pref_rank_aux(const int from_owned, const int rank, const int comm_size) {
    return from_owned ? (large_idx_t)rank : (large_idx_t)rank + (large_idx_t)comm_size;
}

#if !defined(SMESH_DISTRIBUTED_UNIQUE_MINIMAL)
/// WARNING: unused in mesh_transforms_distributed (SS unique-map only).
static large_idx_t owned_pref_eid_aux(const int from_owned, const large_idx_t eid, const ptrdiff_t n_elem_global) {
    return from_owned ? eid : eid + (large_idx_t)n_elem_global;
}

/// WARNING: unused in mesh_transforms_distributed (SS unique-map only).
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
#endif

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
#if !defined(SMESH_DISTRIBUTED_UNIQUE_MINIMAL)
/// WARNING: unused in mesh_transforms_distributed (SS unique-map only).
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
#endif

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
#if !defined(SMESH_DISTRIBUTED_UNIQUE_MINIMAL)
/// WARNING: unused in mesh_transforms_distributed (SS unique-map only).
static void fill_flat_node_gids(const ptrdiff_t          n_uniq,
                                const int                n_per,
                                const large_idx_t        base,
                                const large_idx_t *const gid,
                                large_idx_t             *const out) {
    for (ptrdiff_t u = 0; u < n_uniq; ++u) {
        for (int t = 0; t < n_per; ++t) {
            out[u * (ptrdiff_t)n_per + t] = base + gid[u] * (large_idx_t)n_per + (large_idx_t)t;
        }
    }
}
#endif
static int unique_inc_tuples(MPI_Comm           comm,
                             const ptrdiff_t    n_space,
                             const ptrdiff_t    n_inc,
                             large_idx_t       *keys,
                             large_idx_t       *aux,
                             idx_t             *loc,
                             const ptrdiff_t    n_index,
                             ptrdiff_t         *const n_uniq,
                             ptrdiff_t        **const inc_to_uniq,
                             large_idx_t      **const gid,
                             int              **const owner,
                             int              **const shared,
                             ptrdiff_t         *const n_global) {
    large_idx_t *uniq_keys = nullptr;
    large_idx_t *uniq_aux  = nullptr;
    if (local_unique_by_index(n_inc, keys, aux, loc, n_index, n_uniq, &uniq_keys, &uniq_aux, inc_to_uniq) !=
        SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    SMESH_FREE(keys);
    SMESH_FREE(aux);
    SMESH_FREE(loc);
    *gid    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(*n_uniq, 1) * sizeof(large_idx_t));
    *owner  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(*n_uniq, 1) * sizeof(int));
    *shared = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(*n_uniq, 1) * sizeof(int));
    if (unique_tuples(comm, n_space, *n_uniq, uniq_keys, uniq_aux, *gid, *owner, *shared, n_global) !=
        SMESH_SUCCESS) {
        SMESH_FREE(uniq_keys);
        SMESH_FREE(uniq_aux);
        return SMESH_FAILURE;
    }
    SMESH_FREE(uniq_keys);
    SMESH_FREE(uniq_aux);
    return SMESH_SUCCESS;
}

static large_idx_t *alloc_entity_node_gids(const ptrdiff_t n_uniq, const int n_per) {
    if (n_uniq <= 0 || n_per <= 0) {
        return nullptr;
    }
    return (large_idx_t *)SMESH_ALLOC((size_t)n_uniq * (size_t)n_per * sizeof(large_idx_t));
}

static void pack_entity_nodes(const ptrdiff_t          n_uniq,
                              const int                n_per,
                              const large_idx_t *const node_gid,
                              const int         *const owner,
                              const int         *const shared,
                              const int         *const uo,
                              const int         *const ua,
                              const int                rank,
                              ptrdiff_t                cur[4],
                              large_idx_t             *const nmap,
                              int                     *const nown,
                              idx_t                   *const ss_local) {
    for (ptrdiff_t u = 0; u < n_uniq; ++u) {
        if (!uo[u] && !ua[u]) {
            continue;
        }
        const int b = node_bucket(rank, owner[u], shared[u], uo[u], ua[u]);
        for (int t = 0; t < n_per; ++t) {
            const ptrdiff_t w = cur[b]++;
            nmap[w]           = node_gid[u * (ptrdiff_t)n_per + t];
            nown[w]           = owner[u];
            ss_local[u * (ptrdiff_t)n_per + t] = (idx_t)w;
        }
    }
}


