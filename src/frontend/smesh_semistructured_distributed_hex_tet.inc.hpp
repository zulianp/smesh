std::shared_ptr<Mesh> to_semistructured_distributed_hex_tet(const int                    level,
                                                                   const std::shared_ptr<Mesh> &mesh,
                                                                   const bool                   hierarchical_ordering) {
    SMESH_TRACE_SCOPE("to_semistructured_distributed_hex_tet");

    if (level < 1) {
        fprintf(stderr, "to_semistructured: mixed HEX+TET SS requires level >= 1\n");
        return nullptr;
    }

    auto               dist            = mesh->distributed();
    auto               comm            = mesh->comm();
    const int          rank            = comm->rank();
    const int          size            = comm->size();
    const ptrdiff_t    n_coarse_global = dist->n_nodes_global();
    const ptrdiff_t    n_elem_global   = dist->n_elements_global();
    const large_idx_t *coarse_nmap     = dist->node_mapping()->data();
    const int         *coarse_owner    = dist->node_owner()->data();
    const ptrdiff_t    n_coarse_local  = dist->n_nodes_local();
    const ptrdiff_t    n_coarse_owned  = dist->n_nodes_owned();
    const ptrdiff_t    n_coarse_ons    = dist->n_nodes_owned_not_shared();

    const int nxedge    = level > 1 ? (level - 1) : 0;
    const int nxf_hex   = (level - 1) * (level - 1);
    const int nxf_tet   = sstet4_nxface(level);
    const int nxvol_hex = (level - 1) * (level - 1) * (level - 1);
    const int nxvol_tet = sstet4_nxvol(level);

    const ptrdiff_t n_blocks = (ptrdiff_t)mesh->n_blocks();
    enum ElemType *block_types = (enum ElemType *)SMESH_ALLOC((size_t)n_blocks * sizeof(enum ElemType));
    ptrdiff_t     *n_owned_b   = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    ptrdiff_t     *n_global_b  = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto                block  = mesh->block((size_t)b);
        const enum ElemType family = ss_source_family(mesh->element_type(static_cast<block_idx_t>(b)));
        block_types[b]             = family;
        n_owned_b[b]               = block->n_elements_owned();
        const int n_macro          = (family == HEX8) ? 8 : 4;
        if (block->n_nodes_per_element() != n_macro) {
            fprintf(stderr,
                    "to_semistructured: block '%s' does not have %d nodes per element\n",
                    block->name().c_str(),
                    n_macro);
            return nullptr;
        }
    }
    SMESH_MPI_CATCH(MPI_Allreduce(n_owned_b,
                                  n_global_b,
                                  (int)n_blocks,
                                  mpi_type<ptrdiff_t>(),
                                  MPI_SUM,
                                  comm->get()));

    large_idx_t *hex_concat0 = (large_idx_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(large_idx_t));
    large_idx_t *tet_concat0 = (large_idx_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(large_idx_t));
    large_idx_t  hex_acc     = 0;
    large_idx_t  tet_acc     = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        if (block_types[b] == HEX8) {
            hex_concat0[b] = hex_acc;
            hex_acc += (large_idx_t)n_global_b[b];
            tet_concat0[b] = 0;
        } else {
            tet_concat0[b] = tet_acc;
            tet_acc += (large_idx_t)n_global_b[b];
            hex_concat0[b] = 0;
        }
    }
    const ptrdiff_t n_hex_elem_global = (ptrdiff_t)hex_acc;
    const ptrdiff_t n_tet_elem_global = (ptrdiff_t)tet_acc;
    SMESH_ASSERT(n_hex_elem_global + n_tet_elem_global == n_elem_global);

    ptrdiff_t           *n_e        = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    const idx_t *const **coarse_soa = (const idx_t *const **)SMESH_ALLOC((size_t)n_blocks * sizeof(const idx_t *const *));
    ptrdiff_t            n_e_hex    = 0;
    ptrdiff_t            n_e_tet    = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto block    = mesh->block((size_t)b);
        n_e[b]        = block->n_elements();
        coarse_soa[b] = block->elements()->data();
        if (block_types[b] == HEX8) {
            n_e_hex += n_e[b];
        } else {
            n_e_tet += n_e[b];
        }
    }

    static const int hex_lagr_conn[8][3] = {{1, 3, 4}, {0, 2, 5}, {1, 3, 6}, {0, 2, 7}, {0, 5, 7}, {1, 4, 6}, {2, 5, 7}, {3, 4, 6}};
    static const int tet_lagr_conn[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

    LocalSideTable lst_hex;
    LocalSideTable lst_tet;
    lst_hex.fill(HEX8);
    lst_tet.fill(TET4);

    int hex_corners[8];
    int tet_corners[4];
    hex_corners[0] = sshex8_lidx(level, 0, 0, 0);
    hex_corners[1] = sshex8_lidx(level, level, 0, 0);
    hex_corners[2] = sshex8_lidx(level, level, level, 0);
    hex_corners[3] = sshex8_lidx(level, 0, level, 0);
    hex_corners[4] = sshex8_lidx(level, 0, 0, level);
    hex_corners[5] = sshex8_lidx(level, level, 0, level);
    hex_corners[6] = sshex8_lidx(level, level, level, level);
    hex_corners[7] = sshex8_lidx(level, 0, level, level);
    tet_corners[0] = sstet4_lidx(level, 0, 0, 0);
    tet_corners[1] = sstet4_lidx(level, level, 0, 0);
    tet_corners[2] = sstet4_lidx(level, 0, level, 0);
    tet_corners[3] = sstet4_lidx(level, 0, 0, level);

    int *c_uo = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));
    int *c_ua = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));

    ptrdiff_t n_edge_inc      = 0;
    ptrdiff_t n_hex_face_inc  = 0;
    ptrdiff_t n_tet_face_inc  = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const int n_macro = (block_types[b] == HEX8) ? 8 : 4;
        const int nsides  = (block_types[b] == HEX8) ? 6 : 4;
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    const int *conn = (block_types[b] == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                    for (int k = 0; k < 3; ++k) {
                        if (gc[d1] <= gc[conn[k]]) {
                            n_edge_inc++;
                        }
                    }
                }
            }
            if (block_types[b] == HEX8) {
                n_hex_face_inc += (nxf_hex > 0) ? nsides : 0;
            } else {
                n_tet_face_inc += (nxf_tet > 0) ? nsides : 0;
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
    large_idx_t *hex_face_keys = nullptr;
    large_idx_t *hex_face_aux  = nullptr;
    idx_t       *hex_face_loc  = nullptr;
    if (nxf_hex > 0 && n_hex_face_inc > 0) {
        hex_face_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_hex_face_inc * 4 * sizeof(large_idx_t));
        hex_face_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_hex_face_inc * sizeof(large_idx_t));
        hex_face_loc  = (idx_t *)SMESH_ALLOC((size_t)n_hex_face_inc * sizeof(idx_t));
    }
    large_idx_t *tet_face_keys = nullptr;
    large_idx_t *tet_face_aux  = nullptr;
    idx_t       *tet_face_loc  = nullptr;
    if (nxf_tet > 0 && n_tet_face_inc > 0) {
        tet_face_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_tet_face_inc * 4 * sizeof(large_idx_t));
        tet_face_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_tet_face_inc * sizeof(large_idx_t));
        tet_face_loc  = (idx_t *)SMESH_ALLOC((size_t)n_tet_face_inc * sizeof(idx_t));
    }
    large_idx_t *hex_vol_ids = nullptr;
    large_idx_t *hex_vol_aux = nullptr;
    if (nxvol_hex > 0 && n_e_hex > 0) {
        hex_vol_ids = (large_idx_t *)SMESH_ALLOC((size_t)n_e_hex * sizeof(large_idx_t));
        hex_vol_aux = (large_idx_t *)SMESH_ALLOC((size_t)n_e_hex * sizeof(large_idx_t));
    }
    large_idx_t *tet_vol_ids = nullptr;
    large_idx_t *tet_vol_aux = nullptr;
    if (nxvol_tet > 0 && n_e_tet > 0) {
        tet_vol_ids = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tet * sizeof(large_idx_t));
        tet_vol_aux = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tet * sizeof(large_idx_t));
    }

    ptrdiff_t ie  = 0;
    ptrdiff_t ihf = 0;
    ptrdiff_t itf = 0;
    ptrdiff_t ihv = 0;
    ptrdiff_t itv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block   = mesh->block((size_t)b);
        const ptrdiff_t n_owned = block->n_elements_owned();
        const int       n_macro = (block_types[b] == HEX8) ? 8 : 4;
        const int       nsides  = (block_types[b] == HEX8) ? 6 : 4;
        const int       nnxs    = (block_types[b] == HEX8) ? 4 : 3;
        const LocalSideTable &lst = (block_types[b] == HEX8) ? lst_hex : lst_tet;
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
                    const int *conn = (block_types[b] == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                    for (int k = 0; k < 3; ++k) {
                        const int d2 = conn[k];
                        if (gc[d1] > gc[d2]) {
                            continue;
                        }
                        edge_keys[ie * 4 + 0] = gc[d1];
                        edge_keys[ie * 4 + 1] = gc[d2];
                        edge_keys[ie * 4 + 2] = k_key_pad;
                        edge_keys[ie * 4 + 3] = k_key_pad;
                        edge_aux[ie]          = owned_pref_rank_aux(from_owned, rank, size);
                        edge_loc[ie]          = lc[d1];
                        ie++;
                    }
                }
            }
            if (block_types[b] == HEX8 && nxf_hex > 0) {
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
                    hex_face_keys[ihf * 4 + 0] = fk[0];
                    hex_face_keys[ihf * 4 + 1] = fk[1];
                    hex_face_keys[ihf * 4 + 2] = fk[2];
                    hex_face_keys[ihf * 4 + 3] = fk[3];
                    hex_face_aux[ihf]          = owned_pref_eid_aux(from_owned, element_gid(*block, hex_concat0[b], e), n_elem_global);
                    hex_face_loc[ihf]          = loc_min;
                    ihf++;
                }
            } else if (block_types[b] == TET4 && nxf_tet > 0) {
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
                    tet_face_keys[itf * 4 + 0] = fk[0];
                    tet_face_keys[itf * 4 + 1] = fk[1];
                    tet_face_keys[itf * 4 + 2] = fk[2];
                    tet_face_keys[itf * 4 + 3] = fk[3];
                    tet_face_aux[itf]          = owned_pref_eid_aux(from_owned, element_gid(*block, tet_concat0[b], e), n_elem_global);
                    tet_face_loc[itf]          = loc_min;
                    itf++;
                }
            }
            if (block_types[b] == HEX8 && nxvol_hex > 0) {
                hex_vol_ids[ihv] = element_gid(*block, hex_concat0[b], e);
                hex_vol_aux[ihv] = from_owned ? 0 : 1;
                ihv++;
            } else if (block_types[b] == TET4 && nxvol_tet > 0) {
                tet_vol_ids[itv] = element_gid(*block, tet_concat0[b], e);
                tet_vol_aux[itv] = from_owned ? 0 : 1;
                itv++;
            }
        }
    }
    SMESH_ASSERT(ie == n_edge_inc);
    SMESH_ASSERT(ihf == n_hex_face_inc);
    SMESH_ASSERT(itf == n_tet_face_inc);
    SMESH_ASSERT(ihv == ((nxvol_hex > 0) ? n_e_hex : 0));
    SMESH_ASSERT(itv == ((nxvol_tet > 0) ? n_e_tet : 0));

    ptrdiff_t    n_edge_uniq      = 0;
    ptrdiff_t   *edge_inc_to_uniq = nullptr;
    large_idx_t *edge_gid         = nullptr;
    int         *edge_owner       = nullptr;
    int         *edge_shared      = nullptr;
    ptrdiff_t    n_edges_global   = 0;
    if (nxedge > 0) {
        if (unique_inc_tuples(comm->get(),
                              n_coarse_global,
                              n_edge_inc,
                              edge_keys,
                              edge_aux,
                              edge_loc,
                              n_coarse_local,
                              &n_edge_uniq,
                              &edge_inc_to_uniq,
                              &edge_gid,
                              &edge_owner,
                              &edge_shared,
                              &n_edges_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        edge_keys = nullptr;
        edge_aux  = nullptr;
        edge_loc  = nullptr;
    }

    ptrdiff_t    n_hex_face_uniq      = 0;
    ptrdiff_t   *hex_face_inc_to_uniq = nullptr;
    large_idx_t *hex_face_gid         = nullptr;
    int         *hex_face_owner       = nullptr;
    int         *hex_face_shared      = nullptr;
    ptrdiff_t    n_hex_faces_global   = 0;
    if (nxf_hex > 0) {
        if (unique_inc_tuples(comm->get(),
                              n_coarse_global,
                              n_hex_face_inc,
                              hex_face_keys,
                              hex_face_aux,
                              hex_face_loc,
                              n_coarse_local,
                              &n_hex_face_uniq,
                              &hex_face_inc_to_uniq,
                              &hex_face_gid,
                              &hex_face_owner,
                              &hex_face_shared,
                              &n_hex_faces_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        hex_face_keys = nullptr;
        hex_face_aux  = nullptr;
        hex_face_loc  = nullptr;
    }

    ptrdiff_t    n_tet_face_uniq      = 0;
    ptrdiff_t   *tet_face_inc_to_uniq = nullptr;
    large_idx_t *tet_face_gid         = nullptr;
    int         *tet_face_owner       = nullptr;
    int         *tet_face_shared      = nullptr;
    ptrdiff_t    n_tet_faces_global   = 0;
    if (nxf_tet > 0) {
        if (unique_inc_tuples(comm->get(),
                              n_coarse_global,
                              n_tet_face_inc,
                              tet_face_keys,
                              tet_face_aux,
                              tet_face_loc,
                              n_coarse_local,
                              &n_tet_face_uniq,
                              &tet_face_inc_to_uniq,
                              &tet_face_gid,
                              &tet_face_owner,
                              &tet_face_shared,
                              &n_tet_faces_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        tet_face_keys = nullptr;
        tet_face_aux  = nullptr;
        tet_face_loc  = nullptr;
    }

    large_idx_t *hex_vol_gid    = nullptr;
    int         *hex_vol_owner  = nullptr;
    int         *hex_vol_shared = nullptr;
    if (nxvol_hex > 0) {
        hex_vol_gid    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_hex, 1) * sizeof(large_idx_t));
        hex_vol_owner  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_hex, 1) * sizeof(int));
        hex_vol_shared = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_hex, 1) * sizeof(int));
        if (unique_by_id(comm->get(),
                         n_hex_elem_global,
                         n_e_hex,
                         hex_vol_ids,
                         hex_vol_aux,
                         hex_vol_gid,
                         hex_vol_owner,
                         hex_vol_shared) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(hex_vol_ids);
        SMESH_FREE(hex_vol_aux);
        hex_vol_ids = nullptr;
        hex_vol_aux = nullptr;
    }

    large_idx_t *tet_vol_gid    = nullptr;
    int         *tet_vol_owner  = nullptr;
    int         *tet_vol_shared = nullptr;
    if (nxvol_tet > 0) {
        tet_vol_gid    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_tet, 1) * sizeof(large_idx_t));
        tet_vol_owner  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_tet, 1) * sizeof(int));
        tet_vol_shared = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_tet, 1) * sizeof(int));
        if (unique_by_id(comm->get(),
                         n_tet_elem_global,
                         n_e_tet,
                         tet_vol_ids,
                         tet_vol_aux,
                         tet_vol_gid,
                         tet_vol_owner,
                         tet_vol_shared) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(tet_vol_ids);
        SMESH_FREE(tet_vol_aux);
        tet_vol_ids = nullptr;
        tet_vol_aux = nullptr;
    }

    const large_idx_t edge_base      = (large_idx_t)n_coarse_global;
    const large_idx_t hex_face_base  = edge_base + (large_idx_t)n_edges_global * (large_idx_t)nxedge;
    const large_idx_t tet_face_base  = hex_face_base + (large_idx_t)n_hex_faces_global * (large_idx_t)nxf_hex;
    const large_idx_t hex_vol_base   = tet_face_base + (large_idx_t)n_tet_faces_global * (large_idx_t)nxf_tet;
    const large_idx_t tet_vol_base   = hex_vol_base + (large_idx_t)n_hex_elem_global * (large_idx_t)nxvol_hex;
    const ptrdiff_t   n_ss_global =
            (ptrdiff_t)tet_vol_base + n_tet_elem_global * (ptrdiff_t)nxvol_tet;
    if (n_ss_global < size) {
        fprintf(stderr, "to_semistructured: SS node count smaller than communicator size\n");
        return nullptr;
    }

    const int nlevels = hierarchical_ordering ? sshex8_hierarchical_n_levels(level) : 0;
    int *levels           = nullptr;
    int *edge_layer       = nullptr;
    int *edge_trank       = nullptr;
    int *hex_face_layer   = nullptr;
    int *hex_face_trank   = nullptr;
    int *tet_face_layer   = nullptr;
    int *tet_face_trank   = nullptr;
    int *hex_vol_layer    = nullptr;
    int *hex_vol_trank    = nullptr;
    int *tet_vol_layer    = nullptr;
    int *tet_vol_trank    = nullptr;
    int *n_edge_t         = nullptr;
    int *n_hex_face_t     = nullptr;
    int *n_tet_face_t     = nullptr;
    int *n_hex_vol_t      = nullptr;
    int *n_tet_vol_t      = nullptr;
    large_idx_t *layer_base = nullptr;
    if (hierarchical_ordering) {
        if (nlevels < 1) {
            fprintf(stderr, "to_semistructured: hierarchical mesh levels cannot be formed\n");
            return nullptr;
        }
        levels       = (int *)SMESH_ALLOC((size_t)nlevels * sizeof(int));
        n_edge_t     = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_hex_face_t = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_tet_face_t = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_hex_vol_t  = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_tet_vol_t  = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        layer_base   = (large_idx_t *)SMESH_ALLOC((size_t)(nlevels + 1) * sizeof(large_idx_t));
        sshex8_hierarchical_mesh_levels(level, nlevels, levels);
        if (nxedge > 0) {
            edge_layer = (int *)SMESH_ALLOC((size_t)nxedge * sizeof(int));
            edge_trank = (int *)SMESH_ALLOC((size_t)nxedge * sizeof(int));
            hier_fill_edge_layers(level, nxedge, nlevels, levels, edge_layer);
            hier_slot_ranks(nxedge, edge_layer, nlevels, edge_trank, n_edge_t);
        }
        if (nxf_hex > 0) {
            hex_face_layer = (int *)SMESH_ALLOC((size_t)nxf_hex * sizeof(int));
            hex_face_trank = (int *)SMESH_ALLOC((size_t)nxf_hex * sizeof(int));
            hier_fill_face_layers(HEX8, level, nxf_hex, nlevels, levels, hex_face_layer);
            hier_slot_ranks(nxf_hex, hex_face_layer, nlevels, hex_face_trank, n_hex_face_t);
        }
        if (nxf_tet > 0) {
            tet_face_layer = (int *)SMESH_ALLOC((size_t)nxf_tet * sizeof(int));
            tet_face_trank = (int *)SMESH_ALLOC((size_t)nxf_tet * sizeof(int));
            hier_fill_face_layers(TET4, level, nxf_tet, nlevels, levels, tet_face_layer);
            hier_slot_ranks(nxf_tet, tet_face_layer, nlevels, tet_face_trank, n_tet_face_t);
        }
        if (nxvol_hex > 0) {
            hex_vol_layer = (int *)SMESH_ALLOC((size_t)nxvol_hex * sizeof(int));
            hex_vol_trank = (int *)SMESH_ALLOC((size_t)nxvol_hex * sizeof(int));
            hier_fill_vol_layers(HEX8, level, nlevels, levels, hex_vol_layer);
            hier_slot_ranks(nxvol_hex, hex_vol_layer, nlevels, hex_vol_trank, n_hex_vol_t);
        }
        if (nxvol_tet > 0) {
            tet_vol_layer = (int *)SMESH_ALLOC((size_t)nxvol_tet * sizeof(int));
            tet_vol_trank = (int *)SMESH_ALLOC((size_t)nxvol_tet * sizeof(int));
            hier_fill_vol_layers(TET4, level, nlevels, levels, tet_vol_layer);
            hier_slot_ranks(nxvol_tet, tet_vol_layer, nlevels, tet_vol_trank, n_tet_vol_t);
        }
        layer_base[0] = 0;
        layer_base[1] = (large_idx_t)n_coarse_global;
        for (int k = 1; k < nlevels; ++k) {
            layer_base[k + 1] = layer_base[k] + (large_idx_t)n_edges_global * (large_idx_t)n_edge_t[k] +
                                (large_idx_t)n_hex_faces_global * (large_idx_t)n_hex_face_t[k] +
                                (large_idx_t)n_tet_faces_global * (large_idx_t)n_tet_face_t[k] +
                                (large_idx_t)n_hex_elem_global * (large_idx_t)n_hex_vol_t[k] +
                                (large_idx_t)n_tet_elem_global * (large_idx_t)n_tet_vol_t[k];
        }
        if ((ptrdiff_t)layer_base[nlevels] != n_ss_global) {
            fprintf(stderr,
                    "to_semistructured: hierarchical SS gid count %ld does not match A8 count %ld\n",
                    (long)layer_base[nlevels],
                    (long)n_ss_global);
            return nullptr;
        }
    }

    int *edge_uo         = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *edge_ua         = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *hex_face_uo     = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_hex_face_uniq, 1), sizeof(int));
    int *hex_face_ua     = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_hex_face_uniq, 1), sizeof(int));
    int *tet_face_uo     = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_tet_face_uniq, 1), sizeof(int));
    int *tet_face_ua     = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_tet_face_uniq, 1), sizeof(int));
    int *hex_vol_uo      = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_hex, 1), sizeof(int));
    int *hex_vol_ua      = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_hex, 1), sizeof(int));
    int *tet_vol_uo      = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tet, 1), sizeof(int));
    int *tet_vol_ua      = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tet, 1), sizeof(int));

    ie  = 0;
    ihf = 0;
    itf = 0;
    ihv = 0;
    itv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const ptrdiff_t n_owned = mesh->block((size_t)b)->n_elements_owned();
        const int       n_macro = (block_types[b] == HEX8) ? 8 : 4;
        const int       nsides  = (block_types[b] == HEX8) ? 6 : 4;
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    const int *conn = (block_types[b] == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
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
            if (block_types[b] == HEX8 && nxf_hex > 0) {
                for (int f = 0; f < nsides; ++f) {
                    const ptrdiff_t u = hex_face_inc_to_uniq[ihf++];
                    if (from_owned) {
                        hex_face_uo[u] = 1;
                    } else {
                        hex_face_ua[u] = 1;
                    }
                }
            } else if (block_types[b] == TET4 && nxf_tet > 0) {
                for (int f = 0; f < nsides; ++f) {
                    const ptrdiff_t u = tet_face_inc_to_uniq[itf++];
                    if (from_owned) {
                        tet_face_uo[u] = 1;
                    } else {
                        tet_face_ua[u] = 1;
                    }
                }
            }
            if (block_types[b] == HEX8 && nxvol_hex > 0) {
                if (from_owned) {
                    hex_vol_uo[ihv] = 1;
                } else {
                    hex_vol_ua[ihv] = 1;
                }
                ihv++;
            } else if (block_types[b] == TET4 && nxvol_tet > 0) {
                if (from_owned) {
                    tet_vol_uo[itv] = 1;
                } else {
                    tet_vol_ua[itv] = 1;
                }
                itv++;
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
    if (nxf_hex > 0) {
        count_entity_nodes(n_hex_face_uniq, nxf_hex, hex_face_owner, hex_face_shared, hex_face_uo, hex_face_ua, rank, n_bkt);
    }
    if (nxf_tet > 0) {
        count_entity_nodes(n_tet_face_uniq, nxf_tet, tet_face_owner, tet_face_shared, tet_face_uo, tet_face_ua, rank, n_bkt);
    }
    if (nxvol_hex > 0) {
        count_entity_nodes(n_e_hex, nxvol_hex, hex_vol_owner, hex_vol_shared, hex_vol_uo, hex_vol_ua, rank, n_bkt);
    }
    if (nxvol_tet > 0) {
        count_entity_nodes(n_e_tet, nxvol_tet, tet_vol_owner, tet_vol_shared, tet_vol_uo, tet_vol_ua, rank, n_bkt);
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

    auto         node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto         node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap         = node_mapping->data();
    int         *nown         = node_owner->data();
    auto         points       = create_host_buffer<geom_t>((size_t)mesh->spatial_dimension(), (size_t)n_local);
    auto         coarse_p     = mesh->points()->data();
    auto         p            = points->data();
    const int    sdim         = mesh->spatial_dimension();

    idx_t *corner_ss    = (idx_t *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(idx_t));
    idx_t *edge_ss      = (nxedge > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_edge_uniq * (size_t)nxedge * sizeof(idx_t)) : nullptr;
    idx_t *hex_face_ss  = (nxf_hex > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_hex_face_uniq * (size_t)nxf_hex * sizeof(idx_t)) : nullptr;
    idx_t *tet_face_ss  = (nxf_tet > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_tet_face_uniq * (size_t)nxf_tet * sizeof(idx_t)) : nullptr;
    idx_t *hex_vol_ss   = (nxvol_hex > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_e_hex * (size_t)nxvol_hex * sizeof(idx_t)) : nullptr;
    idx_t *tet_vol_ss   = (nxvol_tet > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_e_tet * (size_t)nxvol_tet * sizeof(idx_t)) : nullptr;

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

    large_idx_t *edge_node_gid     = alloc_entity_node_gids(n_edge_uniq, nxedge);
    large_idx_t *hex_face_node_gid = alloc_entity_node_gids(n_hex_face_uniq, nxf_hex);
    large_idx_t *tet_face_node_gid = alloc_entity_node_gids(n_tet_face_uniq, nxf_tet);
    large_idx_t *hex_vol_node_gid  = alloc_entity_node_gids(n_e_hex, nxvol_hex);
    large_idx_t *tet_vol_node_gid  = alloc_entity_node_gids(n_e_tet, nxvol_tet);
    if (nxedge > 0) {
        if (hierarchical_ordering) {
            fill_hier_mixed_node_gids(n_edge_uniq,
                                      nxedge,
                                      0,
                                      edge_gid,
                                      edge_layer,
                                      edge_trank,
                                      n_edge_t,
                                      n_hex_face_t,
                                      n_tet_face_t,
                                      n_hex_vol_t,
                                      n_tet_vol_t,
                                      layer_base,
                                      n_edges_global,
                                      n_hex_faces_global,
                                      n_tet_faces_global,
                                      n_hex_elem_global,
                                      edge_node_gid);
        } else {
            fill_flat_node_gids(n_edge_uniq, nxedge, edge_base, edge_gid, edge_node_gid);
        }
        pack_entity_nodes(n_edge_uniq, nxedge, edge_node_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    }
    if (nxf_hex > 0) {
        if (hierarchical_ordering) {
            fill_hier_mixed_node_gids(n_hex_face_uniq,
                                      nxf_hex,
                                      1,
                                      hex_face_gid,
                                      hex_face_layer,
                                      hex_face_trank,
                                      n_edge_t,
                                      n_hex_face_t,
                                      n_tet_face_t,
                                      n_hex_vol_t,
                                      n_tet_vol_t,
                                      layer_base,
                                      n_edges_global,
                                      n_hex_faces_global,
                                      n_tet_faces_global,
                                      n_hex_elem_global,
                                      hex_face_node_gid);
        } else {
            fill_flat_node_gids(n_hex_face_uniq, nxf_hex, hex_face_base, hex_face_gid, hex_face_node_gid);
        }
        pack_entity_nodes(n_hex_face_uniq,
                          nxf_hex,
                          hex_face_node_gid,
                          hex_face_owner,
                          hex_face_shared,
                          hex_face_uo,
                          hex_face_ua,
                          rank,
                          cur,
                          nmap,
                          nown,
                          hex_face_ss);
    }
    if (nxf_tet > 0) {
        if (hierarchical_ordering) {
            fill_hier_mixed_node_gids(n_tet_face_uniq,
                                      nxf_tet,
                                      2,
                                      tet_face_gid,
                                      tet_face_layer,
                                      tet_face_trank,
                                      n_edge_t,
                                      n_hex_face_t,
                                      n_tet_face_t,
                                      n_hex_vol_t,
                                      n_tet_vol_t,
                                      layer_base,
                                      n_edges_global,
                                      n_hex_faces_global,
                                      n_tet_faces_global,
                                      n_hex_elem_global,
                                      tet_face_node_gid);
        } else {
            fill_flat_node_gids(n_tet_face_uniq, nxf_tet, tet_face_base, tet_face_gid, tet_face_node_gid);
        }
        pack_entity_nodes(n_tet_face_uniq,
                          nxf_tet,
                          tet_face_node_gid,
                          tet_face_owner,
                          tet_face_shared,
                          tet_face_uo,
                          tet_face_ua,
                          rank,
                          cur,
                          nmap,
                          nown,
                          tet_face_ss);
    }
    if (nxvol_hex > 0) {
        if (hierarchical_ordering) {
            fill_hier_mixed_node_gids(n_e_hex,
                                      nxvol_hex,
                                      3,
                                      hex_vol_gid,
                                      hex_vol_layer,
                                      hex_vol_trank,
                                      n_edge_t,
                                      n_hex_face_t,
                                      n_tet_face_t,
                                      n_hex_vol_t,
                                      n_tet_vol_t,
                                      layer_base,
                                      n_edges_global,
                                      n_hex_faces_global,
                                      n_tet_faces_global,
                                      n_hex_elem_global,
                                      hex_vol_node_gid);
        } else {
            fill_flat_node_gids(n_e_hex, nxvol_hex, hex_vol_base, hex_vol_gid, hex_vol_node_gid);
        }
        pack_entity_nodes(n_e_hex, nxvol_hex, hex_vol_node_gid, hex_vol_owner, hex_vol_shared, hex_vol_uo, hex_vol_ua, rank, cur, nmap, nown, hex_vol_ss);
    }
    if (nxvol_tet > 0) {
        if (hierarchical_ordering) {
            fill_hier_mixed_node_gids(n_e_tet,
                                      nxvol_tet,
                                      4,
                                      tet_vol_gid,
                                      tet_vol_layer,
                                      tet_vol_trank,
                                      n_edge_t,
                                      n_hex_face_t,
                                      n_tet_face_t,
                                      n_hex_vol_t,
                                      n_tet_vol_t,
                                      layer_base,
                                      n_edges_global,
                                      n_hex_faces_global,
                                      n_tet_faces_global,
                                      n_hex_elem_global,
                                      tet_vol_node_gid);
        } else {
            fill_flat_node_gids(n_e_tet, nxvol_tet, tet_vol_base, tet_vol_gid, tet_vol_node_gid);
        }
        pack_entity_nodes(n_e_tet, nxvol_tet, tet_vol_node_gid, tet_vol_owner, tet_vol_shared, tet_vol_uo, tet_vol_ua, rank, cur, nmap, nown, tet_vol_ss);
    }
    SMESH_FREE(edge_node_gid);
    SMESH_FREE(hex_face_node_gid);
    SMESH_FREE(tet_face_node_gid);
    SMESH_FREE(hex_vol_node_gid);
    SMESH_FREE(tet_vol_node_gid);
    SMESH_FREE(levels);
    SMESH_FREE(edge_layer);
    SMESH_FREE(edge_trank);
    SMESH_FREE(hex_face_layer);
    SMESH_FREE(hex_face_trank);
    SMESH_FREE(tet_face_layer);
    SMESH_FREE(tet_face_trank);
    SMESH_FREE(hex_vol_layer);
    SMESH_FREE(hex_vol_trank);
    SMESH_FREE(tet_vol_layer);
    SMESH_FREE(tet_vol_trank);
    SMESH_FREE(n_edge_t);
    SMESH_FREE(n_hex_face_t);
    SMESH_FREE(n_tet_face_t);
    SMESH_FREE(n_hex_vol_t);
    SMESH_FREE(n_tet_vol_t);
    SMESH_FREE(layer_base);

    int *hex_coords[3] = {nullptr, nullptr, nullptr};
    int *tet_coords[3] = {nullptr, nullptr, nullptr};
    hex_fill_lattice(level, hex_coords);
    tet_fill_lattice(level, tet_coords);

    std::vector<std::shared_ptr<Mesh::Block>> ss_blocks((size_t)n_blocks);
    ie  = 0;
    ihf = 0;
    itf = 0;
    ihv = 0;
    itv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            coarse_block = mesh->block((size_t)b);
        const int       n_macro      = (block_types[b] == HEX8) ? 8 : 4;
        const int       nsides       = (block_types[b] == HEX8) ? 6 : 4;
        const int       nxe          = (block_types[b] == HEX8) ? sshex8_nxe(level) : sstet4_nxe(level);
        const int      *lagr         = (block_types[b] == HEX8) ? hex_corners : tet_corners;
        auto            ss_elems     = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_e[b]);
        idx_t         **out          = ss_elems->data();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                const idx_t local = coarse_soa[b][d][e];
                gc[d]             = coarse_nmap[local];
                out[lagr[d]][e]   = corner_ss[local];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    const int *conn = (block_types[b] == HEX8) ? hex_lagr_conn[d1] : tet_lagr_conn[d1];
                    for (int k = 0; k < 3; ++k) {
                        const int d2 = conn[k];
                        if (gc[d1] > gc[d2]) {
                            continue;
                        }
                        const idx_t estart = edge_ss[edge_inc_to_uniq[ie++] * nxedge];
                        if (block_types[b] == HEX8) {
                            write_hex_edge(level, hex_corners, hex_coords, d1, d2, estart, e, out);
                        } else {
                            write_tet_edge(level, tet_corners, tet_coords, d1, d2, estart, e, out);
                        }
                    }
                }
            }
            if (block_types[b] == HEX8 && nxf_hex > 0) {
                for (int f = 0; f < nsides; ++f) {
                    const idx_t foff = hex_face_ss[hex_face_inc_to_uniq[ihf++] * nxf_hex];
                    write_hex_face(level, gc, lst_hex, hex_corners, hex_coords, foff, e, f, out);
                }
            } else if (block_types[b] == TET4 && nxf_tet > 0) {
                for (int f = 0; f < nsides; ++f) {
                    const idx_t foff = tet_face_ss[tet_face_inc_to_uniq[itf++] * nxf_tet];
                    write_tet_face(level, gc, lst_tet, tet_corners, tet_coords, foff, e, f, out);
                }
            }
            if (block_types[b] == HEX8 && nxvol_hex > 0) {
                const idx_t voff = hex_vol_ss[ihv * nxvol_hex];
                const int   Lm1  = level - 1;
                for (int zi = 1; zi < level; ++zi) {
                    for (int yi = 1; yi < level; ++yi) {
                        for (int xi = 1; xi < level; ++xi) {
                            const int lidx = sshex8_lidx(level, xi, yi, zi);
                            const int en   = (zi - 1) * Lm1 * Lm1 + (yi - 1) * Lm1 + (xi - 1);
                            out[lidx][e]   = voff + (idx_t)en;
                        }
                    }
                }
                ihv++;
            } else if (block_types[b] == TET4 && nxvol_tet > 0) {
                const idx_t voff = tet_vol_ss[itv * nxvol_tet];
                for (int z = 1; z <= level - 3; ++z) {
                    for (int y = 1; y <= level - 2 - z; ++y) {
                        for (int x = 1; x <= level - 1 - z - y; ++x) {
                            const int lidx = sstet4_lidx(level, x, y, z);
                            const int en   = sstet4_lidx(level - 4, x - 1, y - 1, z - 1);
                            out[lidx][e]   = voff + (idx_t)en;
                        }
                    }
                }
                itv++;
            }
        }
        auto ss_block = std::make_shared<Mesh::Block>();
        ss_block->set_name(coarse_block->name());
        ss_block->set_element_type(semistructured_type(block_types[b], level));
        ss_block->set_elements(ss_elems);
        ss_block->set_distributed_elements(coarse_block->n_elements_owned(),
                                           coarse_block->n_elements_shared(),
                                           coarse_block->n_elements_ghosts(),
                                           coarse_block->element_mapping(),
                                           coarse_block->aura_element_mapping());
        ss_blocks[(size_t)b] = ss_block;
    }
    for (int d = 0; d < 3; ++d) {
        SMESH_FREE(hex_coords[d]);
        SMESH_FREE(tet_coords[d]);
    }

    auto ghosts_and_aura = create_host_buffer<idx_t>((size_t)(n_ghosts + n_aura));
    auto node_offsets    = create_host_buffer<ptrdiff_t>((size_t)size + 1);
    node_ownership_ranges(comm->get(), n_owned, node_offsets->data());
    SMESH_ASSERT(node_offsets->data()[size] == n_ss_global);

    if (n_ghosts + n_aura > 0) {
        const ptrdiff_t n_id        = rank_split(n_ss_global, size, rank);
        idx_t          *global2owned = (idx_t *)SMESH_CALLOC((size_t)n_id, sizeof(idx_t));
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

    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        idx_t **el = ss_blocks[(size_t)b]->elements()->data();
        if (block_types[b] == HEX8) {
            sshex8_fill_points(level, n_e[b], el, p, p);
        } else {
            sstet4_fill_points(level, n_e[b], el, p, p);
        }
    }

    SMESH_FREE(block_types);
    SMESH_FREE(n_owned_b);
    SMESH_FREE(n_global_b);
    SMESH_FREE(hex_concat0);
    SMESH_FREE(tet_concat0);
    SMESH_FREE(n_e);
    SMESH_FREE(coarse_soa);
    SMESH_FREE(c_uo);
    SMESH_FREE(c_ua);
    SMESH_FREE(edge_inc_to_uniq);
    SMESH_FREE(edge_gid);
    SMESH_FREE(edge_owner);
    SMESH_FREE(edge_shared);
    SMESH_FREE(hex_face_inc_to_uniq);
    SMESH_FREE(hex_face_gid);
    SMESH_FREE(hex_face_owner);
    SMESH_FREE(hex_face_shared);
    SMESH_FREE(tet_face_inc_to_uniq);
    SMESH_FREE(tet_face_gid);
    SMESH_FREE(tet_face_owner);
    SMESH_FREE(tet_face_shared);
    SMESH_FREE(hex_vol_gid);
    SMESH_FREE(hex_vol_owner);
    SMESH_FREE(hex_vol_shared);
    SMESH_FREE(tet_vol_gid);
    SMESH_FREE(tet_vol_owner);
    SMESH_FREE(tet_vol_shared);
    SMESH_FREE(edge_uo);
    SMESH_FREE(edge_ua);
    SMESH_FREE(hex_face_uo);
    SMESH_FREE(hex_face_ua);
    SMESH_FREE(tet_face_uo);
    SMESH_FREE(tet_face_ua);
    SMESH_FREE(hex_vol_uo);
    SMESH_FREE(hex_vol_ua);
    SMESH_FREE(tet_vol_uo);
    SMESH_FREE(tet_vol_ua);
    SMESH_FREE(corner_ss);
    SMESH_FREE(edge_ss);
    SMESH_FREE(hex_face_ss);
    SMESH_FREE(tet_face_ss);
    SMESH_FREE(hex_vol_ss);
    SMESH_FREE(tet_vol_ss);

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
