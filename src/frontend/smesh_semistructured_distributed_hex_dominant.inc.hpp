std::shared_ptr<Mesh> to_semistructured_distributed_hex_dominant(const int                    level,
                                                                 const std::shared_ptr<Mesh> &mesh,
                                                                 const bool                   hierarchical_ordering) {
    SMESH_TRACE_SCOPE("to_semistructured_distributed_hex_dominant");
    if (level < 1) {
        fprintf(stderr, "to_semistructured: mixed HEX-dominant SS requires level >= 1\n");
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

    const int nxedge     = (level > 1) ? (level - 1) : 0;
    const int nxf_tri    = sstet4_nxface(level);
    const int nxf_quad   = (level > 1) ? ((level - 1) * (level - 1)) : 0;
    const int nxvol_hex  = (level - 1) * (level - 1) * (level - 1);
    const int nxvol_tet  = sstet4_nxvol(level);
    const int nxvol_wedge = sswedge_nxvol(level);
    const int nxvol_pyr  = sspyramid_nxvol(level);

    static const int hex_conn[8][3]   = {{1, 3, 4}, {0, 2, 5}, {1, 3, 6}, {0, 2, 7}, {0, 5, 7}, {1, 4, 6}, {2, 5, 7}, {3, 4, 6}};
    static const int tet_conn[4][3]   = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
    static const int wedge_conn[6][3] = {{1, 2, 3}, {0, 2, 4}, {0, 1, 5}, {0, 4, 5}, {1, 3, 5}, {2, 3, 4}};
    static const int pyr_nneigh[5]    = {3, 3, 3, 3, 4};
    static const int pyr_conn[5][4]   = {{1, 3, 4, -1}, {0, 2, 4, -1}, {1, 3, 4, -1}, {0, 2, 4, -1}, {0, 1, 2, 3}};
    static const int wedge_xyz[6][3]  = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
    static const int pyr_ijk[5][3]    = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}};
    static const int tet_xyz[4][3]    = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    auto n_macro_of = [](const enum ElemType f) {
        if (f == HEX8) return 8;
        if (f == TET4) return 4;
        if (f == WEDGE6) return 6;
        return 5;
    };
    auto nsides_of = [](const enum ElemType f) {
        if (f == HEX8) return 6;
        if (f == TET4) return 4;
        return 5;
    };
    auto nxe_of = [&](const enum ElemType f) {
        if (f == HEX8) return sshex8_nxe(level);
        if (f == TET4) return sstet4_nxe(level);
        if (f == WEDGE6) return sswedge_nxe(level);
        return sspyramid_nxe(level);
    };
    auto nxvol_of = [&](const enum ElemType f) {
        if (f == HEX8) return nxvol_hex;
        if (f == TET4) return nxvol_tet;
        if (f == WEDGE6) return nxvol_wedge;
        return nxvol_pyr;
    };
    auto lidx_of = [](const enum ElemType f) {
        if (f == HEX8) return lidx_hex;
        if (f == TET4) return lidx_tet;
        if (f == WEDGE6) return lidx_wedge;
        return lidx_pyr;
    };
    auto nneigh = [&](const enum ElemType f, const int d1) {
        if (f == PYRAMID5) return pyr_nneigh[d1];
        return 3;
    };
    auto neigh = [&](const enum ElemType f, const int d1, const int k) {
        if (f == HEX8) return hex_conn[d1][k];
        if (f == TET4) return tet_conn[d1][k];
        if (f == WEDGE6) return wedge_conn[d1][k];
        return pyr_conn[d1][k];
    };
    auto side_nnxs = [&](const enum ElemType f, const int s) {
        if (f == HEX8) return 4;
        if (f == TET4) return 3;
        if (f == WEDGE6) return s < 3 ? 4 : 3;
        return s < 4 ? 3 : 4;
    };

    const ptrdiff_t n_blocks = (ptrdiff_t)mesh->n_blocks();
    enum ElemType *block_types = (enum ElemType *)SMESH_ALLOC((size_t)n_blocks * sizeof(enum ElemType));
    ptrdiff_t     *n_owned_b   = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    ptrdiff_t     *n_global_b  = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto                block  = mesh->block((size_t)b);
        const enum ElemType family = ss_source_family(mesh->element_type(static_cast<block_idx_t>(b)));
        block_types[b]             = family;
        n_owned_b[b]               = block->n_elements_owned();
        if (block->n_nodes_per_element() != n_macro_of(family)) {
            fprintf(stderr, "to_semistructured: block '%s' node count mismatch\n", block->name().c_str());
            return nullptr;
        }
    }
    SMESH_MPI_CATCH(MPI_Allreduce(n_owned_b, n_global_b, (int)n_blocks, mpi_type<ptrdiff_t>(), MPI_SUM, comm->get()));

    large_idx_t concat0[4] = {0, 0, 0, 0};
    large_idx_t acc[4]     = {0, 0, 0, 0};
    auto fam_slot = [](const enum ElemType f) {
        if (f == HEX8) return 0;
        if (f == TET4) return 1;
        if (f == WEDGE6) return 2;
        return 3;
    };
    large_idx_t *block_concat = (large_idx_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(large_idx_t));
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const int s      = fam_slot(block_types[b]);
        block_concat[b]  = acc[s];
        acc[s]          += (large_idx_t)n_global_b[b];
    }
    const ptrdiff_t n_hex_elem_global   = (ptrdiff_t)acc[0];
    const ptrdiff_t n_tet_elem_global   = (ptrdiff_t)acc[1];
    const ptrdiff_t n_wedge_elem_global = (ptrdiff_t)acc[2];
    const ptrdiff_t n_pyr_elem_global   = (ptrdiff_t)acc[3];
    SMESH_UNUSED(concat0);

    ptrdiff_t           *n_e        = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    const idx_t *const **coarse_soa = (const idx_t *const **)SMESH_ALLOC((size_t)n_blocks * sizeof(const idx_t *const *));
    ptrdiff_t            n_e_fam[4] = {0, 0, 0, 0};
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        n_e[b]        = mesh->block((size_t)b)->n_elements();
        coarse_soa[b] = mesh->block((size_t)b)->elements()->data();
        n_e_fam[fam_slot(block_types[b])] += n_e[b];
    }

    LocalSideTable lst[4];
    lst[0].fill(HEX8);
    lst[1].fill(TET4);
    lst[2].fill(WEDGE6);
    lst[3].fill(PYRAMID5);

    int hex_corners[8] = {sshex8_lidx(level, 0, 0, 0),
                          sshex8_lidx(level, level, 0, 0),
                          sshex8_lidx(level, level, level, 0),
                          sshex8_lidx(level, 0, level, 0),
                          sshex8_lidx(level, 0, 0, level),
                          sshex8_lidx(level, level, 0, level),
                          sshex8_lidx(level, level, level, level),
                          sshex8_lidx(level, 0, level, level)};
    int tet_corners[4];
    int wedge_corners[6];
    int pyr_corners[5];
    for (int d = 0; d < 4; ++d) {
        tet_corners[d] = sstet4_lidx(level, tet_xyz[d][0] * level, tet_xyz[d][1] * level, tet_xyz[d][2] * level);
    }
    for (int d = 0; d < 6; ++d) {
        wedge_corners[d] = sswedge_lidx(level, wedge_xyz[d][0] * level, wedge_xyz[d][1] * level, wedge_xyz[d][2] * level);
    }
    for (int d = 0; d < 5; ++d) {
        pyr_corners[d] = (d == 4) ? sspyramid_lidx(level, 0, 0, level)
                                  : sspyramid_lidx(level, pyr_ijk[d][0] * level, pyr_ijk[d][1] * level, pyr_ijk[d][2] * level);
    }
    auto corners_of = [&](const enum ElemType f) -> int * {
        if (f == HEX8) return hex_corners;
        if (f == TET4) return tet_corners;
        if (f == WEDGE6) return wedge_corners;
        return pyr_corners;
    };

    int *c_uo = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));
    int *c_ua = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));

    ptrdiff_t n_edge_inc = 0, n_tri_inc = 0, n_quad_inc = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const enum ElemType f = block_types[b];
        const int n_macro     = n_macro_of(f);
        const int nsides      = nsides_of(f);
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    for (int k = 0; k < nneigh(f, d1); ++k) {
                        if (gc[d1] <= gc[neigh(f, d1, k)]) {
                            n_edge_inc++;
                        }
                    }
                }
            }
            for (int s = 0; s < nsides; ++s) {
                if (side_nnxs(f, s) == 3) {
                    n_tri_inc += (nxf_tri > 0) ? 1 : 0;
                } else {
                    n_quad_inc += (nxf_quad > 0) ? 1 : 0;
                }
            }
        }
    }

    large_idx_t *edge_keys = nullptr, *edge_aux = nullptr;
    idx_t       *edge_loc  = nullptr;
    if (nxedge > 0 && n_edge_inc > 0) {
        edge_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_edge_inc * 4 * sizeof(large_idx_t));
        edge_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_edge_inc * sizeof(large_idx_t));
        edge_loc  = (idx_t *)SMESH_ALLOC((size_t)n_edge_inc * sizeof(idx_t));
    }
    large_idx_t *tri_keys = nullptr, *tri_aux = nullptr;
    idx_t       *tri_loc  = nullptr;
    if (nxf_tri > 0 && n_tri_inc > 0) {
        tri_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_tri_inc * 4 * sizeof(large_idx_t));
        tri_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_tri_inc * sizeof(large_idx_t));
        tri_loc  = (idx_t *)SMESH_ALLOC((size_t)n_tri_inc * sizeof(idx_t));
    }
    large_idx_t *quad_keys = nullptr, *quad_aux = nullptr;
    idx_t       *quad_loc  = nullptr;
    if (nxf_quad > 0 && n_quad_inc > 0) {
        quad_keys = (large_idx_t *)SMESH_ALLOC((size_t)n_quad_inc * 4 * sizeof(large_idx_t));
        quad_aux  = (large_idx_t *)SMESH_ALLOC((size_t)n_quad_inc * sizeof(large_idx_t));
        quad_loc  = (idx_t *)SMESH_ALLOC((size_t)n_quad_inc * sizeof(idx_t));
    }
    large_idx_t *vol_ids[4] = {nullptr, nullptr, nullptr, nullptr};
    large_idx_t *vol_aux[4] = {nullptr, nullptr, nullptr, nullptr};
    for (int s = 0; s < 4; ++s) {
        const int nxv = (s == 0) ? nxvol_hex : (s == 1) ? nxvol_tet : (s == 2) ? nxvol_wedge : nxvol_pyr;
        if (nxv > 0 && n_e_fam[s] > 0) {
            vol_ids[s] = (large_idx_t *)SMESH_ALLOC((size_t)n_e_fam[s] * sizeof(large_idx_t));
            vol_aux[s] = (large_idx_t *)SMESH_ALLOC((size_t)n_e_fam[s] * sizeof(large_idx_t));
        }
    }

    ptrdiff_t ie = 0, itri = 0, iq = 0, iv[4] = {0, 0, 0, 0};
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto                block   = mesh->block((size_t)b);
        const enum ElemType f       = block_types[b];
        const int           n_macro = n_macro_of(f);
        const int           nsides  = nsides_of(f);
        const int           fs      = fam_slot(f);
        const ptrdiff_t     n_owned = block->n_elements_owned();
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
                    for (int k = 0; k < nneigh(f, d1); ++k) {
                        const int d2 = neigh(f, d1, k);
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
            for (int s = 0; s < nsides; ++s) {
                const int nnxs = side_nnxs(f, s);
                if (nnxs == 3 && nxf_tri > 0) {
                    large_idx_t fk[4] = {k_key_pad, k_key_pad, k_key_pad, k_key_pad};
                    idx_t       loc_min = lc[lst[fs](s, 0)];
                    large_idx_t gmin    = gc[lst[fs](s, 0)];
                    for (int i = 0; i < 3; ++i) {
                        fk[i] = gc[lst[fs](s, i)];
                        if (gc[lst[fs](s, i)] < gmin) {
                            gmin    = gc[lst[fs](s, i)];
                            loc_min = lc[lst[fs](s, i)];
                        }
                    }
                    sort4(fk, 3);
                    for (int i = 0; i < 4; ++i) {
                        tri_keys[itri * 4 + i] = fk[i];
                    }
                    tri_aux[itri] = owned_pref_eid_aux(from_owned, element_gid(*block, block_concat[b], e), n_elem_global);
                    tri_loc[itri] = loc_min;
                    itri++;
                } else if (nnxs == 4 && nxf_quad > 0) {
                    large_idx_t fk[4] = {k_key_pad, k_key_pad, k_key_pad, k_key_pad};
                    idx_t       loc_min = lc[lst[fs](s, 0)];
                    large_idx_t gmin    = gc[lst[fs](s, 0)];
                    for (int i = 0; i < 4; ++i) {
                        fk[i] = gc[lst[fs](s, i)];
                        if (gc[lst[fs](s, i)] < gmin) {
                            gmin    = gc[lst[fs](s, i)];
                            loc_min = lc[lst[fs](s, i)];
                        }
                    }
                    sort4(fk, 4);
                    for (int i = 0; i < 4; ++i) {
                        quad_keys[iq * 4 + i] = fk[i];
                    }
                    quad_aux[iq] = owned_pref_eid_aux(from_owned, element_gid(*block, block_concat[b], e), n_elem_global);
                    quad_loc[iq] = loc_min;
                    iq++;
                }
            }
            if (nxvol_of(f) > 0) {
                vol_ids[fs][iv[fs]] = element_gid(*block, block_concat[b], e);
                vol_aux[fs][iv[fs]] = from_owned ? 0 : 1;
                iv[fs]++;
            }
        }
    }

    ptrdiff_t    n_edge_uniq = 0, n_tri_uniq = 0, n_quad_uniq = 0;
    ptrdiff_t   *edge_inc_to_uniq = nullptr, *tri_inc_to_uniq = nullptr, *quad_inc_to_uniq = nullptr;
    large_idx_t *edge_gid = nullptr, *tri_gid = nullptr, *quad_gid = nullptr;
    int         *edge_owner = nullptr, *tri_owner = nullptr, *quad_owner = nullptr;
    int         *edge_shared = nullptr, *tri_shared = nullptr, *quad_shared = nullptr;
    ptrdiff_t    n_edges_global = 0, n_tri_global = 0, n_quad_global = 0;
    if (nxedge > 0) {
        if (unique_inc_tuples(comm->get(), n_coarse_global, n_edge_inc, edge_keys, edge_aux, edge_loc, n_coarse_local,
                              &n_edge_uniq, &edge_inc_to_uniq, &edge_gid, &edge_owner, &edge_shared, &n_edges_global) !=
            SMESH_SUCCESS) {
            return nullptr;
        }
    }
    if (nxf_tri > 0) {
        if (unique_inc_tuples(comm->get(), n_coarse_global, n_tri_inc, tri_keys, tri_aux, tri_loc, n_coarse_local,
                              &n_tri_uniq, &tri_inc_to_uniq, &tri_gid, &tri_owner, &tri_shared, &n_tri_global) !=
            SMESH_SUCCESS) {
            return nullptr;
        }
    }
    if (nxf_quad > 0) {
        if (unique_inc_tuples(comm->get(), n_coarse_global, n_quad_inc, quad_keys, quad_aux, quad_loc, n_coarse_local,
                              &n_quad_uniq, &quad_inc_to_uniq, &quad_gid, &quad_owner, &quad_shared, &n_quad_global) !=
            SMESH_SUCCESS) {
            return nullptr;
        }
    }
    large_idx_t *vol_gid[4]    = {nullptr, nullptr, nullptr, nullptr};
    int         *vol_owner[4]  = {nullptr, nullptr, nullptr, nullptr};
    int         *vol_shared[4] = {nullptr, nullptr, nullptr, nullptr};
    const ptrdiff_t n_elem_g[4] = {n_hex_elem_global, n_tet_elem_global, n_wedge_elem_global, n_pyr_elem_global};
    const int       nxvol[4]    = {nxvol_hex, nxvol_tet, nxvol_wedge, nxvol_pyr};
    for (int s = 0; s < 4; ++s) {
        if (nxvol[s] > 0) {
            vol_gid[s]    = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_fam[s], 1) * sizeof(large_idx_t));
            vol_owner[s]  = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_fam[s], 1) * sizeof(int));
            vol_shared[s] = (int *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_e_fam[s], 1) * sizeof(int));
            if (unique_by_id(comm->get(), n_elem_g[s], n_e_fam[s], vol_ids[s], vol_aux[s], vol_gid[s], vol_owner[s], vol_shared[s]) !=
                SMESH_SUCCESS) {
                return nullptr;
            }
            SMESH_FREE(vol_ids[s]);
            SMESH_FREE(vol_aux[s]);
        }
    }

    const large_idx_t edge_base = (large_idx_t)n_coarse_global;
    const large_idx_t tri_base  = edge_base + (large_idx_t)n_edges_global * (large_idx_t)nxedge;
    const large_idx_t quad_base = tri_base + (large_idx_t)n_tri_global * (large_idx_t)nxf_tri;
    large_idx_t       vol_base[4];
    vol_base[0] = quad_base + (large_idx_t)n_quad_global * (large_idx_t)nxf_quad;
    vol_base[1] = vol_base[0] + (large_idx_t)n_hex_elem_global * (large_idx_t)nxvol_hex;
    vol_base[2] = vol_base[1] + (large_idx_t)n_tet_elem_global * (large_idx_t)nxvol_tet;
    vol_base[3] = vol_base[2] + (large_idx_t)n_wedge_elem_global * (large_idx_t)nxvol_wedge;
    const ptrdiff_t n_ss_global =
            (ptrdiff_t)vol_base[3] + n_pyr_elem_global * (ptrdiff_t)nxvol_pyr;
    if (n_ss_global < size) {
        fprintf(stderr, "to_semistructured: SS node count smaller than communicator size\n");
        return nullptr;
    }

    const int nlevels = hierarchical_ordering ? sshex8_hierarchical_n_levels(level) : 0;
    int *levels = nullptr, *edge_layer = nullptr, *edge_trank = nullptr;
    int *tri_layer = nullptr, *tri_trank = nullptr, *quad_layer = nullptr, *quad_trank = nullptr;
    int *vol_layer[4] = {nullptr, nullptr, nullptr, nullptr};
    int *vol_trank[4] = {nullptr, nullptr, nullptr, nullptr};
    int *n_edge_t = nullptr, *n_tri_t = nullptr, *n_quad_t = nullptr;
    int *n_vol_t[4] = {nullptr, nullptr, nullptr, nullptr};
    large_idx_t *layer_base = nullptr;
    if (hierarchical_ordering) {
        levels   = (int *)SMESH_ALLOC((size_t)nlevels * sizeof(int));
        n_edge_t = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_tri_t  = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_quad_t = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        for (int s = 0; s < 4; ++s) {
            n_vol_t[s] = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        }
        layer_base = (large_idx_t *)SMESH_ALLOC((size_t)(nlevels + 1) * sizeof(large_idx_t));
        sshex8_hierarchical_mesh_levels(level, nlevels, levels);
        if (nxedge > 0) {
            edge_layer = (int *)SMESH_ALLOC((size_t)nxedge * sizeof(int));
            edge_trank = (int *)SMESH_ALLOC((size_t)nxedge * sizeof(int));
            hier_fill_edge_layers(level, nxedge, nlevels, levels, edge_layer);
            hier_slot_ranks(nxedge, edge_layer, nlevels, edge_trank, n_edge_t);
        }
        if (nxf_tri > 0) {
            tri_layer = (int *)SMESH_ALLOC((size_t)nxf_tri * sizeof(int));
            tri_trank = (int *)SMESH_ALLOC((size_t)nxf_tri * sizeof(int));
            int t     = 0;
            for (int tt = 1; tt <= level - 2; ++tt) {
                for (int s = 1; s <= level - 1 - tt; ++s) {
                    tri_layer[t++] = hier_first_layer(level, nlevels, levels, level - s - tt, s, tt);
                }
            }
            hier_slot_ranks(nxf_tri, tri_layer, nlevels, tri_trank, n_tri_t);
        }
        if (nxf_quad > 0) {
            quad_layer = (int *)SMESH_ALLOC((size_t)nxf_quad * sizeof(int));
            quad_trank = (int *)SMESH_ALLOC((size_t)nxf_quad * sizeof(int));
            int t      = 0;
            for (int tt = 1; tt < level; ++tt) {
                for (int s = 1; s < level; ++s) {
                    quad_layer[t++] = hier_first_layer(level, nlevels, levels, s, tt, 0);
                }
            }
            hier_slot_ranks(nxf_quad, quad_layer, nlevels, quad_trank, n_quad_t);
        }
        if (nxvol_hex > 0) {
            vol_layer[0] = (int *)SMESH_ALLOC((size_t)nxvol_hex * sizeof(int));
            vol_trank[0] = (int *)SMESH_ALLOC((size_t)nxvol_hex * sizeof(int));
            hier_fill_vol_layers(HEX8, level, nlevels, levels, vol_layer[0]);
            hier_slot_ranks(nxvol_hex, vol_layer[0], nlevels, vol_trank[0], n_vol_t[0]);
        }
        if (nxvol_tet > 0) {
            vol_layer[1] = (int *)SMESH_ALLOC((size_t)nxvol_tet * sizeof(int));
            vol_trank[1] = (int *)SMESH_ALLOC((size_t)nxvol_tet * sizeof(int));
            hier_fill_vol_layers(TET4, level, nlevels, levels, vol_layer[1]);
            hier_slot_ranks(nxvol_tet, vol_layer[1], nlevels, vol_trank[1], n_vol_t[1]);
        }
        if (nxvol_wedge > 0) {
            vol_layer[2] = (int *)SMESH_ALLOC((size_t)nxvol_wedge * sizeof(int));
            vol_trank[2] = (int *)SMESH_ALLOC((size_t)nxvol_wedge * sizeof(int));
            int t        = 0;
            for (int z = 1; z <= level - 1; ++z) {
                for (int y = 1; y <= level - 2; ++y) {
                    for (int x = 1; x <= level - 1 - y; ++x) {
                        vol_layer[2][t++] = hier_first_layer(level, nlevels, levels, x, y, z);
                    }
                }
            }
            hier_slot_ranks(nxvol_wedge, vol_layer[2], nlevels, vol_trank[2], n_vol_t[2]);
        }
        if (nxvol_pyr > 0) {
            vol_layer[3] = (int *)SMESH_ALLOC((size_t)nxvol_pyr * sizeof(int));
            vol_trank[3] = (int *)SMESH_ALLOC((size_t)nxvol_pyr * sizeof(int));
            int t        = 0;
            for (int k = 1; k <= level - 2; ++k) {
                for (int j = 1; j <= level - k - 1; ++j) {
                    for (int i = 1; i <= level - k - 1; ++i) {
                        vol_layer[3][t++] = hier_first_layer(level, nlevels, levels, i, j, k);
                    }
                }
            }
            hier_slot_ranks(nxvol_pyr, vol_layer[3], nlevels, vol_trank[3], n_vol_t[3]);
        }
        layer_base[0] = 0;
        layer_base[1] = (large_idx_t)n_coarse_global;
        for (int k = 1; k < nlevels; ++k) {
            layer_base[k + 1] = layer_base[k] + (large_idx_t)n_edges_global * (large_idx_t)n_edge_t[k] +
                                (large_idx_t)n_tri_global * (large_idx_t)n_tri_t[k] +
                                (large_idx_t)n_quad_global * (large_idx_t)n_quad_t[k] +
                                (large_idx_t)n_hex_elem_global * (large_idx_t)n_vol_t[0][k] +
                                (large_idx_t)n_tet_elem_global * (large_idx_t)n_vol_t[1][k] +
                                (large_idx_t)n_wedge_elem_global * (large_idx_t)n_vol_t[2][k] +
                                (large_idx_t)n_pyr_elem_global * (large_idx_t)n_vol_t[3][k];
        }
        if ((ptrdiff_t)layer_base[nlevels] != n_ss_global) {
            fprintf(stderr, "to_semistructured: hierarchical HEX-dominant gid count mismatch\n");
            return nullptr;
        }
    }

    int *edge_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *edge_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *tri_uo  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_tri_uniq, 1), sizeof(int));
    int *tri_ua  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_tri_uniq, 1), sizeof(int));
    int *quad_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_quad_uniq, 1), sizeof(int));
    int *quad_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_quad_uniq, 1), sizeof(int));
    int *vol_uo[4], *vol_ua[4];
    for (int s = 0; s < 4; ++s) {
        vol_uo[s] = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_fam[s], 1), sizeof(int));
        vol_ua[s] = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_fam[s], 1), sizeof(int));
    }

    ie = 0;
    itri = 0;
    iq = 0;
    iv[0] = iv[1] = iv[2] = iv[3] = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const enum ElemType f       = block_types[b];
        const int           n_macro = n_macro_of(f);
        const int           nsides  = nsides_of(f);
        const int           fs      = fam_slot(f);
        const ptrdiff_t     n_owned = mesh->block((size_t)b)->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    for (int k = 0; k < nneigh(f, d1); ++k) {
                        if (gc[d1] > gc[neigh(f, d1, k)]) continue;
                        const ptrdiff_t u = edge_inc_to_uniq[ie++];
                        if (from_owned) edge_uo[u] = 1;
                        else edge_ua[u] = 1;
                    }
                }
            }
            for (int s = 0; s < nsides; ++s) {
                if (side_nnxs(f, s) == 3 && nxf_tri > 0) {
                    const ptrdiff_t u = tri_inc_to_uniq[itri++];
                    if (from_owned) tri_uo[u] = 1;
                    else tri_ua[u] = 1;
                } else if (side_nnxs(f, s) == 4 && nxf_quad > 0) {
                    const ptrdiff_t u = quad_inc_to_uniq[iq++];
                    if (from_owned) quad_uo[u] = 1;
                    else quad_ua[u] = 1;
                }
            }
            if (nxvol_of(f) > 0) {
                if (from_owned) vol_uo[fs][iv[fs]] = 1;
                else vol_ua[fs][iv[fs]] = 1;
                iv[fs]++;
            }
        }
    }

    ptrdiff_t n_bkt[4] = {0, 0, 0, 0};
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        if (!c_uo[i] && !c_ua[i]) continue;
        const int sh = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
        n_bkt[node_bucket(rank, coarse_owner[i], sh, c_uo[i], c_ua[i])]++;
    }
    if (nxedge > 0) count_entity_nodes(n_edge_uniq, nxedge, edge_owner, edge_shared, edge_uo, edge_ua, rank, n_bkt);
    if (nxf_tri > 0) count_entity_nodes(n_tri_uniq, nxf_tri, tri_owner, tri_shared, tri_uo, tri_ua, rank, n_bkt);
    if (nxf_quad > 0) count_entity_nodes(n_quad_uniq, nxf_quad, quad_owner, quad_shared, quad_uo, quad_ua, rank, n_bkt);
    for (int s = 0; s < 4; ++s) {
        if (nxvol[s] > 0) count_entity_nodes(n_e_fam[s], nxvol[s], vol_owner[s], vol_shared[s], vol_uo[s], vol_ua[s], rank, n_bkt);
    }

    ptrdiff_t off[5] = {0, 0, 0, 0, 0};
    for (int k = 0; k < 4; ++k) off[k + 1] = off[k] + n_bkt[k];
    const ptrdiff_t n_owned  = off[2];
    const ptrdiff_t n_shared = n_bkt[1];
    const ptrdiff_t n_ghosts = n_bkt[2];
    const ptrdiff_t n_aura   = n_bkt[3];
    const ptrdiff_t n_local  = off[4];

    auto         node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto         node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap         = node_mapping->data();
    int         *nown         = node_owner->data();
    const int    sdim         = mesh->spatial_dimension();
    auto         points       = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_local);
    auto         coarse_p     = mesh->points()->data();
    auto         p            = points->data();
    idx_t *corner_ss = (idx_t *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(idx_t));
    idx_t *edge_ss = (nxedge > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_edge_uniq * (size_t)nxedge * sizeof(idx_t)) : nullptr;
    idx_t *tri_ss  = (nxf_tri > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_tri_uniq * (size_t)nxf_tri * sizeof(idx_t)) : nullptr;
    idx_t *quad_ss = (nxf_quad > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_quad_uniq * (size_t)nxf_quad * sizeof(idx_t)) : nullptr;
    idx_t *vol_ss[4] = {nullptr, nullptr, nullptr, nullptr};
    for (int s = 0; s < 4; ++s) {
        if (nxvol[s] > 0 && n_e_fam[s] > 0) {
            vol_ss[s] = (idx_t *)SMESH_ALLOC((size_t)n_e_fam[s] * (size_t)nxvol[s] * sizeof(idx_t));
        }
    }

    ptrdiff_t cur[4] = {off[0], off[1], off[2], off[3]};
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        if (!c_uo[i] && !c_ua[i]) continue;
        const int sh = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
        const int bk = node_bucket(rank, coarse_owner[i], sh, c_uo[i], c_ua[i]);
        const ptrdiff_t w = cur[bk]++;
        nmap[w] = coarse_nmap[i];
        nown[w] = coarse_owner[i];
        corner_ss[i] = (idx_t)w;
        for (int d = 0; d < sdim; ++d) p[d][w] = coarse_p[d][i];
    }

    large_idx_t *edge_node_gid = alloc_entity_node_gids(n_edge_uniq, nxedge);
    large_idx_t *tri_node_gid  = alloc_entity_node_gids(n_tri_uniq, nxf_tri);
    large_idx_t *quad_node_gid = alloc_entity_node_gids(n_quad_uniq, nxf_quad);
    large_idx_t *vol_node_gid[4];
    for (int s = 0; s < 4; ++s) vol_node_gid[s] = alloc_entity_node_gids(n_e_fam[s], nxvol[s]);
    if (nxedge > 0) {
        if (hierarchical_ordering) {
            fill_hier_hexdom_node_gids(n_edge_uniq, nxedge, 0, edge_gid, edge_layer, edge_trank, n_edge_t, n_tri_t, n_quad_t,
                                       n_vol_t[0], n_vol_t[1], n_vol_t[2], n_vol_t[3], layer_base, n_edges_global, n_tri_global,
                                       n_quad_global, n_hex_elem_global, n_tet_elem_global, n_wedge_elem_global, edge_node_gid);
        } else {
            fill_flat_node_gids(n_edge_uniq, nxedge, edge_base, edge_gid, edge_node_gid);
        }
        pack_entity_nodes(n_edge_uniq, nxedge, edge_node_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    }
    if (nxf_tri > 0) {
        if (hierarchical_ordering) {
            fill_hier_hexdom_node_gids(n_tri_uniq, nxf_tri, 1, tri_gid, tri_layer, tri_trank, n_edge_t, n_tri_t, n_quad_t,
                                       n_vol_t[0], n_vol_t[1], n_vol_t[2], n_vol_t[3], layer_base, n_edges_global, n_tri_global,
                                       n_quad_global, n_hex_elem_global, n_tet_elem_global, n_wedge_elem_global, tri_node_gid);
        } else {
            fill_flat_node_gids(n_tri_uniq, nxf_tri, tri_base, tri_gid, tri_node_gid);
        }
        pack_entity_nodes(n_tri_uniq, nxf_tri, tri_node_gid, tri_owner, tri_shared, tri_uo, tri_ua, rank, cur, nmap, nown, tri_ss);
    }
    if (nxf_quad > 0) {
        if (hierarchical_ordering) {
            fill_hier_hexdom_node_gids(n_quad_uniq, nxf_quad, 2, quad_gid, quad_layer, quad_trank, n_edge_t, n_tri_t, n_quad_t,
                                       n_vol_t[0], n_vol_t[1], n_vol_t[2], n_vol_t[3], layer_base, n_edges_global, n_tri_global,
                                       n_quad_global, n_hex_elem_global, n_tet_elem_global, n_wedge_elem_global, quad_node_gid);
        } else {
            fill_flat_node_gids(n_quad_uniq, nxf_quad, quad_base, quad_gid, quad_node_gid);
        }
        pack_entity_nodes(n_quad_uniq, nxf_quad, quad_node_gid, quad_owner, quad_shared, quad_uo, quad_ua, rank, cur, nmap, nown, quad_ss);
    }
    for (int s = 0; s < 4; ++s) {
        if (nxvol[s] <= 0) continue;
        if (hierarchical_ordering) {
            fill_hier_hexdom_node_gids(n_e_fam[s], nxvol[s], 3 + s, vol_gid[s], vol_layer[s], vol_trank[s], n_edge_t, n_tri_t, n_quad_t,
                                       n_vol_t[0], n_vol_t[1], n_vol_t[2], n_vol_t[3], layer_base, n_edges_global, n_tri_global,
                                       n_quad_global, n_hex_elem_global, n_tet_elem_global, n_wedge_elem_global, vol_node_gid[s]);
        } else {
            fill_flat_node_gids(n_e_fam[s], nxvol[s], vol_base[s], vol_gid[s], vol_node_gid[s]);
        }
        pack_entity_nodes(n_e_fam[s], nxvol[s], vol_node_gid[s], vol_owner[s], vol_shared[s], vol_uo[s], vol_ua[s], rank, cur, nmap, nown, vol_ss[s]);
    }
    SMESH_FREE(edge_node_gid);
    SMESH_FREE(tri_node_gid);
    SMESH_FREE(quad_node_gid);
    for (int s = 0; s < 4; ++s) SMESH_FREE(vol_node_gid[s]);

    int *hex_coords[3], *tet_coords[3], *wedge_coords[3], *pyr_coords[3];
    for (int d = 0; d < 3; ++d) {
        hex_coords[d]   = (int *)SMESH_ALLOC((size_t)sshex8_nxe(level) * sizeof(int));
        tet_coords[d]   = (int *)SMESH_ALLOC((size_t)sstet4_nxe(level) * sizeof(int));
        wedge_coords[d] = (int *)SMESH_ALLOC((size_t)sswedge_nxe(level) * sizeof(int));
        pyr_coords[d]   = (int *)SMESH_ALLOC((size_t)sspyramid_nxe(level) * sizeof(int));
    }
    for (int zi = 0; zi <= level; ++zi)
        for (int yi = 0; yi <= level; ++yi)
            for (int xi = 0; xi <= level; ++xi) {
                const int id = sshex8_lidx(level, xi, yi, zi);
                hex_coords[0][id] = xi; hex_coords[1][id] = yi; hex_coords[2][id] = zi;
            }
    for (int z = 0; z <= level; ++z)
        for (int y = 0; y <= level - z; ++y)
            for (int x = 0; x <= level - z - y; ++x) {
                const int id = sstet4_lidx(level, x, y, z);
                tet_coords[0][id] = x; tet_coords[1][id] = y; tet_coords[2][id] = z;
            }
    for (int z = 0; z <= level; ++z)
        for (int y = 0; y <= level; ++y)
            for (int x = 0; x <= level - y; ++x) {
                const int id = sswedge_lidx(level, x, y, z);
                wedge_coords[0][id] = x; wedge_coords[1][id] = y; wedge_coords[2][id] = z;
            }
    for (int k = 0; k <= level; ++k)
        for (int j = 0; j <= level - k; ++j)
            for (int i = 0; i <= level - k; ++i) {
                const int id = sspyramid_lidx(level, i, j, k);
                pyr_coords[0][id] = i; pyr_coords[1][id] = j; pyr_coords[2][id] = k;
            }
    auto coords_of = [&](const enum ElemType f) -> int ** {
        if (f == HEX8) return hex_coords;
        if (f == TET4) return tet_coords;
        if (f == WEDGE6) return wedge_coords;
        return pyr_coords;
    };

    std::vector<std::shared_ptr<Mesh::Block>> ss_blocks((size_t)n_blocks);
    ie = 0; itri = 0; iq = 0;
    iv[0] = iv[1] = iv[2] = iv[3] = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto                coarse_block = mesh->block((size_t)b);
        const enum ElemType f            = block_types[b];
        const int           n_macro      = n_macro_of(f);
        const int           nsides       = nsides_of(f);
        const int           fs           = fam_slot(f);
        const int           nxe          = nxe_of(f);
        int                *corners      = corners_of(f);
        int               **coords       = coords_of(f);
        ss_lidx_fn          lidx         = lidx_of(f);
        auto                ss_elems     = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_e[b]);
        idx_t             **out          = ss_elems->data();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[8];
            for (int d = 0; d < n_macro; ++d) {
                const idx_t local = coarse_soa[b][d][e];
                gc[d]             = coarse_nmap[local];
                out[corners[d]][e] = corner_ss[local];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    for (int k = 0; k < nneigh(f, d1); ++k) {
                        const int d2 = neigh(f, d1, k);
                        if (gc[d1] > gc[d2]) continue;
                        const idx_t estart = edge_ss[edge_inc_to_uniq[ie++] * nxedge];
                        write_lin_edge(level, corners, coords, lidx, d1, d2, estart, e, out);
                    }
                }
            }
            for (int s = 0; s < nsides; ++s) {
                if (side_nnxs(f, s) == 3 && nxf_tri > 0) {
                    const idx_t foff = tri_ss[tri_inc_to_uniq[itri++] * nxf_tri];
                    write_bary_tri_face(level, gc, lst[fs], corners, coords, lidx, foff, e, s, out);
                } else if (side_nnxs(f, s) == 4 && nxf_quad > 0) {
                    const idx_t foff = quad_ss[quad_inc_to_uniq[iq++] * nxf_quad];
                    write_bilinear_quad_face(level, gc, lst[fs], corners, coords, lidx, foff, e, s, out);
                }
            }
            if (nxvol_of(f) > 0) {
                const idx_t voff = vol_ss[fs][iv[fs] * nxvol_of(f)];
                int         en   = 0;
                if (f == HEX8) {
                    const int Lm1 = level - 1;
                    for (int zi = 1; zi < level; ++zi)
                        for (int yi = 1; yi < level; ++yi)
                            for (int xi = 1; xi < level; ++xi) {
                                en = (zi - 1) * Lm1 * Lm1 + (yi - 1) * Lm1 + (xi - 1);
                                out[sshex8_lidx(level, xi, yi, zi)][e] = voff + (idx_t)en;
                            }
                } else if (f == TET4) {
                    for (int z = 1; z <= level - 3; ++z)
                        for (int y = 1; y <= level - 2 - z; ++y)
                            for (int x = 1; x <= level - 1 - z - y; ++x)
                                out[sstet4_lidx(level, x, y, z)][e] = voff + (idx_t)sstet4_lidx(level - 4, x - 1, y - 1, z - 1);
                } else if (f == WEDGE6) {
                    for (int z = 1; z <= level - 1; ++z)
                        for (int y = 1; y <= level - 2; ++y)
                            for (int x = 1; x <= level - 1 - y; ++x)
                                out[sswedge_lidx(level, x, y, z)][e] = voff + (idx_t)en++;
                } else {
                    for (int k = 1; k <= level - 2; ++k)
                        for (int j = 1; j <= level - k - 1; ++j)
                            for (int i = 1; i <= level - k - 1; ++i)
                                out[sspyramid_lidx(level, i, j, k)][e] = voff + (idx_t)en++;
                }
                iv[fs]++;
            }
        }
        auto ss_block = std::make_shared<Mesh::Block>();
        ss_block->set_name(coarse_block->name());
        ss_block->set_element_type(semistructured_type(f, level));
        ss_block->set_elements(ss_elems);
        ss_block->set_distributed_elements(coarse_block->n_elements_owned(),
                                           coarse_block->n_elements_shared(),
                                           coarse_block->n_elements_ghosts(),
                                           coarse_block->element_mapping(),
                                           coarse_block->aura_element_mapping());
        ss_blocks[(size_t)b] = ss_block;
    }

    auto ghosts_and_aura = create_host_buffer<idx_t>((size_t)(n_ghosts + n_aura));
    auto node_offsets    = create_host_buffer<ptrdiff_t>((size_t)size + 1);
    node_ownership_ranges(comm->get(), n_owned, node_offsets->data());
    SMESH_ASSERT(node_offsets->data()[size] == n_ss_global);
    if (n_ghosts + n_aura > 0) {
        const ptrdiff_t n_id         = rank_split(n_ss_global, size, rank);
        idx_t          *global2owned = (idx_t *)SMESH_CALLOC((size_t)n_id, sizeof(idx_t));
        if (prepare_node_renumbering(comm->get(), n_ss_global, node_offsets->data()[rank], n_owned, nmap, global2owned) !=
            SMESH_SUCCESS) {
            return nullptr;
        }
        if (collect_ghost_and_aura_import_indices(comm->get(), n_owned, n_ghosts, n_aura, n_ss_global, nmap, global2owned,
                                                  node_offsets->data(), ghosts_and_aura->data()) != SMESH_SUCCESS) {
            return nullptr;
        }
        SMESH_FREE(global2owned);
    }

    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        idx_t **el = ss_blocks[(size_t)b]->elements()->data();
        const enum ElemType f = block_types[b];
        if (f == HEX8) sshex8_fill_points(level, n_e[b], el, p, p);
        else if (f == TET4) sstet4_fill_points(level, n_e[b], el, p, p);
        else if (f == WEDGE6) sswedge_fill_points(level, n_e[b], el, p, p);
        else sspyramid_fill_points(level, n_e[b], el, p, p);
    }

    for (int d = 0; d < 3; ++d) {
        SMESH_FREE(hex_coords[d]); SMESH_FREE(tet_coords[d]); SMESH_FREE(wedge_coords[d]); SMESH_FREE(pyr_coords[d]);
    }
    SMESH_FREE(block_types); SMESH_FREE(n_owned_b); SMESH_FREE(n_global_b); SMESH_FREE(block_concat);
    SMESH_FREE(n_e); SMESH_FREE(coarse_soa); SMESH_FREE(c_uo); SMESH_FREE(c_ua);
    SMESH_FREE(edge_inc_to_uniq); SMESH_FREE(edge_gid); SMESH_FREE(edge_owner); SMESH_FREE(edge_shared);
    SMESH_FREE(tri_inc_to_uniq); SMESH_FREE(tri_gid); SMESH_FREE(tri_owner); SMESH_FREE(tri_shared);
    SMESH_FREE(quad_inc_to_uniq); SMESH_FREE(quad_gid); SMESH_FREE(quad_owner); SMESH_FREE(quad_shared);
    for (int s = 0; s < 4; ++s) {
        SMESH_FREE(vol_gid[s]); SMESH_FREE(vol_owner[s]); SMESH_FREE(vol_shared[s]);
        SMESH_FREE(vol_uo[s]); SMESH_FREE(vol_ua[s]); SMESH_FREE(vol_ss[s]);
        SMESH_FREE(vol_layer[s]); SMESH_FREE(vol_trank[s]); SMESH_FREE(n_vol_t[s]);
    }
    SMESH_FREE(edge_uo); SMESH_FREE(edge_ua); SMESH_FREE(tri_uo); SMESH_FREE(tri_ua); SMESH_FREE(quad_uo); SMESH_FREE(quad_ua);
    SMESH_FREE(corner_ss); SMESH_FREE(edge_ss); SMESH_FREE(tri_ss); SMESH_FREE(quad_ss);
    SMESH_FREE(levels); SMESH_FREE(edge_layer); SMESH_FREE(edge_trank);
    SMESH_FREE(tri_layer); SMESH_FREE(tri_trank); SMESH_FREE(quad_layer); SMESH_FREE(quad_trank);
    SMESH_FREE(n_edge_t); SMESH_FREE(n_tri_t); SMESH_FREE(n_quad_t); SMESH_FREE(layer_base);

    auto ret     = std::make_shared<Mesh>(comm, ss_blocks, points);
    auto ss_dist = std::make_shared<Distributed>();
    ss_dist->set_nodes(n_ss_global, n_owned, n_shared, n_ghosts, n_aura, node_mapping, node_owner, node_offsets, ghosts_and_aura);
    ss_dist->set_elements(dist->n_elements_global(), dist->n_elements_owned(), dist->n_elements_shared(),
                          dist->n_elements_ghosts(), dist->element_mapping(), dist->aura_element_mapping());
    ret->set_distributed(ss_dist);
    return ret;
}
