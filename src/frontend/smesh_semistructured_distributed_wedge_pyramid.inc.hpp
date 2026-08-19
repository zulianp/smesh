std::shared_ptr<Mesh> to_semistructured_distributed_wedge_pyramid(const int                    level,
                                                                  const std::shared_ptr<Mesh> &mesh,
                                                                  const bool                   hierarchical_ordering,
                                                                  const enum ElemType          family) {
    SMESH_TRACE_SCOPE("to_semistructured_distributed_wedge_pyramid");
    if (level < 1) {
        fprintf(stderr, "to_semistructured: WEDGE/PYRAMID SS requires level >= 1\n");
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

    const bool is_wedge = (family == WEDGE6);
    const int  n_macro  = is_wedge ? 6 : 5;
    const int  nsides   = 5;
    const int  nxe      = is_wedge ? sswedge_nxe(level) : sspyramid_nxe(level);
    const int  nxedge   = is_wedge ? sswedge_nxedge(level) : sspyramid_nxedge(level);
    const int  nxf_tri  = is_wedge ? sswedge_nx_tri_face(level) : sspyramid_nx_tri_face(level);
    const int  nxf_quad = is_wedge ? sswedge_nx_quad_face(level) : sspyramid_nx_quad_face(level);
    const int  nxvol    = is_wedge ? sswedge_nxvol(level) : sspyramid_nxvol(level);
    ss_lidx_fn lidx     = is_wedge ? lidx_wedge : lidx_pyr;

    static const int wedge_conn[6][3] = {{1, 2, 3}, {0, 2, 4}, {0, 1, 5}, {0, 4, 5}, {1, 3, 5}, {2, 3, 4}};
    static const int pyr_nneigh[5]    = {3, 3, 3, 3, 4};
    static const int pyr_conn[5][4]   = {{1, 3, 4, -1}, {0, 2, 4, -1}, {1, 3, 4, -1}, {0, 2, 4, -1}, {0, 1, 2, 3}};
    static const int wedge_xyz[6][3]  = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
    static const int pyr_ijk[5][3]    = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}};

    const ptrdiff_t n_blocks = (ptrdiff_t)mesh->n_blocks();
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        if (mesh->block((size_t)b)->n_nodes_per_element() != n_macro) {
            fprintf(stderr,
                    "to_semistructured: block '%s' does not have %d nodes per element\n",
                    mesh->block((size_t)b)->name().c_str(),
                    n_macro);
            return nullptr;
        }
    }

    ptrdiff_t *n_owned_b  = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    ptrdiff_t *n_global_b = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        n_owned_b[b] = mesh->block((size_t)b)->n_elements_owned();
    }
    SMESH_MPI_CATCH(MPI_Allreduce(n_owned_b, n_global_b, (int)n_blocks, mpi_type<ptrdiff_t>(), MPI_SUM, comm->get()));
    large_idx_t *concat0 = (large_idx_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(large_idx_t));
    {
        large_idx_t acc = 0;
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            concat0[b] = acc;
            acc += (large_idx_t)n_global_b[b];
        }
        SMESH_ASSERT((ptrdiff_t)acc == n_elem_global);
    }

    ptrdiff_t           *n_e        = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    const idx_t *const **coarse_soa = (const idx_t *const **)SMESH_ALLOC((size_t)n_blocks * sizeof(const idx_t *const *));
    ptrdiff_t            n_e_tot    = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto block    = mesh->block((size_t)b);
        n_e[b]        = block->n_elements();
        coarse_soa[b] = block->elements()->data();
        n_e_tot += n_e[b];
    }

    LocalSideTable lst;
    lst.fill(is_wedge ? WEDGE6 : PYRAMID5);
    int corners[6];
    if (is_wedge) {
        for (int d = 0; d < 6; ++d) {
            corners[d] = sswedge_lidx(level, wedge_xyz[d][0] * level, wedge_xyz[d][1] * level, wedge_xyz[d][2] * level);
        }
    } else {
        for (int d = 0; d < 5; ++d) {
            corners[d] = (d == 4) ? sspyramid_lidx(level, 0, 0, level)
                                  : sspyramid_lidx(level, pyr_ijk[d][0] * level, pyr_ijk[d][1] * level, pyr_ijk[d][2] * level);
        }
    }

    int *c_uo = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));
    int *c_ua = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));

    auto nneigh = [&](const int d1) { return is_wedge ? 3 : pyr_nneigh[d1]; };
    auto neigh  = [&](const int d1, const int k) { return is_wedge ? wedge_conn[d1][k] : pyr_conn[d1][k]; };
    auto side_nnxs = [&](const int s) {
        if (is_wedge) {
            return s < 3 ? 4 : 3;
        }
        return s < 4 ? 3 : 4;
    };

    ptrdiff_t n_edge_inc = 0, n_tri_inc = 0, n_quad_inc = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[6];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    for (int k = 0; k < nneigh(d1); ++k) {
                        if (gc[d1] <= gc[neigh(d1, k)]) {
                            n_edge_inc++;
                        }
                    }
                }
            }
            for (int s = 0; s < nsides; ++s) {
                if (side_nnxs(s) == 3) {
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
    large_idx_t *vol_ids = nullptr, *vol_aux = nullptr;
    if (nxvol > 0 && n_e_tot > 0) {
        vol_ids = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tot * sizeof(large_idx_t));
        vol_aux = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tot * sizeof(large_idx_t));
    }

    ptrdiff_t ie = 0, it = 0, iq = 0, iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block   = mesh->block((size_t)b);
        const ptrdiff_t n_owned = block->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[6];
            idx_t       lc[6];
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
                    for (int k = 0; k < nneigh(d1); ++k) {
                        const int d2 = neigh(d1, k);
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
                const int nnxs = side_nnxs(s);
                if (nnxs == 3 && nxf_tri > 0) {
                    large_idx_t fk[4] = {k_key_pad, k_key_pad, k_key_pad, k_key_pad};
                    idx_t       loc_min = lc[lst(s, 0)];
                    large_idx_t gmin    = gc[lst(s, 0)];
                    for (int i = 0; i < 3; ++i) {
                        fk[i] = gc[lst(s, i)];
                        if (gc[lst(s, i)] < gmin) {
                            gmin    = gc[lst(s, i)];
                            loc_min = lc[lst(s, i)];
                        }
                    }
                    sort4(fk, 3);
                    tri_keys[it * 4 + 0] = fk[0];
                    tri_keys[it * 4 + 1] = fk[1];
                    tri_keys[it * 4 + 2] = fk[2];
                    tri_keys[it * 4 + 3] = fk[3];
                    tri_aux[it]          = owned_pref_eid_aux(from_owned, element_gid(*block, concat0[b], e), n_elem_global);
                    tri_loc[it]          = loc_min;
                    it++;
                } else if (nnxs == 4 && nxf_quad > 0) {
                    large_idx_t fk[4] = {k_key_pad, k_key_pad, k_key_pad, k_key_pad};
                    idx_t       loc_min = lc[lst(s, 0)];
                    large_idx_t gmin    = gc[lst(s, 0)];
                    for (int i = 0; i < 4; ++i) {
                        fk[i] = gc[lst(s, i)];
                        if (gc[lst(s, i)] < gmin) {
                            gmin    = gc[lst(s, i)];
                            loc_min = lc[lst(s, i)];
                        }
                    }
                    sort4(fk, 4);
                    quad_keys[iq * 4 + 0] = fk[0];
                    quad_keys[iq * 4 + 1] = fk[1];
                    quad_keys[iq * 4 + 2] = fk[2];
                    quad_keys[iq * 4 + 3] = fk[3];
                    quad_aux[iq]          = owned_pref_eid_aux(from_owned, element_gid(*block, concat0[b], e), n_elem_global);
                    quad_loc[iq]          = loc_min;
                    iq++;
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
    SMESH_ASSERT(it == n_tri_inc);
    SMESH_ASSERT(iq == n_quad_inc);

    ptrdiff_t    n_edge_uniq = 0, n_tri_uniq = 0, n_quad_uniq = 0;
    ptrdiff_t   *edge_inc_to_uniq = nullptr, *tri_inc_to_uniq = nullptr, *quad_inc_to_uniq = nullptr;
    large_idx_t *edge_gid = nullptr, *tri_gid = nullptr, *quad_gid = nullptr;
    int         *edge_owner = nullptr, *tri_owner = nullptr, *quad_owner = nullptr;
    int         *edge_shared = nullptr, *tri_shared = nullptr, *quad_shared = nullptr;
    ptrdiff_t    n_edges_global = 0, n_tri_global = 0, n_quad_global = 0;
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
    if (nxf_tri > 0) {
        if (unique_inc_tuples(comm->get(),
                              n_coarse_global,
                              n_tri_inc,
                              tri_keys,
                              tri_aux,
                              tri_loc,
                              n_coarse_local,
                              &n_tri_uniq,
                              &tri_inc_to_uniq,
                              &tri_gid,
                              &tri_owner,
                              &tri_shared,
                              &n_tri_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        tri_keys = nullptr;
        tri_aux  = nullptr;
        tri_loc  = nullptr;
    }
    if (nxf_quad > 0) {
        if (unique_inc_tuples(comm->get(),
                              n_coarse_global,
                              n_quad_inc,
                              quad_keys,
                              quad_aux,
                              quad_loc,
                              n_coarse_local,
                              &n_quad_uniq,
                              &quad_inc_to_uniq,
                              &quad_gid,
                              &quad_owner,
                              &quad_shared,
                              &n_quad_global) != SMESH_SUCCESS) {
            return nullptr;
        }
        quad_keys = nullptr;
        quad_aux  = nullptr;
        quad_loc  = nullptr;
    }
    large_idx_t *vol_gid = nullptr;
    int         *vol_owner = nullptr, *vol_shared = nullptr;
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
    const large_idx_t tri_base  = edge_base + (large_idx_t)n_edges_global * (large_idx_t)nxedge;
    const large_idx_t quad_base = tri_base + (large_idx_t)n_tri_global * (large_idx_t)nxf_tri;
    const large_idx_t vol_base  = quad_base + (large_idx_t)n_quad_global * (large_idx_t)nxf_quad;
    const ptrdiff_t   n_ss_global = (ptrdiff_t)vol_base + n_elem_global * (ptrdiff_t)nxvol;
    if (n_ss_global < size) {
        fprintf(stderr, "to_semistructured: SS node count smaller than communicator size\n");
        return nullptr;
    }

    const int nlevels = hierarchical_ordering ? sshex8_hierarchical_n_levels(level) : 0;
    int *levels = nullptr, *edge_layer = nullptr, *edge_trank = nullptr;
    int *tri_layer = nullptr, *tri_trank = nullptr, *quad_layer = nullptr, *quad_trank = nullptr;
    int *vol_layer = nullptr, *vol_trank = nullptr;
    int *n_edge_t = nullptr, *n_tri_t = nullptr, *n_quad_t = nullptr, *n_vol_t = nullptr;
    large_idx_t *layer_base = nullptr;
    if (hierarchical_ordering) {
        if (nlevels < 1) {
            fprintf(stderr, "to_semistructured: hierarchical mesh levels cannot be formed\n");
            return nullptr;
        }
        levels     = (int *)SMESH_ALLOC((size_t)nlevels * sizeof(int));
        n_edge_t   = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_tri_t    = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_quad_t   = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_vol_t    = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
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
                    const int w   = level - s - tt;
                    tri_layer[t++] = hier_first_layer(level, nlevels, levels, w, s, tt);
                }
            }
            SMESH_ASSERT(t == nxf_tri);
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
            SMESH_ASSERT(t == nxf_quad);
            hier_slot_ranks(nxf_quad, quad_layer, nlevels, quad_trank, n_quad_t);
        }
        if (nxvol > 0) {
            vol_layer = (int *)SMESH_ALLOC((size_t)nxvol * sizeof(int));
            vol_trank = (int *)SMESH_ALLOC((size_t)nxvol * sizeof(int));
            int t     = 0;
            if (is_wedge) {
                for (int z = 1; z <= level - 1; ++z) {
                    for (int y = 1; y <= level - 2; ++y) {
                        for (int x = 1; x <= level - 1 - y; ++x) {
                            vol_layer[t++] = hier_first_layer(level, nlevels, levels, x, y, z);
                        }
                    }
                }
            } else {
                for (int k = 1; k <= level - 2; ++k) {
                    for (int j = 1; j <= level - k - 1; ++j) {
                        for (int i = 1; i <= level - k - 1; ++i) {
                            vol_layer[t++] = hier_first_layer(level, nlevels, levels, i, j, k);
                        }
                    }
                }
            }
            SMESH_ASSERT(t == nxvol);
            hier_slot_ranks(nxvol, vol_layer, nlevels, vol_trank, n_vol_t);
        }
        layer_base[0] = 0;
        layer_base[1] = (large_idx_t)n_coarse_global;
        for (int k = 1; k < nlevels; ++k) {
            layer_base[k + 1] = layer_base[k] + (large_idx_t)n_edges_global * (large_idx_t)n_edge_t[k] +
                                (large_idx_t)n_tri_global * (large_idx_t)n_tri_t[k] +
                                (large_idx_t)n_quad_global * (large_idx_t)n_quad_t[k] +
                                (large_idx_t)n_elem_global * (large_idx_t)n_vol_t[k];
        }
        if ((ptrdiff_t)layer_base[nlevels] != n_ss_global) {
            fprintf(stderr,
                    "to_semistructured: hierarchical SS gid count %ld does not match A8 count %ld\n",
                    (long)layer_base[nlevels],
                    (long)n_ss_global);
            return nullptr;
        }
    }

    int *edge_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *edge_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *tri_uo  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_tri_uniq, 1), sizeof(int));
    int *tri_ua  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_tri_uniq, 1), sizeof(int));
    int *quad_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_quad_uniq, 1), sizeof(int));
    int *quad_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_quad_uniq, 1), sizeof(int));
    int *vol_uo  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1), sizeof(int));
    int *vol_ua  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1), sizeof(int));

    ie = 0;
    it = 0;
    iq = 0;
    iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const ptrdiff_t n_owned = mesh->block((size_t)b)->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int   from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[6];
            for (int d = 0; d < n_macro; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    for (int k = 0; k < nneigh(d1); ++k) {
                        if (gc[d1] > gc[neigh(d1, k)]) {
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
            for (int s = 0; s < nsides; ++s) {
                const int nnxs = side_nnxs(s);
                if (nnxs == 3 && nxf_tri > 0) {
                    const ptrdiff_t u = tri_inc_to_uniq[it++];
                    if (from_owned) {
                        tri_uo[u] = 1;
                    } else {
                        tri_ua[u] = 1;
                    }
                } else if (nnxs == 4 && nxf_quad > 0) {
                    const ptrdiff_t u = quad_inc_to_uniq[iq++];
                    if (from_owned) {
                        quad_uo[u] = 1;
                    } else {
                        quad_ua[u] = 1;
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
    if (nxf_tri > 0) {
        count_entity_nodes(n_tri_uniq, nxf_tri, tri_owner, tri_shared, tri_uo, tri_ua, rank, n_bkt);
    }
    if (nxf_quad > 0) {
        count_entity_nodes(n_quad_uniq, nxf_quad, quad_owner, quad_shared, quad_uo, quad_ua, rank, n_bkt);
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

    auto         node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto         node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap         = node_mapping->data();
    int         *nown         = node_owner->data();
    const int    sdim         = mesh->spatial_dimension();
    auto         points       = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_local);
    auto         coarse_p     = mesh->points()->data();
    auto         p            = points->data();

    idx_t *corner_ss = (idx_t *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(idx_t));
    idx_t *edge_ss   = (nxedge > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_edge_uniq * (size_t)nxedge * sizeof(idx_t)) : nullptr;
    idx_t *tri_ss    = (nxf_tri > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_tri_uniq * (size_t)nxf_tri * sizeof(idx_t)) : nullptr;
    idx_t *quad_ss   = (nxf_quad > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_quad_uniq * (size_t)nxf_quad * sizeof(idx_t)) : nullptr;
    idx_t *vol_ss    = (nxvol > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_e_tot * (size_t)nxvol * sizeof(idx_t)) : nullptr;

    ptrdiff_t cur[4] = {off[0], off[1], off[2], off[3]};
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        if (!c_uo[i] && !c_ua[i]) {
            continue;
        }
        const int       sh = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
        const int       bk = node_bucket(rank, coarse_owner[i], sh, c_uo[i], c_ua[i]);
        const ptrdiff_t w  = cur[bk]++;
        nmap[w]            = coarse_nmap[i];
        nown[w]            = coarse_owner[i];
        corner_ss[i]       = (idx_t)w;
        for (int d = 0; d < sdim; ++d) {
            p[d][w] = coarse_p[d][i];
        }
    }

    large_idx_t *edge_node_gid = alloc_entity_node_gids(n_edge_uniq, nxedge);
    large_idx_t *tri_node_gid  = alloc_entity_node_gids(n_tri_uniq, nxf_tri);
    large_idx_t *quad_node_gid = alloc_entity_node_gids(n_quad_uniq, nxf_quad);
    large_idx_t *vol_node_gid  = alloc_entity_node_gids(n_e_tot, nxvol);
    if (nxedge > 0) {
        if (hierarchical_ordering) {
            fill_hier_prism_node_gids(n_edge_uniq, nxedge, 0, edge_gid, edge_layer, edge_trank, n_edge_t, n_tri_t, n_quad_t, n_vol_t, layer_base, n_edges_global, n_tri_global, n_quad_global, edge_node_gid);
        } else {
            fill_flat_node_gids(n_edge_uniq, nxedge, edge_base, edge_gid, edge_node_gid);
        }
        pack_entity_nodes(n_edge_uniq, nxedge, edge_node_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    }
    if (nxf_tri > 0) {
        if (hierarchical_ordering) {
            fill_hier_prism_node_gids(n_tri_uniq, nxf_tri, 1, tri_gid, tri_layer, tri_trank, n_edge_t, n_tri_t, n_quad_t, n_vol_t, layer_base, n_edges_global, n_tri_global, n_quad_global, tri_node_gid);
        } else {
            fill_flat_node_gids(n_tri_uniq, nxf_tri, tri_base, tri_gid, tri_node_gid);
        }
        pack_entity_nodes(n_tri_uniq, nxf_tri, tri_node_gid, tri_owner, tri_shared, tri_uo, tri_ua, rank, cur, nmap, nown, tri_ss);
    }
    if (nxf_quad > 0) {
        if (hierarchical_ordering) {
            fill_hier_prism_node_gids(n_quad_uniq, nxf_quad, 2, quad_gid, quad_layer, quad_trank, n_edge_t, n_tri_t, n_quad_t, n_vol_t, layer_base, n_edges_global, n_tri_global, n_quad_global, quad_node_gid);
        } else {
            fill_flat_node_gids(n_quad_uniq, nxf_quad, quad_base, quad_gid, quad_node_gid);
        }
        pack_entity_nodes(n_quad_uniq, nxf_quad, quad_node_gid, quad_owner, quad_shared, quad_uo, quad_ua, rank, cur, nmap, nown, quad_ss);
    }
    if (nxvol > 0) {
        if (hierarchical_ordering) {
            fill_hier_prism_node_gids(n_e_tot, nxvol, 3, vol_gid, vol_layer, vol_trank, n_edge_t, n_tri_t, n_quad_t, n_vol_t, layer_base, n_edges_global, n_tri_global, n_quad_global, vol_node_gid);
        } else {
            fill_flat_node_gids(n_e_tot, nxvol, vol_base, vol_gid, vol_node_gid);
        }
        pack_entity_nodes(n_e_tot, nxvol, vol_node_gid, vol_owner, vol_shared, vol_uo, vol_ua, rank, cur, nmap, nown, vol_ss);
    }
    SMESH_FREE(edge_node_gid);
    SMESH_FREE(tri_node_gid);
    SMESH_FREE(quad_node_gid);
    SMESH_FREE(vol_node_gid);
    SMESH_FREE(levels);
    SMESH_FREE(edge_layer);
    SMESH_FREE(edge_trank);
    SMESH_FREE(tri_layer);
    SMESH_FREE(tri_trank);
    SMESH_FREE(quad_layer);
    SMESH_FREE(quad_trank);
    SMESH_FREE(vol_layer);
    SMESH_FREE(vol_trank);
    SMESH_FREE(n_edge_t);
    SMESH_FREE(n_tri_t);
    SMESH_FREE(n_quad_t);
    SMESH_FREE(n_vol_t);
    SMESH_FREE(layer_base);

    int *coords[3] = {nullptr, nullptr, nullptr};
    for (int d = 0; d < 3; ++d) {
        coords[d] = (int *)SMESH_ALLOC((size_t)nxe * sizeof(int));
    }
    if (is_wedge) {
        for (int z = 0; z <= level; ++z) {
            for (int y = 0; y <= level; ++y) {
                for (int x = 0; x <= level - y; ++x) {
                    const int id = sswedge_lidx(level, x, y, z);
                    coords[0][id] = x;
                    coords[1][id] = y;
                    coords[2][id] = z;
                }
            }
        }
    } else {
        for (int k = 0; k <= level; ++k) {
            for (int j = 0; j <= level - k; ++j) {
                for (int i = 0; i <= level - k; ++i) {
                    const int id = sspyramid_lidx(level, i, j, k);
                    coords[0][id] = i;
                    coords[1][id] = j;
                    coords[2][id] = k;
                }
            }
        }
    }

    std::vector<std::shared_ptr<Mesh::Block>> ss_blocks((size_t)n_blocks);
    ie = 0;
    it = 0;
    iq = 0;
    iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto     coarse_block = mesh->block((size_t)b);
        auto     ss_elems     = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_e[b]);
        idx_t **out           = ss_elems->data();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[6];
            for (int d = 0; d < n_macro; ++d) {
                const idx_t local = coarse_soa[b][d][e];
                gc[d]             = coarse_nmap[local];
                out[corners[d]][e] = corner_ss[local];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < n_macro; ++d1) {
                    for (int k = 0; k < nneigh(d1); ++k) {
                        const int d2 = neigh(d1, k);
                        if (gc[d1] > gc[d2]) {
                            continue;
                        }
                        const idx_t estart = edge_ss[edge_inc_to_uniq[ie++] * nxedge];
                        write_lin_edge(level, corners, coords, lidx, d1, d2, estart, e, out);
                    }
                }
            }
            for (int s = 0; s < nsides; ++s) {
                const int nnxs = side_nnxs(s);
                if (nnxs == 3 && nxf_tri > 0) {
                    const idx_t foff = tri_ss[tri_inc_to_uniq[it++] * nxf_tri];
                    write_bary_tri_face(level, gc, lst, corners, coords, lidx, foff, e, s, out);
                } else if (nnxs == 4 && nxf_quad > 0) {
                    const idx_t foff = quad_ss[quad_inc_to_uniq[iq++] * nxf_quad];
                    write_bilinear_quad_face(level, gc, lst, corners, coords, lidx, foff, e, s, out);
                }
            }
            if (nxvol > 0) {
                const idx_t voff = vol_ss[iv * nxvol];
                int         en   = 0;
                if (is_wedge) {
                    for (int z = 1; z <= level - 1; ++z) {
                        for (int y = 1; y <= level - 2; ++y) {
                            for (int x = 1; x <= level - 1 - y; ++x) {
                                out[sswedge_lidx(level, x, y, z)][e] = voff + (idx_t)en++;
                            }
                        }
                    }
                } else {
                    for (int k = 1; k <= level - 2; ++k) {
                        for (int j = 1; j <= level - k - 1; ++j) {
                            for (int i = 1; i <= level - k - 1; ++i) {
                                out[sspyramid_lidx(level, i, j, k)][e] = voff + (idx_t)en++;
                            }
                        }
                    }
                }
                iv++;
            }
        }
        auto ss_block = std::make_shared<Mesh::Block>();
        ss_block->set_name(coarse_block->name());
        ss_block->set_element_type(semistructured_type(family, level));
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
        const ptrdiff_t n_id         = rank_split(n_ss_global, size, rank);
        idx_t          *global2owned = (idx_t *)SMESH_CALLOC((size_t)n_id, sizeof(idx_t));
        if (prepare_node_renumbering(comm->get(), n_ss_global, node_offsets->data()[rank], n_owned, nmap, global2owned) !=
            SMESH_SUCCESS) {
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
        if (is_wedge) {
            sswedge_fill_points(level, n_e[b], el, p, p);
        } else {
            sspyramid_fill_points(level, n_e[b], el, p, p);
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
    SMESH_FREE(tri_inc_to_uniq);
    SMESH_FREE(tri_gid);
    SMESH_FREE(tri_owner);
    SMESH_FREE(tri_shared);
    SMESH_FREE(quad_inc_to_uniq);
    SMESH_FREE(quad_gid);
    SMESH_FREE(quad_owner);
    SMESH_FREE(quad_shared);
    SMESH_FREE(vol_gid);
    SMESH_FREE(vol_owner);
    SMESH_FREE(vol_shared);
    SMESH_FREE(edge_uo);
    SMESH_FREE(edge_ua);
    SMESH_FREE(tri_uo);
    SMESH_FREE(tri_ua);
    SMESH_FREE(quad_uo);
    SMESH_FREE(quad_ua);
    SMESH_FREE(vol_uo);
    SMESH_FREE(vol_ua);
    SMESH_FREE(corner_ss);
    SMESH_FREE(edge_ss);
    SMESH_FREE(tri_ss);
    SMESH_FREE(quad_ss);
    SMESH_FREE(vol_ss);

    auto ret     = std::make_shared<Mesh>(comm, ss_blocks, points);
    auto ss_dist = std::make_shared<Distributed>();
    ss_dist->set_nodes(n_ss_global, n_owned, n_shared, n_ghosts, n_aura, node_mapping, node_owner, node_offsets, ghosts_and_aura);
    ss_dist->set_elements(dist->n_elements_global(),
                          dist->n_elements_owned(),
                          dist->n_elements_shared(),
                          dist->n_elements_ghosts(),
                          dist->element_mapping(),
                          dist->aura_element_mapping());
    ret->set_distributed(ss_dist);
    return ret;
}

std::shared_ptr<Mesh> to_semistructured_distributed_wedge(const int                    level,
                                                          const std::shared_ptr<Mesh> &mesh,
                                                          const bool                   hierarchical_ordering) {
    return to_semistructured_distributed_wedge_pyramid(level, mesh, hierarchical_ordering, WEDGE6);
}

std::shared_ptr<Mesh> to_semistructured_distributed_pyramid(const int                    level,
                                                            const std::shared_ptr<Mesh> &mesh,
                                                            const bool                   hierarchical_ordering) {
    return to_semistructured_distributed_wedge_pyramid(level, mesh, hierarchical_ordering, PYRAMID5);
}
