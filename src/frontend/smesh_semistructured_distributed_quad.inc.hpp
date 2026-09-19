std::shared_ptr<Mesh> to_semistructured_distributed_quad(const int                    level,
                                                         const std::shared_ptr<Mesh> &mesh,
                                                         const bool                   hierarchical_ordering) {
    SMESH_TRACE_SCOPE("to_semistructured_distributed_quad");

    if (level < 1) {
        fprintf(stderr, "to_semistructured: QUAD SS requires level >= 1\n");
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

    const int nxe    = ssquad4_nxe(level);
    const int nxedge = ssquad4_nxedge(level);
    const int nxf    = ssquad4_nxface(level);

    static const int quad_lagr_conn[4][2] = {{1, 3}, {0, 2}, {1, 3}, {0, 2}};

    const ptrdiff_t n_blocks = (ptrdiff_t)mesh->n_blocks();
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        if (mesh->block((size_t)b)->n_nodes_per_element() != 4) {
            fprintf(stderr,
                    "to_semistructured: block '%s' does not have 4 nodes per element\n",
                    mesh->block((size_t)b)->name().c_str());
            return nullptr;
        }
    }

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

    ptrdiff_t           *n_e        = (ptrdiff_t *)SMESH_ALLOC((size_t)n_blocks * sizeof(ptrdiff_t));
    const idx_t *const **coarse_soa = (const idx_t *const **)SMESH_ALLOC((size_t)n_blocks * sizeof(const idx_t *const *));
    ptrdiff_t            n_e_tot    = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto block    = mesh->block((size_t)b);
        n_e[b]        = block->n_elements();
        coarse_soa[b] = block->elements()->data();
        n_e_tot += n_e[b];
    }

    int quad_corners[4];
    quad_corners[0] = ssquad4_lidx(level, 0, 0);
    quad_corners[1] = ssquad4_lidx(level, level, 0);
    quad_corners[2] = ssquad4_lidx(level, level, level);
    quad_corners[3] = ssquad4_lidx(level, 0, level);

    int *c_uo = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));
    int *c_ua = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));

    ptrdiff_t n_edge_inc = 0;
    if (nxedge > 0) {
        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
                large_idx_t gc[4];
                for (int d = 0; d < 4; ++d) {
                    gc[d] = coarse_nmap[coarse_soa[b][d][e]];
                }
                for (int d1 = 0; d1 < 4; ++d1) {
                    for (int k = 0; k < 2; ++k) {
                        if (gc[d1] <= gc[quad_lagr_conn[d1][k]]) {
                            n_edge_inc++;
                        }
                    }
                }
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
    large_idx_t *vol_ids = nullptr;
    large_idx_t *vol_aux = nullptr;
    if (nxf > 0 && n_e_tot > 0) {
        vol_ids = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tot * sizeof(large_idx_t));
        vol_aux = (large_idx_t *)SMESH_ALLOC((size_t)n_e_tot * sizeof(large_idx_t));
    }

    ptrdiff_t ie = 0;
    ptrdiff_t iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block   = mesh->block((size_t)b);
        const ptrdiff_t n_owned = block->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int   from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[4];
            idx_t       lc[4];
            for (int d = 0; d < 4; ++d) {
                lc[d] = coarse_soa[b][d][e];
                gc[d] = coarse_nmap[lc[d]];
                if (from_owned) {
                    c_uo[lc[d]] = 1;
                } else {
                    c_ua[lc[d]] = 1;
                }
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < 4; ++d1) {
                    for (int k = 0; k < 2; ++k) {
                        const int d2 = quad_lagr_conn[d1][k];
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
            if (nxf > 0) {
                vol_ids[iv] = element_gid(*block, concat0[b], e);
                vol_aux[iv] = from_owned ? 0 : 1;
                iv++;
            }
        }
    }
    SMESH_ASSERT(ie == n_edge_inc);
    SMESH_ASSERT(iv == ((nxf > 0) ? n_e_tot : 0));

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

    large_idx_t *vol_gid    = nullptr;
    int         *vol_owner  = nullptr;
    int         *vol_shared = nullptr;
    if (nxf > 0) {
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
    const large_idx_t vol_base  = edge_base + (large_idx_t)n_edges_global * (large_idx_t)nxedge;
    const ptrdiff_t   n_ss_global =
            (ptrdiff_t)vol_base + n_elem_global * (ptrdiff_t)nxf;
    if (n_ss_global < size) {
        fprintf(stderr, "to_semistructured: SS node count smaller than communicator size\n");
        return nullptr;
    }

    const int nlevels = hierarchical_ordering ? sshex8_hierarchical_n_levels(level) : 0;
    int      *levels     = nullptr;
    int      *edge_layer = nullptr;
    int      *edge_trank = nullptr;
    int      *vol_layer  = nullptr;
    int      *vol_trank  = nullptr;
    int      *n_edge_t   = nullptr;
    int      *n_face_t   = nullptr;
    int      *n_vol_t    = nullptr;
    large_idx_t *layer_base = nullptr;
    if (hierarchical_ordering) {
        if (nlevels < 1) {
            fprintf(stderr, "to_semistructured: hierarchical mesh levels cannot be formed\n");
            return nullptr;
        }
        levels     = (int *)SMESH_ALLOC((size_t)nlevels * sizeof(int));
        n_edge_t   = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_face_t   = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        n_vol_t    = (int *)SMESH_CALLOC((size_t)nlevels, sizeof(int));
        layer_base = (large_idx_t *)SMESH_ALLOC((size_t)(nlevels + 1) * sizeof(large_idx_t));
        sshex8_hierarchical_mesh_levels(level, nlevels, levels);
        if (nxedge > 0) {
            edge_layer = (int *)SMESH_ALLOC((size_t)nxedge * sizeof(int));
            edge_trank = (int *)SMESH_ALLOC((size_t)nxedge * sizeof(int));
            hier_fill_edge_layers(level, nxedge, nlevels, levels, edge_layer);
            hier_slot_ranks(nxedge, edge_layer, nlevels, edge_trank, n_edge_t);
        }
        if (nxf > 0) {
            vol_layer = (int *)SMESH_ALLOC((size_t)nxf * sizeof(int));
            vol_trank = (int *)SMESH_ALLOC((size_t)nxf * sizeof(int));
            int t     = 0;
            for (int yi = 1; yi < level; ++yi) {
                for (int xi = 1; xi < level; ++xi) {
                    vol_layer[t++] = hier_first_layer(level, nlevels, levels, xi, yi, 0);
                }
            }
            SMESH_ASSERT(t == nxf);
            hier_slot_ranks(nxf, vol_layer, nlevels, vol_trank, n_vol_t);
        }
        layer_base[0] = 0;
        layer_base[1] = (large_idx_t)n_coarse_global;
        for (int k = 1; k < nlevels; ++k) {
            layer_base[k + 1] = layer_base[k] + (large_idx_t)n_edges_global * (large_idx_t)n_edge_t[k] +
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
    int *vol_uo  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1), sizeof(int));
    int *vol_ua  = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_e_tot, 1), sizeof(int));

    ie = 0;
    iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        const ptrdiff_t n_owned = mesh->block((size_t)b)->n_elements_owned();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            const int   from_owned = e < n_owned ? 1 : 0;
            large_idx_t gc[4];
            for (int d = 0; d < 4; ++d) {
                gc[d] = coarse_nmap[coarse_soa[b][d][e]];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < 4; ++d1) {
                    for (int k = 0; k < 2; ++k) {
                        if (gc[d1] > gc[quad_lagr_conn[d1][k]]) {
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
        count_entity_nodes(n_e_tot, nxf, vol_owner, vol_shared, vol_uo, vol_ua, rank, n_bkt);
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
    idx_t *vol_ss    = (nxf > 0) ? (idx_t *)SMESH_ALLOC((size_t)n_e_tot * (size_t)nxf * sizeof(idx_t)) : nullptr;

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

    large_idx_t *edge_node_gid = alloc_entity_node_gids(n_edge_uniq, nxedge);
    large_idx_t *vol_node_gid  = alloc_entity_node_gids(n_e_tot, nxf);
    if (nxedge > 0) {
        if (hierarchical_ordering) {
            fill_hier_node_gids(n_edge_uniq,
                                nxedge,
                                0,
                                edge_gid,
                                edge_layer,
                                edge_trank,
                                n_edge_t,
                                n_face_t,
                                n_vol_t,
                                layer_base,
                                n_edges_global,
                                0,
                                edge_node_gid);
        } else {
            fill_flat_node_gids(n_edge_uniq, nxedge, edge_base, edge_gid, edge_node_gid);
        }
        pack_entity_nodes(n_edge_uniq, nxedge, edge_node_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    }
    if (nxf > 0) {
        if (hierarchical_ordering) {
            fill_hier_node_gids(n_e_tot,
                                nxf,
                                2,
                                vol_gid,
                                vol_layer,
                                vol_trank,
                                n_edge_t,
                                n_face_t,
                                n_vol_t,
                                layer_base,
                                n_edges_global,
                                0,
                                vol_node_gid);
        } else {
            fill_flat_node_gids(n_e_tot, nxf, vol_base, vol_gid, vol_node_gid);
        }
        pack_entity_nodes(n_e_tot, nxf, vol_node_gid, vol_owner, vol_shared, vol_uo, vol_ua, rank, cur, nmap, nown, vol_ss);
    }
    SMESH_FREE(edge_node_gid);
    SMESH_FREE(vol_node_gid);
    SMESH_FREE(levels);
    SMESH_FREE(edge_layer);
    SMESH_FREE(edge_trank);
    SMESH_FREE(vol_layer);
    SMESH_FREE(vol_trank);
    SMESH_FREE(n_edge_t);
    SMESH_FREE(n_face_t);
    SMESH_FREE(n_vol_t);
    SMESH_FREE(layer_base);

    int *coords[2] = {nullptr, nullptr};
    {
        for (int d = 0; d < 2; ++d) {
            coords[d] = (int *)SMESH_ALLOC((size_t)nxe * sizeof(int));
        }
        for (int yi = 0; yi <= level; ++yi) {
            for (int xi = 0; xi <= level; ++xi) {
                const int lidx    = ssquad4_lidx(level, xi, yi);
                coords[0][lidx] = xi;
                coords[1][lidx] = yi;
            }
        }
    }

    std::vector<std::shared_ptr<Mesh::Block>> ss_blocks((size_t)n_blocks);
    ie = 0;
    iv = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto     coarse_block = mesh->block((size_t)b);
        auto     ss_elems     = create_host_buffer<idx_t>((size_t)nxe, (size_t)n_e[b]);
        idx_t **out           = ss_elems->data();
        for (ptrdiff_t e = 0; e < n_e[b]; ++e) {
            large_idx_t gc[4];
            for (int d = 0; d < 4; ++d) {
                const idx_t local = coarse_soa[b][d][e];
                gc[d]             = coarse_nmap[local];
                out[quad_corners[d]][e] = corner_ss[local];
            }
            if (nxedge > 0) {
                for (int d1 = 0; d1 < 4; ++d1) {
                    for (int k = 0; k < 2; ++k) {
                        const int d2 = quad_lagr_conn[d1][k];
                        if (gc[d1] > gc[d2]) {
                            continue;
                        }
                        const idx_t estart = edge_ss[edge_inc_to_uniq[ie++] * nxedge];
                        const int   lid1   = quad_corners[d1];
                        const int   lid2   = quad_corners[d2];
                        int         P1[2], P2[2];
                        for (int d = 0; d < 2; ++d) {
                            P1[d] = coords[d][lid1];
                            P2[d] = coords[d][lid2];
                        }
                        for (int t = 1; t < level; ++t) {
                            const int xi   = (P1[0] * (level - t) + P2[0] * t) / level;
                            const int yi   = (P1[1] * (level - t) + P2[1] * t) / level;
                            const int lidx = ssquad4_lidx(level, xi, yi);
                            out[lidx][e]   = estart + (idx_t)(t - 1);
                        }
                    }
                }
            }
            if (nxf > 0) {
                const idx_t voff = vol_ss[iv * nxf];
                const int   Lm1  = level - 1;
                for (int yi = 1; yi < level; ++yi) {
                    for (int xi = 1; xi < level; ++xi) {
                        const int lidx = ssquad4_lidx(level, xi, yi);
                        const int en   = (yi - 1) * Lm1 + (xi - 1);
                        out[lidx][e]   = voff + (idx_t)en;
                    }
                }
                iv++;
            }
        }
        auto ss_block = std::make_shared<Mesh::Block>();
        ss_block->set_name(coarse_block->name());
        ss_block->set_element_type(semistructured_type(coarse_block->element_type(), level));
        ss_block->set_elements(ss_elems);
        ss_block->set_distributed_elements(coarse_block->n_elements_owned(),
                                           coarse_block->n_elements_shared(),
                                           coarse_block->n_elements_ghosts(),
                                           coarse_block->element_mapping(),
                                           coarse_block->aura_element_mapping());
        ss_blocks[(size_t)b] = ss_block;
    }
    for (int d = 0; d < 2; ++d) {
        SMESH_FREE(coords[d]);
    }

    auto ghosts_and_aura = create_host_buffer<idx_t>((size_t)(n_ghosts + n_aura));
    auto node_offsets    = create_host_buffer<ptrdiff_t>((size_t)size + 1);
    node_ownership_ranges(comm->get(), n_owned, node_offsets->data());
    SMESH_ASSERT(node_offsets->data()[size] == n_ss_global);

    if (n_ghosts + n_aura > 0) {
        const ptrdiff_t n_id         = rank_split(n_ss_global, size, rank);
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
        ssquad4_fill_points(level, n_e[b], sdim, el, p, p);
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
    SMESH_FREE(vol_gid);
    SMESH_FREE(vol_owner);
    SMESH_FREE(vol_shared);
    SMESH_FREE(edge_uo);
    SMESH_FREE(edge_ua);
    SMESH_FREE(vol_uo);
    SMESH_FREE(vol_ua);
    SMESH_FREE(corner_ss);
    SMESH_FREE(edge_ss);
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

