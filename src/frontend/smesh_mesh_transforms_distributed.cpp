#include "smesh_mesh.hpp"

#ifdef SMESH_ENABLE_MPI

#include "smesh_alloc.hpp"
#include "smesh_alltoallv.impl.hpp"
#include "smesh_common.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

namespace smesh {

namespace {

// Transforms only need edge/tuple helpers; skip SS face/vol unique helpers.
#define SMESH_DISTRIBUTED_UNIQUE_MINIMAL 1
#include "smesh_distributed_unique.inc.hpp"

static const int tet4_refine_pattern[8][4] = {{0, 4, 6, 7},
                                              {4, 1, 5, 8},
                                              {6, 5, 2, 9},
                                              {7, 8, 9, 3},
                                              {4, 5, 6, 8},
                                              {7, 4, 6, 8},
                                              {6, 5, 9, 8},
                                              {7, 6, 9, 8}};

static const int tri3_refine_pattern[4][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};

static const int tet4_edges[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
static const int tri3_edges[3][2] = {{0, 1}, {1, 2}, {0, 2}};

template <typename T>
static SharedBuffer<T> copy_host_buffer(const SharedBuffer<T> &src) {
    if (!src) {
        return nullptr;
    }
    auto dst = create_host_buffer<T>(src->size());
    if (src->size() > 0) {
        std::memcpy(dst->data(), src->data(), src->size() * sizeof(T));
    }
    return dst;
}

static SharedBuffer<large_idx_t> expand_mapping(const SharedBuffer<large_idx_t> &src,
                                                const ptrdiff_t                  n,
                                                const int                        factor) {
    const ptrdiff_t n_out = n * (ptrdiff_t)factor;
    auto            dst   = create_host_buffer<large_idx_t>((size_t)std::max<ptrdiff_t>(n_out, 0));
    if (!src || n <= 0 || factor <= 0) {
        return dst;
    }
    const large_idx_t *s = src->data();
    large_idx_t       *d = dst->data();
    for (ptrdiff_t i = 0; i < n; ++i) {
        for (int c = 0; c < factor; ++c) {
            d[i * (ptrdiff_t)factor + c] = s[i] * (large_idx_t)factor + (large_idx_t)c;
        }
    }
    return dst;
}

static int finalize_ghosts_and_aura(const std::shared_ptr<Communicator> &comm,
                                    const ptrdiff_t                      n_global,
                                    const ptrdiff_t                      n_owned,
                                    const ptrdiff_t                      n_ghosts,
                                    const ptrdiff_t                      n_aura,
                                    const large_idx_t                   *nmap,
                                    SharedBuffer<ptrdiff_t>             &node_offsets,
                                    SharedBuffer<idx_t>                 &ghosts_and_aura) {
    const int size = comm->size();
    const int rank = comm->rank();
    node_offsets    = create_host_buffer<ptrdiff_t>((size_t)size + 1);
    node_ownership_ranges(comm->get(), n_owned, node_offsets->data());
    ghosts_and_aura = create_host_buffer<idx_t>((size_t)(n_ghosts + n_aura));
    if (n_ghosts + n_aura <= 0) {
        return SMESH_SUCCESS;
    }
    if (n_global < (ptrdiff_t)size) {
        fprintf(stderr, "mesh transform: node count smaller than communicator size\n");
        return SMESH_FAILURE;
    }
    const ptrdiff_t n_id        = rank_split(n_global, size, rank);
    idx_t          *global2owned = (idx_t *)SMESH_CALLOC((size_t)n_id, sizeof(idx_t));
    if (prepare_node_renumbering(comm->get(),
                                 n_global,
                                 node_offsets->data()[rank],
                                 n_owned,
                                 nmap,
                                 global2owned) != SMESH_SUCCESS) {
        SMESH_FREE(global2owned);
        return SMESH_FAILURE;
    }
    if (collect_ghost_and_aura_import_indices(comm->get(),
                                              n_owned,
                                              n_ghosts,
                                              n_aura,
                                              n_global,
                                              nmap,
                                              global2owned,
                                              node_offsets->data(),
                                              ghosts_and_aura->data()) != SMESH_SUCCESS) {
        SMESH_FREE(global2owned);
        return SMESH_FAILURE;
    }
    SMESH_FREE(global2owned);
    return SMESH_SUCCESS;
}

static std::shared_ptr<Mesh>
finish_mesh(const std::shared_ptr<Communicator>              &comm,
            const std::vector<std::shared_ptr<Mesh::Block>> &blocks,
            const SharedBuffer<geom_t *>                     &points,
            const ptrdiff_t                                   n_global,
            const ptrdiff_t                                   n_owned,
            const ptrdiff_t                                   n_shared,
            const ptrdiff_t                                   n_ghosts,
            const ptrdiff_t                                   n_aura,
            SharedBuffer<large_idx_t>                         nmap,
            SharedBuffer<int>                                 nown,
            const ptrdiff_t                                   n_elem_global,
            const ptrdiff_t                                   n_elem_owned,
            const ptrdiff_t                                   n_elem_shared,
            const ptrdiff_t                                   n_elem_ghosts,
            SharedBuffer<large_idx_t>                         emap,
            SharedBuffer<large_idx_t>                         aura_emap) {
    auto dist = MeshTransformsDistributed::make_nodal_distributed(comm,
                                                                  n_global,
                                                                  n_owned,
                                                                  n_shared,
                                                                  n_ghosts,
                                                                  n_aura,
                                                                  std::move(nmap),
                                                                  std::move(nown),
                                                                  n_elem_global,
                                                                  n_elem_owned,
                                                                  n_elem_shared,
                                                                  n_elem_ghosts,
                                                                  std::move(emap),
                                                                  std::move(aura_emap));
    if (!dist) {
        return nullptr;
    }
    return MeshTransformsDistributed::make_distributed_mesh(comm, blocks, points, dist);
}

static void expand_block_from(const Mesh::Block &src, Mesh::Block &dst, const int factor) {
    dst.set_distributed_elements(src.n_elements_owned() * (ptrdiff_t)factor,
                                 src.n_elements_shared() * (ptrdiff_t)factor,
                                 src.n_elements_ghosts() * (ptrdiff_t)factor,
                                 expand_mapping(src.element_mapping(), src.n_elements_owned(), factor),
                                 expand_mapping(src.aura_element_mapping(), src.n_elements_ghosts(), factor));
}

static std::shared_ptr<Mesh> refine_edges_once(const std::shared_ptr<Mesh> &mesh) {
    const enum ElemType et            = mesh->element_type(0);
    const int           n_macro       = (et == TET4) ? 4 : 3;
    const int           n_edges       = (et == TET4) ? 6 : 3;
    const int           refine_factor = (et == TET4) ? 8 : 4;

    auto              comm          = mesh->comm();
    const int         rank          = comm->rank();
    const int         size          = comm->size();
    auto              dist          = mesh->distributed();
    const ptrdiff_t   n_coarse_local = mesh->n_nodes();
    const ptrdiff_t   n_coarse_owned = dist->n_nodes_owned();
    const ptrdiff_t   n_coarse_ons   = dist->n_nodes_owned_not_shared();
    const ptrdiff_t   n_coarse_global = dist->n_nodes_global();
    const large_idx_t *coarse_nmap   = dist->node_mapping()->data();
    const int         *coarse_owner  = dist->node_owner()->data();
    const int          sdim          = mesh->spatial_dimension();
    auto               coarse_p      = mesh->points()->data();
    const ptrdiff_t    n_blocks      = (ptrdiff_t)mesh->n_blocks();

    int *c_uo = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));
    int *c_ua = (int *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(int));

    ptrdiff_t n_e_tot   = 0;
    ptrdiff_t n_edge_inc = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block = mesh->block((size_t)b);
        const ptrdiff_t n_e   = block->n_elements();
        n_e_tot += n_e;
        n_edge_inc += n_e * (ptrdiff_t)n_edges;
    }
    SMESH_UNUSED(n_e_tot);

    large_idx_t *edge_keys = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_inc, 1) * 4 * sizeof(large_idx_t));
    large_idx_t *edge_aux  = (large_idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_inc, 1) * sizeof(large_idx_t));
    idx_t       *edge_loc  = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_inc, 1) * sizeof(idx_t));
    idx_t       *edge_a    = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_inc, 1) * sizeof(idx_t));
    idx_t       *edge_b    = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_inc, 1) * sizeof(idx_t));

    ptrdiff_t ie = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block   = mesh->block((size_t)b);
        const ptrdiff_t n_owned = block->n_elements_owned();
        auto            soa     = block->elements()->data();
        for (ptrdiff_t e = 0; e < block->n_elements(); ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            idx_t     lc[4];
            large_idx_t gc[4];
            for (int d = 0; d < n_macro; ++d) {
                lc[d] = soa[d][e];
                gc[d] = coarse_nmap[lc[d]];
                if (from_owned) {
                    c_uo[lc[d]] = 1;
                } else {
                    c_ua[lc[d]] = 1;
                }
            }
            for (int ed = 0; ed < n_edges; ++ed) {
                const int d1 = (et == TET4) ? tet4_edges[ed][0] : tri3_edges[ed][0];
                const int d2 = (et == TET4) ? tet4_edges[ed][1] : tri3_edges[ed][1];
                int       a  = d1;
                int       b2 = d2;
                if (gc[a] > gc[b2]) {
                    const int tmp = a;
                    a             = b2;
                    b2            = tmp;
                }
                edge_keys[ie * 4 + 0] = gc[a];
                edge_keys[ie * 4 + 1] = gc[b2];
                edge_keys[ie * 4 + 2] = k_key_pad;
                edge_keys[ie * 4 + 3] = k_key_pad;
                edge_aux[ie]          = owned_pref_rank_aux(from_owned, rank, size);
                edge_loc[ie]          = lc[a];
                edge_a[ie]            = lc[a];
                edge_b[ie]            = lc[b2];
                ie++;
            }
        }
    }

    ptrdiff_t    n_edge_uniq     = 0;
    ptrdiff_t   *edge_inc_to_uniq = nullptr;
    large_idx_t *edge_gid        = nullptr;
    int         *edge_owner      = nullptr;
    int         *edge_shared     = nullptr;
    ptrdiff_t    n_edges_global  = 0;
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
        SMESH_FREE(c_uo);
        SMESH_FREE(c_ua);
        SMESH_FREE(edge_a);
        SMESH_FREE(edge_b);
        return nullptr;
    }

    int *edge_uo = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    int *edge_ua = (int *)SMESH_CALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1), sizeof(int));
    idx_t *uniq_a = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1) * sizeof(idx_t));
    idx_t *uniq_b = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1) * sizeof(idx_t));
    ie            = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            block   = mesh->block((size_t)b);
        const ptrdiff_t n_owned = block->n_elements_owned();
        for (ptrdiff_t e = 0; e < block->n_elements(); ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            for (int ed = 0; ed < n_edges; ++ed) {
                const ptrdiff_t u = edge_inc_to_uniq[ie];
                uniq_a[u]         = edge_a[ie];
                uniq_b[u]         = edge_b[ie];
                if (from_owned) {
                    edge_uo[u] = 1;
                } else {
                    edge_ua[u] = 1;
                }
                ie++;
            }
        }
    }
    SMESH_FREE(edge_a);
    SMESH_FREE(edge_b);

    ptrdiff_t n_bkt[4] = {0, 0, 0, 0};
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        if (!c_uo[i] && !c_ua[i]) {
            continue;
        }
        const int sh = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
        n_bkt[node_bucket(rank, coarse_owner[i], sh, c_uo[i], c_ua[i])]++;
    }
    count_entity_nodes(n_edge_uniq, 1, edge_owner, edge_shared, edge_uo, edge_ua, rank, n_bkt);

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
        fprintf(stderr, "refine: local node count exceeds idx_t\n");
        return nullptr;
    }

    auto         node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto         node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap         = node_mapping->data();
    int         *nown         = node_owner->data();
    auto         points       = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_local);
    auto         p            = points->data();
    idx_t       *corner_ss    = (idx_t *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(idx_t));
    idx_t       *edge_ss      = (idx_t *)SMESH_ALLOC((size_t)std::max<ptrdiff_t>(n_edge_uniq, 1) * sizeof(idx_t));

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

    large_idx_t *edge_node_gid = alloc_entity_node_gids(n_edge_uniq, 1);
    for (ptrdiff_t u = 0; u < n_edge_uniq; ++u) {
        edge_node_gid[u] = (large_idx_t)n_coarse_global + edge_gid[u];
    }
    pack_entity_nodes(n_edge_uniq, 1, edge_node_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    for (ptrdiff_t u = 0; u < n_edge_uniq; ++u) {
        if (!edge_uo[u] && !edge_ua[u]) {
            continue;
        }
        const idx_t w = edge_ss[u];
        for (int d = 0; d < sdim; ++d) {
            p[d][w] = (geom_t)0.5 * (coarse_p[d][uniq_a[u]] + coarse_p[d][uniq_b[u]]);
        }
    }

    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.reserve((size_t)n_blocks);
    ie = 0;
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        auto            src = mesh->block((size_t)b);
        const ptrdiff_t n_e = src->n_elements();
        auto            soa = src->elements()->data();
        auto            out = create_host_buffer<idx_t>(n_macro, n_e * (ptrdiff_t)refine_factor);
        auto            o   = out->data();
        for (ptrdiff_t e = 0; e < n_e; ++e) {
            idx_t macro_element[10];
            for (int d = 0; d < n_macro; ++d) {
                macro_element[d] = corner_ss[soa[d][e]];
            }
            for (int ed = 0; ed < n_edges; ++ed) {
                macro_element[n_macro + ed] = edge_ss[edge_inc_to_uniq[ie++]];
            }
            const ptrdiff_t element_offset = e * (ptrdiff_t)refine_factor;
            if (et == TET4) {
                for (int k = 0; k < 4; ++k) {
                    for (int sub_e = 0; sub_e < 8; ++sub_e) {
                        o[k][element_offset + sub_e] = macro_element[tet4_refine_pattern[sub_e][k]];
                    }
                }
            } else {
                for (int k = 0; k < 3; ++k) {
                    for (int sub_e = 0; sub_e < 4; ++sub_e) {
                        o[k][element_offset + sub_e] = macro_element[tri3_refine_pattern[sub_e][k]];
                    }
                }
            }
        }
        auto new_block = std::make_shared<Mesh::Block>();
        new_block->set_name(src->name());
        new_block->set_element_type(et);
        new_block->set_elements(out);
        expand_block_from(*src, *new_block, refine_factor);
        blocks.push_back(new_block);
    }

    auto ret = finish_mesh(comm,
                           blocks,
                           points,
                           n_coarse_global + n_edges_global,
                           n_owned,
                           n_shared,
                           n_ghosts,
                           n_aura,
                           node_mapping,
                           node_owner,
                           dist->n_elements_global() * (ptrdiff_t)refine_factor,
                           dist->n_elements_owned() * (ptrdiff_t)refine_factor,
                           dist->n_elements_shared() * (ptrdiff_t)refine_factor,
                           dist->n_elements_ghosts() * (ptrdiff_t)refine_factor,
                           expand_mapping(dist->element_mapping(), dist->n_elements_owned(), refine_factor),
                           expand_mapping(dist->aura_element_mapping(), dist->n_elements_ghosts(), refine_factor));

    SMESH_FREE(c_uo);
    SMESH_FREE(c_ua);
    SMESH_FREE(edge_inc_to_uniq);
    SMESH_FREE(edge_gid);
    SMESH_FREE(edge_owner);
    SMESH_FREE(edge_shared);
    SMESH_FREE(edge_uo);
    SMESH_FREE(edge_ua);
    SMESH_FREE(uniq_a);
    SMESH_FREE(uniq_b);
    SMESH_FREE(corner_ss);
    SMESH_FREE(edge_ss);
    SMESH_FREE(edge_node_gid);
    return ret;
}

}  // namespace

std::shared_ptr<Distributed> MeshTransformsDistributed::make_nodal_distributed(
        const std::shared_ptr<Communicator> &comm,
        ptrdiff_t                            n_global,
        ptrdiff_t                            n_owned,
        ptrdiff_t                            n_shared,
        ptrdiff_t                            n_ghosts,
        ptrdiff_t                            n_aura,
        SharedBuffer<large_idx_t>            node_mapping,
        SharedBuffer<int>                    node_owner,
        ptrdiff_t                            n_elem_global,
        ptrdiff_t                            n_elem_owned,
        ptrdiff_t                            n_elem_shared,
        ptrdiff_t                            n_elem_ghosts,
        SharedBuffer<large_idx_t>            element_mapping,
        SharedBuffer<large_idx_t>            aura_element_mapping) {
    SharedBuffer<ptrdiff_t> node_offsets;
    SharedBuffer<idx_t>     ghosts_and_aura;
    if (finalize_ghosts_and_aura(comm,
                                 n_global,
                                 n_owned,
                                 n_ghosts,
                                 n_aura,
                                 node_mapping->data(),
                                 node_offsets,
                                 ghosts_and_aura) != SMESH_SUCCESS) {
        return nullptr;
    }
    auto dist = std::make_shared<Distributed>();
    dist->set_nodes(n_global,
                    n_owned,
                    n_shared,
                    n_ghosts,
                    n_aura,
                    std::move(node_mapping),
                    std::move(node_owner),
                    std::move(node_offsets),
                    std::move(ghosts_and_aura));
    dist->set_elements(n_elem_global,
                       n_elem_owned,
                       n_elem_shared,
                       n_elem_ghosts,
                       std::move(element_mapping),
                       std::move(aura_element_mapping));
    return dist;
}

std::shared_ptr<Mesh> MeshTransformsDistributed::make_distributed_mesh(
        const std::shared_ptr<Communicator>              &comm,
        const std::vector<std::shared_ptr<Mesh::Block>> &blocks,
        const SharedBuffer<geom_t *>                     &points,
        const std::shared_ptr<Distributed>               &dist) {
    auto ret = std::make_shared<Mesh>(comm, blocks, points);
    ret->set_distributed(dist);
    return ret;
}

std::shared_ptr<Distributed> MeshTransformsDistributed::copy_distributed(const Distributed &src) {
    auto d = std::make_shared<Distributed>();
    d->set_nodes(src.n_nodes_global(),
                 src.n_nodes_owned(),
                 src.n_nodes_shared(),
                 src.n_nodes_ghosts(),
                 src.n_nodes_aura(),
                 copy_host_buffer(src.node_mapping()),
                 copy_host_buffer(src.node_owner()),
                 copy_host_buffer(src.node_offsets()),
                 copy_host_buffer(src.ghosts_and_aura()));
    d->set_elements(src.n_elements_global(),
                    src.n_elements_owned(),
                    src.n_elements_shared(),
                    src.n_elements_ghosts(),
                    copy_host_buffer(src.element_mapping()),
                    copy_host_buffer(src.aura_element_mapping()));
    return d;
}

void MeshTransformsDistributed::expand_element_maps(Distributed &dist, const int factor) {
    dist.set_elements(dist.n_elements_global() * (ptrdiff_t)factor,
                      dist.n_elements_owned() * (ptrdiff_t)factor,
                      dist.n_elements_shared() * (ptrdiff_t)factor,
                      dist.n_elements_ghosts() * (ptrdiff_t)factor,
                      expand_mapping(dist.element_mapping(), dist.n_elements_owned(), factor),
                      expand_mapping(dist.aura_element_mapping(), dist.n_elements_ghosts(), factor));
}

void MeshTransformsDistributed::expand_block_elements(Mesh::Block &block, const int factor) {
    expand_block_from(block, block, factor);
}

int MeshTransformsDistributed::conversion_factor(const enum ElemType from, const enum ElemType to) {
    if (from == to) {
        return 1;
    }
    if (from == HEX8 && to == TET4) {
        return 6;
    }
    if (from == TET15 && to == HEX8) {
        return 4;
    }
    if (from == WEDGE6 && to == TET4) {
        return 3;
    }
    if (from == PYRAMID5 && to == TET4) {
        return 2;
    }
    if (from == QUAD4 && to == TRI3) {
        return 2;
    }
    if (from == HEX8 && to == PROTEUS_HEX8) {
        return 1;
    }
    if (is_hex_ss_family(from) && to == HEX8) {
        return sshex8_txe(semistructured_level(from));
    }
    if (is_quad_ss_family(from) && (to == QUAD4 || to == QUADSHELL4)) {
        const int L = proteus_quad_micro_elements_per_dim(from);
        return L * L;
    }
    return 1;
}

std::shared_ptr<Mesh> MeshTransformsDistributed::refine(const std::shared_ptr<Mesh> &mesh, const int levels) {
    SMESH_TRACE_SCOPE("MeshTransformsDistributed::refine");
    const enum ElemType et = mesh->element_type(0);
    for (size_t b = 1; b < mesh->n_blocks(); ++b) {
        if (mesh->element_type(static_cast<block_idx_t>(b)) != et) {
            SMESH_ERROR("Refinement requires all blocks to share the same element type\n");
            return nullptr;
        }
    }
    if (et == HEX8) {
        auto ss = to_semistructured(1 << levels, mesh);
        if (!ss) {
            return nullptr;
        }
        return sshex_to_hex8(ss);
    }
    if (et != TET4 && et != TRI3) {
        SMESH_ERROR("Refinement is not supported for element type %s\n", type_to_string(et));
        return nullptr;
    }
    auto out = mesh;
    for (int i = 0; i < levels; ++i) {
        out = refine_edges_once(out);
        if (!out) {
            return nullptr;
        }
    }
    return out;
}

std::shared_ptr<Mesh> MeshTransformsDistributed::extrude(const std::shared_ptr<Mesh> &mesh,
                                                         const geom_t                 height,
                                                         const ptrdiff_t              nlayers) {
    SMESH_TRACE_SCOPE("MeshTransformsDistributed::extrude");
    const enum ElemType et = mesh->element_type(0);
    for (size_t b = 1; b < mesh->n_blocks(); ++b) {
        if (mesh->element_type(static_cast<block_idx_t>(b)) != et) {
            SMESH_ERROR("Extrusion requires all blocks to share the same element type\n");
            return nullptr;
        }
    }
    const int n_face = (et == TRI3) ? 3 : 4;
    if (et != QUAD4 && et != QUADSHELL4 && et != TRI3) {
        SMESH_ERROR("Extrusion not supported for element type %s\n", type_to_string(et));
        return nullptr;
    }
    if (nlayers <= 0) {
        SMESH_ERROR("Extrusion requires nlayers > 0\n");
        return nullptr;
    }

    auto              comm            = mesh->comm();
    const int         rank            = comm->rank();
    auto              dist            = mesh->distributed();
    const ptrdiff_t   n_coarse_local  = mesh->n_nodes();
    const ptrdiff_t   n_coarse_owned  = dist->n_nodes_owned();
    const ptrdiff_t   n_coarse_ons    = dist->n_nodes_owned_not_shared();
    const ptrdiff_t   n_coarse_global = dist->n_nodes_global();
    const large_idx_t *coarse_nmap    = dist->node_mapping()->data();
    const int         *coarse_owner   = dist->node_owner()->data();
    auto               coarse_p       = mesh->points()->data();
    const int          n_per          = (int)nlayers + 1;
    const geom_t       dh             = height / (geom_t)nlayers;

    geom_t *nx = (geom_t *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(geom_t));
    geom_t *ny = (geom_t *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(geom_t));
    geom_t *nz = (geom_t *)SMESH_CALLOC((size_t)n_coarse_local, sizeof(geom_t));
    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        auto            block = mesh->block(b);
        auto            soa   = block->elements()->data();
        const ptrdiff_t n_e   = block->n_elements();
        for (ptrdiff_t i = 0; i < n_e; ++i) {
            const idx_t i0 = soa[0][i];
            const idx_t i1 = soa[1][i];
            const idx_t i2 = soa[2][i];
            geom_t      fnx, fny, fnz;
            normal3(coarse_p[0][i0],
                    coarse_p[1][i0],
                    coarse_p[2][i0],
                    coarse_p[0][i1],
                    coarse_p[1][i1],
                    coarse_p[2][i1],
                    coarse_p[0][i2],
                    coarse_p[1][i2],
                    coarse_p[2][i2],
                    &fnx,
                    &fny,
                    &fnz);
            for (int d = 0; d < n_face; ++d) {
                const idx_t n = soa[d][i];
                nx[n] += fnx;
                ny[n] += fny;
                nz[n] += fnz;
            }
        }
    }
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        const geom_t len2 = nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i];
        if (len2 == (geom_t)0) {
            nx[i] = 0;
            ny[i] = 0;
            nz[i] = 1;
        } else {
            normalize3(&nx[i], &ny[i], &nz[i]);
        }
    }

    int *uo = (int *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(int));
    int *ua = (int *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(int));
    int *sh = (int *)SMESH_ALLOC((size_t)n_coarse_local * sizeof(int));
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        uo[i] = i < n_coarse_owned ? 1 : 0;
        ua[i] = i >= n_coarse_owned ? 1 : 0;
        sh[i] = (i >= n_coarse_ons && i < n_coarse_owned) ? 1 : 0;
    }

    ptrdiff_t n_bkt[4] = {0, 0, 0, 0};
    count_entity_nodes(n_coarse_local, n_per, coarse_owner, sh, uo, ua, rank, n_bkt);
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

    auto         node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto         node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap         = node_mapping->data();
    int         *nown         = node_owner->data();
    auto         points       = create_host_buffer<geom_t>(3, (size_t)n_local);
    auto         p            = points->data();
    idx_t       *ss_local     = (idx_t *)SMESH_ALLOC((size_t)n_coarse_local * (size_t)n_per * sizeof(idx_t));
    large_idx_t *layer_gid    = (large_idx_t *)SMESH_ALLOC((size_t)n_coarse_local * (size_t)n_per * sizeof(large_idx_t));
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        for (int l = 0; l < n_per; ++l) {
            layer_gid[i * (ptrdiff_t)n_per + l] =
                    (large_idx_t)l * (large_idx_t)n_coarse_global + coarse_nmap[i];
        }
    }
    ptrdiff_t cur[4] = {off[0], off[1], off[2], off[3]};
    pack_entity_nodes(n_coarse_local, n_per, layer_gid, coarse_owner, sh, uo, ua, rank, cur, nmap, nown, ss_local);
    for (ptrdiff_t i = 0; i < n_coarse_local; ++i) {
        for (int l = 0; l < n_per; ++l) {
            const idx_t w = ss_local[i * (ptrdiff_t)n_per + l];
            p[0][w]       = coarse_p[0][i] + ((geom_t)l * dh * nx[i]);
            p[1][w]       = coarse_p[1][i] + ((geom_t)l * dh * ny[i]);
            p[2][w]       = coarse_p[2][i] + ((geom_t)l * dh * nz[i]);
        }
    }

    const enum ElemType out_type = (et == TRI3) ? WEDGE6 : HEX8;
    const int           nxe_out  = (et == TRI3) ? 6 : 8;
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.reserve(mesh->n_blocks());
    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        auto            src = mesh->block(b);
        const ptrdiff_t n_e = src->n_elements();
        auto            soa = src->elements()->data();
        auto            out = create_host_buffer<idx_t>(nxe_out, n_e * nlayers);
        auto            o   = out->data();
        for (ptrdiff_t e = 0; e < n_e; ++e) {
            for (ptrdiff_t l = 0; l < nlayers; ++l) {
                const ptrdiff_t dst_e = e * nlayers + l;
                for (int d = 0; d < n_face; ++d) {
                    const idx_t node = soa[d][e];
                    o[d][dst_e]           = ss_local[node * (ptrdiff_t)n_per + (int)l];
                    o[n_face + d][dst_e]  = ss_local[node * (ptrdiff_t)n_per + (int)l + 1];
                }
            }
        }
        auto new_block = std::make_shared<Mesh::Block>();
        new_block->set_name(src->name());
        new_block->set_element_type(out_type);
        new_block->set_elements(out);
        expand_block_from(*src, *new_block, (int)nlayers);
        blocks.push_back(new_block);
    }

    auto ret = finish_mesh(comm,
                           blocks,
                           points,
                           n_coarse_global * (ptrdiff_t)n_per,
                           n_owned,
                           n_shared,
                           n_ghosts,
                           n_aura,
                           node_mapping,
                           node_owner,
                           dist->n_elements_global() * nlayers,
                           dist->n_elements_owned() * nlayers,
                           dist->n_elements_shared() * nlayers,
                           dist->n_elements_ghosts() * nlayers,
                           expand_mapping(dist->element_mapping(), dist->n_elements_owned(), (int)nlayers),
                           expand_mapping(dist->aura_element_mapping(), dist->n_elements_ghosts(), (int)nlayers));

    SMESH_FREE(nx);
    SMESH_FREE(ny);
    SMESH_FREE(nz);
    SMESH_FREE(uo);
    SMESH_FREE(ua);
    SMESH_FREE(sh);
    SMESH_FREE(ss_local);
    SMESH_FREE(layer_gid);
    return ret;
}

void MeshTransformsDistributed::clone_distributed(const Mesh &src, Mesh &dst) {
    dst.set_distributed(copy_distributed(*src.distributed()));
    for (size_t b = 0; b < src.n_blocks(); ++b) {
        auto sb = src.block(b);
        auto db = dst.block(b);
        if (!sb || !db || !sb->distributed()) {
            continue;
        }
        db->set_distributed_elements(sb->n_elements_owned(),
                                     sb->n_elements_shared(),
                                     sb->n_elements_ghosts(),
                                     copy_host_buffer(sb->element_mapping()),
                                     copy_host_buffer(sb->aura_element_mapping()));
    }
}

int MeshTransformsDistributed::attach_convert_distributed(const Mesh &src, Mesh &dst) {
    auto dist   = copy_distributed(*src.distributed());
    int  factor = -1;
    for (size_t b = 0; b < src.n_blocks(); ++b) {
        const int f = conversion_factor(src.element_type(static_cast<block_idx_t>(b)),
                                        dst.element_type(static_cast<block_idx_t>(b)));
        expand_block_from(*src.block(b), *dst.block(b), f);
        if (factor < 0) {
            factor = f;
        } else if (factor != f) {
            factor = 0;
        }
    }
    if (factor > 0) {
        expand_element_maps(*dist, factor);
    } else {
        ptrdiff_t n_owned  = 0;
        ptrdiff_t n_shared = 0;
        ptrdiff_t n_ghosts = 0;
        for (size_t b = 0; b < dst.n_blocks(); ++b) {
            auto block = dst.block(b);
            n_owned += block->n_elements_owned();
            n_shared += block->n_elements_shared();
            n_ghosts += block->n_elements_ghosts();
        }
        auto      emap = create_host_buffer<large_idx_t>((size_t)n_owned);
        auto      aura = create_host_buffer<large_idx_t>((size_t)n_ghosts);
        ptrdiff_t io   = 0;
        ptrdiff_t ia   = 0;
        for (size_t b = 0; b < dst.n_blocks(); ++b) {
            auto block = dst.block(b);
            if (block->element_mapping() && block->n_elements_owned() > 0) {
                const ptrdiff_t n = block->n_elements_owned();
                std::memcpy(emap->data() + io, block->element_mapping()->data(), (size_t)n * sizeof(large_idx_t));
                io += n;
            }
            if (block->aura_element_mapping() && block->n_elements_ghosts() > 0) {
                const ptrdiff_t n = block->n_elements_ghosts();
                std::memcpy(aura->data() + ia, block->aura_element_mapping()->data(), (size_t)n * sizeof(large_idx_t));
                ia += n;
            }
        }
        ptrdiff_t n_owned_sum = 0;
        SMESH_MPI_CATCH(MPI_Allreduce(&n_owned, &n_owned_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, src.comm()->get()));
        dist->set_elements(n_owned_sum, n_owned, n_shared, n_ghosts, emap, aura);
    }
    dst.set_distributed(dist);
    return SMESH_SUCCESS;
}

std::shared_ptr<Mesh>
MeshTransformsDistributed::attach_sshex_to_hex8(const std::shared_ptr<Mesh> &ss,
                                                const std::shared_ptr<Mesh> &hex) {
    auto dist = copy_distributed(*ss->distributed());
    int  factor = -1;
    for (size_t b = 0; b < ss->n_blocks(); ++b) {
        const int L = semistructured_level(ss->element_type(static_cast<block_idx_t>(b)));
        const int f = sshex8_txe(L);
        expand_block_from(*ss->block(b), *hex->block(b), f);
        if (factor < 0) {
            factor = f;
        } else if (factor != f) {
            factor = 0;
        }
    }
    if (factor > 0) {
        expand_element_maps(*dist, factor);
    }
    hex->set_distributed(dist);
    return hex;
}

std::shared_ptr<Mesh>
MeshTransformsDistributed::derefine(const std::shared_ptr<Mesh>              &mesh,
                                    std::vector<std::shared_ptr<Mesh::Block>> &blocks) {
    auto            dist           = mesh->distributed();
    const ptrdiff_t n_fine_local   = mesh->n_nodes();
    const ptrdiff_t n_fine_owned   = dist->n_nodes_owned();
    const ptrdiff_t n_fine_ons     = dist->n_nodes_owned_not_shared();
    const large_idx_t *fine_nmap   = dist->node_mapping()->data();
    const int         *fine_owner  = dist->node_owner()->data();
    const int           rank       = mesh->comm()->rank();
    const int           sdim       = mesh->spatial_dimension();
    auto                fine_p     = mesh->points()->data();

    int *used = (int *)SMESH_CALLOC((size_t)n_fine_local, sizeof(int));
    int *uo   = (int *)SMESH_CALLOC((size_t)n_fine_local, sizeof(int));
    int *ua   = (int *)SMESH_CALLOC((size_t)n_fine_local, sizeof(int));
    for (size_t b = 0; b < blocks.size(); ++b) {
        auto            block   = blocks[b];
        auto            src     = mesh->block(b);
        const ptrdiff_t n_owned = src->n_elements_owned();
        auto            soa     = block->elements()->data();
        const int       nxe     = block->n_nodes_per_element();
        for (ptrdiff_t e = 0; e < block->n_elements(); ++e) {
            const int from_owned = e < n_owned ? 1 : 0;
            for (int v = 0; v < nxe; ++v) {
                const idx_t n = soa[v][e];
                used[n]       = 1;
                if (from_owned) {
                    uo[n] = 1;
                } else {
                    ua[n] = 1;
                }
            }
        }
    }

    int *sh = (int *)SMESH_ALLOC((size_t)n_fine_local * sizeof(int));
    for (ptrdiff_t i = 0; i < n_fine_local; ++i) {
        sh[i] = (i >= n_fine_ons && i < n_fine_owned) ? 1 : 0;
        if (!used[i]) {
            uo[i] = 0;
            ua[i] = 0;
        }
    }

    ptrdiff_t n_bkt[4] = {0, 0, 0, 0};
    count_entity_nodes(n_fine_local, 1, fine_owner, sh, uo, ua, rank, n_bkt);
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

    auto         node_mapping = create_host_buffer<large_idx_t>((size_t)n_local);
    auto         node_owner   = create_host_buffer<int>((size_t)n_local);
    large_idx_t *nmap         = node_mapping->data();
    int         *nown         = node_owner->data();
    auto         points       = create_host_buffer<geom_t>((size_t)sdim, (size_t)n_local);
    auto         p            = points->data();
    idx_t       *new_local    = (idx_t *)SMESH_ALLOC((size_t)n_fine_local * sizeof(idx_t));
    ptrdiff_t    cur[4]       = {off[0], off[1], off[2], off[3]};
    pack_entity_nodes(n_fine_local, 1, fine_nmap, fine_owner, sh, uo, ua, rank, cur, nmap, nown, new_local);
    for (ptrdiff_t i = 0; i < n_fine_local; ++i) {
        if (!uo[i] && !ua[i]) {
            continue;
        }
        const idx_t w = new_local[i];
        for (int d = 0; d < sdim; ++d) {
            p[d][w] = fine_p[d][i];
        }
    }

    ptrdiff_t n_owned_sum = 0;
    SMESH_MPI_CATCH(MPI_Allreduce(&n_owned, &n_owned_sum, 1, mpi_type<ptrdiff_t>(), MPI_SUM, mesh->comm()->get()));

    std::vector<std::shared_ptr<Mesh::Block>> out_blocks;
    out_blocks.reserve(blocks.size());
    for (size_t b = 0; b < blocks.size(); ++b) {
        auto            src_view = blocks[b];
        auto            src      = mesh->block(b);
        const int       nxe      = src_view->n_nodes_per_element();
        const ptrdiff_t n_e      = src_view->n_elements();
        auto            in       = src_view->elements()->data();
        auto            elems    = create_host_buffer<idx_t>(nxe, n_e);
        auto            out      = elems->data();
        for (int v = 0; v < nxe; ++v) {
            for (ptrdiff_t e = 0; e < n_e; ++e) {
                out[v][e] = new_local[in[v][e]];
            }
        }
        auto new_block = std::make_shared<Mesh::Block>();
        new_block->set_name(src_view->name());
        new_block->set_element_type(src_view->element_type());
        new_block->set_elements(elems);
        if (src->distributed()) {
            new_block->set_distributed_elements(src->n_elements_owned(),
                                                src->n_elements_shared(),
                                                src->n_elements_ghosts(),
                                                copy_host_buffer(src->element_mapping()),
                                                copy_host_buffer(src->aura_element_mapping()));
        }
        out_blocks.push_back(new_block);
    }

    auto ret = finish_mesh(mesh->comm(),
                           out_blocks,
                           points,
                           n_owned_sum,
                           n_owned,
                           n_shared,
                           n_ghosts,
                           n_aura,
                           node_mapping,
                           node_owner,
                           dist->n_elements_global(),
                           dist->n_elements_owned(),
                           dist->n_elements_shared(),
                           dist->n_elements_ghosts(),
                           copy_host_buffer(dist->element_mapping()),
                           copy_host_buffer(dist->aura_element_mapping()));

    SMESH_FREE(used);
    SMESH_FREE(uo);
    SMESH_FREE(ua);
    SMESH_FREE(sh);
    SMESH_FREE(new_local);
    return ret;
}

}  // namespace smesh

#endif


