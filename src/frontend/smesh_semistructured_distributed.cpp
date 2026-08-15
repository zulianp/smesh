#include "smesh_semistructured.hpp"

#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_line_quadrature.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_mesh.hpp"
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

#include "smesh_distributed_unique.inc.hpp"

static int hier_first_layer(const int L, const int nlevels, const int *const levels, const int x, const int y,
                            const int z) {
    for (int k = 0; k < nlevels; ++k) {
        const int step = L / levels[k];
        if ((x % step) == 0 && (y % step) == 0 && (z % step) == 0) {
            return k;
        }
    }
    return nlevels - 1;
}

static void hier_slot_ranks(const int n_slots, const int *const layer, const int nlevels, int *const rank,
                            int *const n_t) {
    for (int k = 0; k < nlevels; ++k) {
        n_t[k] = 0;
    }
    for (int t = 0; t < n_slots; ++t) {
        rank[t] = n_t[layer[t]]++;
    }
}

static void hier_fill_edge_layers(const int L, const int nxedge, const int nlevels, const int *const levels,
                                  int *const layer) {
    for (int t = 0; t < nxedge; ++t) {
        layer[t] = hier_first_layer(L, nlevels, levels, t + 1, 0, 0);
    }
}

static void hier_fill_face_layers(const enum ElemType family,
                                  const int           L,
                                  const int           nxf,
                                  const int           nlevels,
                                  const int          *const levels,
                                  int                *const layer) {
    int t = 0;
    if (family == HEX8) {
        const int Lm1 = L - 1;
        for (int v = 0; v < Lm1; ++v) {
            for (int u = 0; u < Lm1; ++u) {
                layer[t++] = hier_first_layer(L, nlevels, levels, u + 1, v + 1, 0);
            }
        }
    } else {
        for (int tt = 1; tt <= L - 2; ++tt) {
            for (int s = 1; s <= L - 1 - tt; ++s) {
                const int w = L - s - tt;
                layer[t++]  = hier_first_layer(L, nlevels, levels, w, s, tt);
            }
        }
    }
    SMESH_ASSERT(t == nxf);
}

static void hier_fill_vol_layers(const enum ElemType family,
                                 const int           L,
                                 const int           nlevels,
                                 const int          *const levels,
                                 int                *const layer) {
    if (family == HEX8) {
        const int Lm1 = L - 1;
        for (int zi = 1; zi < L; ++zi) {
            for (int yi = 1; yi < L; ++yi) {
                for (int xi = 1; xi < L; ++xi) {
                    const int t = (zi - 1) * Lm1 * Lm1 + (yi - 1) * Lm1 + (xi - 1);
                    layer[t]    = hier_first_layer(L, nlevels, levels, xi, yi, zi);
                }
            }
        }
    } else {
        for (int z = 1; z <= L - 3; ++z) {
            for (int y = 1; y <= L - 2 - z; ++y) {
                for (int x = 1; x <= L - 1 - z - y; ++x) {
                    const int t = sstet4_lidx(L - 4, x - 1, y - 1, z - 1);
                    layer[t]    = hier_first_layer(L, nlevels, levels, x, y, z);
                }
            }
        }
    }
}

static void fill_hier_node_gids(const ptrdiff_t          n_uniq,
                                const int                n_per,
                                const int                kind,
                                const large_idx_t *const entity_gid,
                                const int         *const slot_layer,
                                const int         *const slot_rank,
                                const int         *const n_edge_t,
                                const int         *const n_face_t,
                                const int         *const n_vol_t,
                                const large_idx_t *const layer_base,
                                const ptrdiff_t          n_edges,
                                const ptrdiff_t          n_faces,
                                large_idx_t             *const out) {
    for (ptrdiff_t u = 0; u < n_uniq; ++u) {
        for (int t = 0; t < n_per; ++t) {
            const int         k  = slot_layer[t];
            const int         tr = slot_rank[t];
            large_idx_t       g  = layer_base[k];
            if (kind == 0) {
                g += entity_gid[u] * (large_idx_t)n_edge_t[k] + (large_idx_t)tr;
            } else if (kind == 1) {
                g += (large_idx_t)n_edges * (large_idx_t)n_edge_t[k] + entity_gid[u] * (large_idx_t)n_face_t[k] +
                     (large_idx_t)tr;
            } else {
                g += (large_idx_t)n_edges * (large_idx_t)n_edge_t[k] +
                     (large_idx_t)n_faces * (large_idx_t)n_face_t[k] + entity_gid[u] * (large_idx_t)n_vol_t[k] +
                     (large_idx_t)tr;
            }
            out[u * (ptrdiff_t)n_per + t] = g;
        }
    }
}

/// kind: 0 edge, 1 HEX face, 2 TET face, 3 HEX vol, 4 TET vol
static void fill_hier_mixed_node_gids(const ptrdiff_t          n_uniq,
                                      const int                n_per,
                                      const int                kind,
                                      const large_idx_t *const entity_gid,
                                      const int         *const slot_layer,
                                      const int         *const slot_rank,
                                      const int         *const n_edge_t,
                                      const int         *const n_hex_face_t,
                                      const int         *const n_tet_face_t,
                                      const int         *const n_hex_vol_t,
                                      const int         *const n_tet_vol_t,
                                      const large_idx_t *const layer_base,
                                      const ptrdiff_t          n_edges,
                                      const ptrdiff_t          n_hex_faces,
                                      const ptrdiff_t          n_tet_faces,
                                      const ptrdiff_t          n_hex_elem,
                                      large_idx_t             *const out) {
    for (ptrdiff_t u = 0; u < n_uniq; ++u) {
        for (int t = 0; t < n_per; ++t) {
            const int   k  = slot_layer[t];
            const int   tr = slot_rank[t];
            large_idx_t g  = layer_base[k];
            if (kind == 0) {
                g += entity_gid[u] * (large_idx_t)n_edge_t[k] + (large_idx_t)tr;
            } else if (kind == 1) {
                g += (large_idx_t)n_edges * (large_idx_t)n_edge_t[k] +
                     entity_gid[u] * (large_idx_t)n_hex_face_t[k] + (large_idx_t)tr;
            } else if (kind == 2) {
                g += (large_idx_t)n_edges * (large_idx_t)n_edge_t[k] +
                     (large_idx_t)n_hex_faces * (large_idx_t)n_hex_face_t[k] +
                     entity_gid[u] * (large_idx_t)n_tet_face_t[k] + (large_idx_t)tr;
            } else if (kind == 3) {
                g += (large_idx_t)n_edges * (large_idx_t)n_edge_t[k] +
                     (large_idx_t)n_hex_faces * (large_idx_t)n_hex_face_t[k] +
                     (large_idx_t)n_tet_faces * (large_idx_t)n_tet_face_t[k] +
                     entity_gid[u] * (large_idx_t)n_hex_vol_t[k] + (large_idx_t)tr;
            } else {
                g += (large_idx_t)n_edges * (large_idx_t)n_edge_t[k] +
                     (large_idx_t)n_hex_faces * (large_idx_t)n_hex_face_t[k] +
                     (large_idx_t)n_tet_faces * (large_idx_t)n_tet_face_t[k] +
                     (large_idx_t)n_hex_elem * (large_idx_t)n_hex_vol_t[k] +
                     entity_gid[u] * (large_idx_t)n_tet_vol_t[k] + (large_idx_t)tr;
            }
            out[u * (ptrdiff_t)n_per + t] = g;
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

#include "smesh_semistructured_distributed_hex_tet.inc.hpp"
#include "smesh_semistructured_distributed_quad.inc.hpp"

std::shared_ptr<Mesh> to_semistructured_distributed(const int                    level,
                                                    const std::shared_ptr<Mesh> &mesh,
                                                    const bool                   use_GLL,
                                                    const bool                   hierarchical_ordering) {
    SMESH_TRACE_SCOPE("to_semistructured_distributed");

    bool          has_hex   = false;
    bool          has_tet   = false;
    bool          has_quad  = false;
    bool          has_other = false;
    enum ElemType other_family = INVALID;
    for (size_t b = 0; b < mesh->n_blocks(); ++b) {
        const enum ElemType f = ss_source_family(mesh->element_type(static_cast<block_idx_t>(b)));
        if (f == HEX8) {
            has_hex = true;
        } else if (f == TET4) {
            has_tet = true;
        } else if (f == QUAD4) {
            has_quad = true;
        } else {
            has_other    = true;
            other_family = f;
        }
    }
    if (has_hex && has_tet) {
        if (has_other || has_quad) {
            fprintf(stderr, "to_semistructured: mixed-family semistructured conversion is not implemented\n");
            return nullptr;
        }
        if (use_GLL) {
            fprintf(stderr, "to_semistructured: GLL nodes are not implemented for mixed HEX+TET SS\n");
            return nullptr;
        }
        return to_semistructured_distributed_hex_tet(level, mesh, hierarchical_ordering);
    }
    if ((has_other || has_quad) && (has_hex || has_tet)) {
        fprintf(stderr, "to_semistructured: mixed-family semistructured conversion is not implemented\n");
        return nullptr;
    }
    if (has_quad) {
        if (has_other) {
            fprintf(stderr, "to_semistructured: mixed-family semistructured conversion is not implemented\n");
            return nullptr;
        }
        if (use_GLL) {
            fprintf(stderr, "to_semistructured: GLL nodes are not implemented for QUAD SS\n");
            return nullptr;
        }
        return to_semistructured_distributed_quad(level, mesh, hierarchical_ordering);
    }

    const enum ElemType family = has_hex ? HEX8 : (has_tet ? TET4 : other_family);
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

    const int nlevels = hierarchical_ordering ? sshex8_hierarchical_n_levels(level) : 0;
    int *levels       = nullptr;
    int *edge_layer   = nullptr;
    int *edge_trank   = nullptr;
    int *face_layer   = nullptr;
    int *face_trank   = nullptr;
    int *vol_layer    = nullptr;
    int *vol_trank    = nullptr;
    int *n_edge_t     = nullptr;
    int *n_face_t     = nullptr;
    int *n_vol_t      = nullptr;
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
            face_layer = (int *)SMESH_ALLOC((size_t)nxf * sizeof(int));
            face_trank = (int *)SMESH_ALLOC((size_t)nxf * sizeof(int));
            hier_fill_face_layers(family, level, nxf, nlevels, levels, face_layer);
            hier_slot_ranks(nxf, face_layer, nlevels, face_trank, n_face_t);
        }
        if (nxvol > 0) {
            vol_layer = (int *)SMESH_ALLOC((size_t)nxvol * sizeof(int));
            vol_trank = (int *)SMESH_ALLOC((size_t)nxvol * sizeof(int));
            hier_fill_vol_layers(family, level, nlevels, levels, vol_layer);
            hier_slot_ranks(nxvol, vol_layer, nlevels, vol_trank, n_vol_t);
        }
        layer_base[0] = 0;
        layer_base[1] = (large_idx_t)n_coarse_global;
        for (int k = 1; k < nlevels; ++k) {
            layer_base[k + 1] = layer_base[k] + (large_idx_t)n_edges_global * (large_idx_t)n_edge_t[k] +
                                (large_idx_t)n_faces_global * (large_idx_t)n_face_t[k] +
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

    large_idx_t *edge_node_gid = alloc_entity_node_gids(n_edge_uniq, nxedge);
    large_idx_t *face_node_gid = alloc_entity_node_gids(n_face_uniq, nxf);
    large_idx_t *vol_node_gid  = alloc_entity_node_gids(n_e_tot, nxvol);
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
                                n_faces_global,
                                edge_node_gid);
        } else {
            fill_flat_node_gids(n_edge_uniq, nxedge, edge_base, edge_gid, edge_node_gid);
        }
        pack_entity_nodes(n_edge_uniq, nxedge, edge_node_gid, edge_owner, edge_shared, edge_uo, edge_ua, rank, cur, nmap, nown, edge_ss);
    }
    if (nxf > 0) {
        if (hierarchical_ordering) {
            fill_hier_node_gids(n_face_uniq,
                                nxf,
                                1,
                                face_gid,
                                face_layer,
                                face_trank,
                                n_edge_t,
                                n_face_t,
                                n_vol_t,
                                layer_base,
                                n_edges_global,
                                n_faces_global,
                                face_node_gid);
        } else {
            fill_flat_node_gids(n_face_uniq, nxf, face_base, face_gid, face_node_gid);
        }
        pack_entity_nodes(n_face_uniq, nxf, face_node_gid, face_owner, face_shared, face_uo, face_ua, rank, cur, nmap, nown, face_ss);
    }
    if (nxvol > 0) {
        if (hierarchical_ordering) {
            fill_hier_node_gids(n_e_tot,
                                nxvol,
                                2,
                                vol_gid,
                                vol_layer,
                                vol_trank,
                                n_edge_t,
                                n_face_t,
                                n_vol_t,
                                layer_base,
                                n_edges_global,
                                n_faces_global,
                                vol_node_gid);
        } else {
            fill_flat_node_gids(n_e_tot, nxvol, vol_base, vol_gid, vol_node_gid);
        }
        pack_entity_nodes(n_e_tot, nxvol, vol_node_gid, vol_owner, vol_shared, vol_uo, vol_ua, rank, cur, nmap, nown, vol_ss);
    }
    SMESH_FREE(edge_node_gid);
    SMESH_FREE(face_node_gid);
    SMESH_FREE(vol_node_gid);
    SMESH_FREE(levels);
    SMESH_FREE(edge_layer);
    SMESH_FREE(edge_trank);
    SMESH_FREE(face_layer);
    SMESH_FREE(face_trank);
    SMESH_FREE(vol_layer);
    SMESH_FREE(vol_trank);
    SMESH_FREE(n_edge_t);
    SMESH_FREE(n_face_t);
    SMESH_FREE(n_vol_t);
    SMESH_FREE(layer_base);

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
                                                    const bool /*use_GLL*/,
                                                    const bool /*hierarchical_ordering*/) {
    return nullptr;
}

#endif

}  // namespace smesh

