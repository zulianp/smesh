#include "smesh_extractions.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_edgeset.hpp"
#include "smesh_edgesets.hpp"
#include "smesh_extract_shape_features.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"

namespace smesh {

static SharedBuffer<idx_t *> edgeset_node_pairs(Mesh &mesh, const Edgeset &es) {
    auto block = mesh.block(es.block_id());
    if (!block) {
        return nullptr;
    }
    LocalEdgeTable let;
    if (let.fill(block->element_type()) != SMESH_SUCCESS) {
        LocalEdgeTable::report_unsupported("edgeset_node_pairs", block->element_type());
        return nullptr;
    }
    const int       nne = let.nnxe > 0 ? let.nnxe : 2;
    const ptrdiff_t n   = es.size();
    auto            soa = create_host_buffer<idx_t>(nne, n);
    if (n > 0) {
        if (extract_edges_from_edgeset(block->element_type(),
                                       block->elements()->data(),
                                       n,
                                       es.parent()->data(),
                                       es.lei()->data(),
                                       soa->data()) != SMESH_SUCCESS) {
            return nullptr;
        }
    }
    return soa;
}

static int assign_parent_lei_for_node_pairs(Mesh                        &mesh,
                                            const block_idx_t            block_id,
                                            const ptrdiff_t              n_pairs,
                                            const idx_t *const           e0,
                                            const idx_t *const           e1,
                                            SharedBuffer<element_idx_t> &parent_out,
                                            SharedBuffer<i16>           &lei_out) {
    auto block = mesh.block(block_id);
    if (!block) {
        return SMESH_FAILURE;
    }
    const enum ElemType et = block->element_type();
    LocalEdgeTable      let;
    if (let.fill(et) != SMESH_SUCCESS) {
        LocalEdgeTable::report_unsupported("extract_sharp_edges", et);
        return SMESH_FAILURE;
    }

    const int           n_lei  = elem_num_edges(et);
    const ptrdiff_t     n_el   = block->n_elements();
    idx_t *const *const elems  = block->elements()->data();
    auto                n2n    = mesh.node_to_node_graph();
    const count_t      *rowptr = n2n->rowptr()->data();
    const idx_t        *colidx = n2n->colidx()->data();
    const ptrdiff_t     n_nodes = mesh.n_nodes();
    const ptrdiff_t     n_n2n   = (ptrdiff_t)rowptr[n_nodes];

    auto owner_e = create_host_buffer<element_idx_t>((size_t)n_n2n);
    auto owner_s = create_host_buffer<i16>((size_t)n_n2n);
    element_idx_t *d_oe = owner_e->data();
    i16           *d_os = owner_s->data();
    for (ptrdiff_t i = 0; i < n_n2n; ++i) {
        d_oe[i] = invalid_idx<element_idx_t>();
        d_os[i] = -1;
    }

    auto find_slot = [&](idx_t a, idx_t b) -> ptrdiff_t {
        if (a > b) {
            const idx_t t = a;
            a             = b;
            b             = t;
        }
        const ptrdiff_t beg = (ptrdiff_t)rowptr[a];
        const ptrdiff_t end = (ptrdiff_t)rowptr[a + 1];
        for (ptrdiff_t k = beg; k < end; ++k) {
            if (colidx[k] == b) {
                return k;
            }
        }
        return -1;
    };

    for (ptrdiff_t e = 0; e < n_el; ++e) {
        for (int le = 0; le < n_lei; ++le) {
            const idx_t     n0   = elems[let(le, 0)][e];
            const idx_t     n1   = elems[let(le, 1)][e];
            const ptrdiff_t slot = find_slot(n0, n1);
            if (slot < 0) {
                continue;
            }
            if (d_oe[slot] == invalid_idx<element_idx_t>() || (element_idx_t)e < d_oe[slot] ||
                ((element_idx_t)e == d_oe[slot] && (i16)le < d_os[slot])) {
                d_oe[slot] = (element_idx_t)e;
                d_os[slot] = (i16)le;
            }
        }
    }

    parent_out             = create_host_buffer<element_idx_t>((size_t)n_pairs);
    lei_out                = create_host_buffer<i16>((size_t)n_pairs);
    element_idx_t *d_p     = n_pairs > 0 ? parent_out->data() : nullptr;
    i16           *d_l     = n_pairs > 0 ? lei_out->data() : nullptr;
    for (ptrdiff_t i = 0; i < n_pairs; ++i) {
        const ptrdiff_t slot = find_slot(e0[i], e1[i]);
        if (slot < 0 || d_oe[slot] == invalid_idx<element_idx_t>()) {
            SMESH_ERROR("extract_sharp_edges: no parent for node pair (%ld, %ld)\n",
                        (long)e0[i],
                        (long)e1[i]);
            return SMESH_FAILURE;
        }
        d_p[i] = d_oe[slot];
        d_l[i] = d_os[slot];
    }
    return SMESH_SUCCESS;
}

std::shared_ptr<Edgeset> extract_sharp_edges(Mesh &mesh, const geom_t cos_angle_threshold) {
    if (mesh.n_blocks() > 1) {
        SMESH_ERROR("extract_sharp_edges is not supported for multiblock meshes");
        return nullptr;
    }

    const block_idx_t block_id = 0;
    idx_t            *out_e0   = nullptr;
    idx_t            *out_e1   = nullptr;
    ptrdiff_t         out_n    = 0;
    auto              n2n      = mesh.node_to_node_graph();
    extract_sharp_edges(mesh.element_type(block_id),
                        mesh.n_elements(block_id),
                        mesh.elements(block_id)->data(),
                        mesh.n_nodes(),
                        mesh.points()->data(),
                        n2n->rowptr()->data(),
                        n2n->colidx()->data(),
                        cos_angle_threshold,
                        &out_n,
                        &out_e0,
                        &out_e1);

    SharedBuffer<element_idx_t> parent;
    SharedBuffer<i16>           lei;
    const int err = assign_parent_lei_for_node_pairs(mesh, block_id, out_n, out_e0, out_e1, parent, lei);
    SMESH_FREE(out_e0);
    SMESH_FREE(out_e1);
    if (err != SMESH_SUCCESS) {
        return nullptr;
    }

    SharedBuffer<large_idx_t> mapping = nullptr;
    if (mesh.is_distributed()) {
        auto block = mesh.block(block_id);
        if (block) {
            mapping = block->element_mapping();
        }
    }
    return Edgeset::create(mesh.comm(), parent, lei, block_id, mapping);
}

std::shared_ptr<Nodeset> extract_sharp_corners(Mesh                     &mesh,
                                               std::shared_ptr<Edgeset> &sharp_edges,
                                               const bool                edge_clean_up) {
    if (!sharp_edges) {
        return nullptr;
    }
    auto extracted = edgeset_node_pairs(mesh, *sharp_edges);
    const ptrdiff_t n_edges = sharp_edges->size();
    if (n_edges > 0 && !extracted) {
        return nullptr;
    }

    idx_t *e0 = n_edges > 0 ? extracted->data()[0] : nullptr;
    idx_t *e1 = n_edges > 0 ? extracted->data()[1] : nullptr;

    idx_t    *out_corners  = nullptr;
    ptrdiff_t out_ncorners = 0;
    const ptrdiff_t out_n_sharp = extract_sharp_corners(
            mesh.n_nodes(), n_edges, e0, e1, &out_ncorners, &out_corners, edge_clean_up);

    if (edge_clean_up && out_n_sharp != n_edges) {
        SharedBuffer<element_idx_t> parent;
        SharedBuffer<i16>           lei;
        if (assign_parent_lei_for_node_pairs(
                    mesh, sharp_edges->block_id(), out_n_sharp, e0, e1, parent, lei) != SMESH_SUCCESS) {
            SMESH_FREE(out_corners);
            return nullptr;
        }
        sharp_edges = Edgeset::create(mesh.comm(),
                                      parent,
                                      lei,
                                      sharp_edges->block_id(),
                                      sharp_edges->element_mapping());
    }

    SharedBuffer<large_idx_t> mapping = nullptr;
    if (mesh.is_distributed() && mesh.distributed()) {
        mapping = mesh.distributed()->node_mapping();
    }
    return Nodeset::create(mesh.comm(), manage_host_buffer<idx_t>(out_ncorners, out_corners), mapping);
}

SharedBuffer<element_idx_t> extract_disconnected_faces(Mesh &mesh, const Edgeset &sharp_edges) {
    if (mesh.n_blocks() > 1) {
        SMESH_ERROR("extract_disconnected_faces is not supported for multiblock meshes");
        return nullptr;
    }

    auto            extracted = edgeset_node_pairs(mesh, sharp_edges);
    const ptrdiff_t n_edges   = sharp_edges.size();
    if (n_edges > 0 && !extracted) {
        return nullptr;
    }

    const block_idx_t block_id = 0;
    ptrdiff_t         out_n    = 0;
    element_idx_t    *out_el   = nullptr;
    extract_disconnected_faces(mesh.element_type(block_id),
                               mesh.n_elements(block_id),
                               mesh.n_nodes(),
                               mesh.elements(block_id)->data(),
                               n_edges,
                               n_edges > 0 ? extracted->data()[0] : nullptr,
                               n_edges > 0 ? extracted->data()[1] : nullptr,
                               &out_n,
                               &out_el);
    return manage_host_buffer<element_idx_t>(out_n, out_el);
}

}  // namespace smesh
