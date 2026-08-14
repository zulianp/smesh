#include "smesh_semistructured.hpp"
#include "smesh_alloc.hpp"

// STL
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

#include "smesh_device_buffer.hpp"
#include "smesh_glob.hpp"
#include "smesh_line_quadrature.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_mesh.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sstet4_graph.hpp"
#include "smesh_sstet4_mesh.hpp"
#include "smesh_tracer.hpp"

namespace smesh {

    int semistructured_hierarchical_renumbering(const enum ElemType          element_type,
                                                const int                    level,
                                                const ptrdiff_t              n_nodes,
                                                const SharedBuffer<idx_t *> &elements,
                                                const bool                   preserve_corner_ordering) {
        const int        nlevels = sshex8_hierarchical_n_levels(level);
        std::vector<int> levels(nlevels);
        sshex8_hierarchical_mesh_levels(level, nlevels, levels.data());

        auto node_mapping = create_host_buffer<idx_t>(n_nodes);

        if (ss_source_family(element_type) == TET4) {
            return sstet4_hierarchical_renumbering(level,
                                                   nlevels,
                                                   levels.data(),
                                                   elements->extent(1),
                                                   n_nodes,
                                                   elements->data(),
                                                   node_mapping->data(),
                                                   preserve_corner_ordering);
        }

        return sshex8_hierarchical_renumbering(level,
                                               nlevels,
                                               levels.data(),
                                               elements->extent(1),
                                               n_nodes,
                                               elements->data(),
                                               node_mapping->data(),
                                               preserve_corner_ordering);
    }

    int semistructured_hierarchical_renumbering(const enum ElemType           element_type,
                                                const int                     level,
                                                const ptrdiff_t               n_nodes,
                                                const SharedBuffer<idx_t *>  &elements,
                                                const SharedBuffer<geom_t *> &points,
                                                const bool                    preserve_corner_ordering) {
        const int        nlevels = sshex8_hierarchical_n_levels(level);
        std::vector<int> levels(nlevels);
        sshex8_hierarchical_mesh_levels(level, nlevels, levels.data());

        auto node_mapping = create_host_buffer<idx_t>(n_nodes);

        const int err = (ss_source_family(element_type) == TET4)
                                ? sstet4_hierarchical_renumbering(level,
                                                                  nlevels,
                                                                  levels.data(),
                                                                  elements->extent(1),
                                                                  n_nodes,
                                                                  elements->data(),
                                                                  node_mapping->data(),
                                                                  preserve_corner_ordering)
                                : sshex8_hierarchical_renumbering(level,
                                                                  nlevels,
                                                                  levels.data(),
                                                                  elements->extent(1),
                                                                  n_nodes,
                                                                  elements->data(),
                                                                  node_mapping->data(),
                                                                  preserve_corner_ordering);
        if (err != SMESH_SUCCESS) {
            return err;
        }

        auto      temp = create_host_buffer<geom_t>(n_nodes);
        const int dims = points->extent(0);

        auto        points_data = points->data();
        auto        temp_data   = temp->data();
        const idx_t invalid     = invalid_idx<idx_t>();
        if (!preserve_corner_ordering) {
            for (int d = 0; d < dims; d++) {
                memcpy(temp_data, points_data[d], n_nodes * sizeof(geom_t));
                for (int i = 0; i < n_nodes; i++) {
                    const idx_t idx = node_mapping->data()[i];
                    if (idx != invalid) {
                        SMESH_ASSERT(idx != invalid);
                        points_data[d][idx] = temp_data[i];
                    }
                }
            }
        } else {
            for (int d = 0; d < dims; d++) {
                memcpy(temp_data, points_data[d], n_nodes * sizeof(geom_t));
                for (int i = 0; i < n_nodes; i++) {
                    const idx_t idx = node_mapping->data()[i];

                    SMESH_ASSERT(idx != invalid);
                    points_data[d][idx] = temp_data[i];
                }
            }
        }

        return SMESH_SUCCESS;
    }

    // FIXME hardcoded for sshex8
    std::shared_ptr<Mesh> to_semistructured(const int                    level,
                                            const std::shared_ptr<Mesh> &mesh,
                                            const bool                   hiearchical_ordering,
                                            const bool                   use_GLL) {
        SMESH_TRACE_SCOPE("to_semistructured");

        if (mesh->n_blocks() == 1) {
            auto              block  = mesh->block(0);
            const enum ElemType family = ss_source_family(block->element_type());

            if (family == TET4) {
                if (use_GLL) {
                    fprintf(stderr, "to_semistructured: GLL nodes are not implemented for TET SS\n");
                    return nullptr;
                }
                if (block->n_nodes_per_element() != 4) {
                    fprintf(stderr,
                            "to_semistructured: TET family block '%s' does not have 4 nodes per element\n",
                            block->name().c_str());
                    return nullptr;
                }

                auto default_block = std::make_shared<Mesh::Block>();
                default_block->set_name(block->name());

                enum ElemType element_type = semistructured_type(TET4, level);
                default_block->set_element_type(element_type);

                const int nxe      = sstet4_nxe(level);
                auto      elements = create_host_buffer<idx_t>(nxe, mesh->n_elements());

                ptrdiff_t n_unique_nodes{-1};
                ptrdiff_t interior_start{-1};
                sstet4_generate_elements(level,
                                         mesh->n_elements(),
                                         mesh->n_nodes(),
                                         mesh->elements(0)->data(),
                                         elements->data(),
                                         &n_unique_nodes,
                                         &interior_start);

                default_block->set_elements(elements);
                std::vector<std::shared_ptr<Mesh::Block>> blocks;
                blocks.push_back(default_block);

                if (hiearchical_ordering) {
                    semistructured_hierarchical_renumbering(element_type, level, n_unique_nodes, elements, true);
                }

                auto p       = smesh::create_host_buffer<geom_t>(mesh->spatial_dimension(), n_unique_nodes);
                auto macro_p = mesh->points()->data();
                sstet4_fill_points(level, mesh->n_elements(), elements->data(), macro_p, p->data());
                return std::make_shared<Mesh>(mesh->comm(), blocks, p);
            }

            if (family != HEX8) {
                fprintf(stderr,
                        "to_semistructured: SS family %s is not implemented\n",
                        family == INVALID ? "INVALID" : type_to_string(family));
                return nullptr;
            }

            auto default_block = std::make_shared<Mesh::Block>();
            default_block->set_name(block->name());

            enum ElemType element_type = semistructured_type(block->element_type(), level);
            default_block->set_element_type(element_type);

            const int nxe      = sshex8_nxe(level);
            auto      elements = create_host_buffer<idx_t>(nxe, mesh->n_elements());

            ptrdiff_t n_unique_nodes{-1};
            ptrdiff_t interior_start{-1};
            sshex8_generate_elements(level,
                                     mesh->n_elements(),
                                     mesh->n_nodes(),
                                     mesh->elements(0)->data(),
                                     elements->data(),
                                     &n_unique_nodes,
                                     &interior_start);

            default_block->set_elements(elements);
            std::vector<std::shared_ptr<Mesh::Block>> blocks;
            blocks.push_back(default_block);

            if (hiearchical_ordering) {
                semistructured_hierarchical_renumbering(element_type, level, n_unique_nodes, elements, true);
            }

            auto p       = smesh::create_host_buffer<geom_t>(mesh->spatial_dimension(), n_unique_nodes);
            auto macro_p = mesh->points()->data();

            if (use_GLL) {
                const double *qx{nullptr};
                switch (level) {
                    case 1: {
                        qx = line_GL_q2_x;
                        break;
                    }
                    case 2: {
                        qx = line_GL_q3_x;
                        break;
                    }
                    case 4: {
                        qx = line_GL_q5_x;
                        break;
                    }
                    case 8: {
                        qx = line_GL_q9_x;
                        break;
                    }
                    case 16: {
                        qx = line_GL_q17_x;
                        break;
                    }
                    default: {
                        SMESH_ERROR("Unsupported order %d!", level);
                    }
                }

                sshex8_fill_points_1D_map(level, mesh->n_elements(), elements->data(), macro_p, qx, p->data());
            } else {
                sshex8_fill_points(level, mesh->n_elements(), elements->data(), macro_p, p->data());
            }

            return std::make_shared<Mesh>(mesh->comm(), blocks, p);
        }

        if (mesh->comm()->size() > 1) {
            SMESH_ERROR("to_semistructured is not supported for distributed meshes\n");
            return nullptr;
        }

        enum ElemType family = INVALID;
        for (size_t b = 0; b < mesh->n_blocks(); ++b) {
            const enum ElemType f = ss_source_family(mesh->element_type(static_cast<block_idx_t>(b)));
            if (family == INVALID) {
                family = f;
            } else if (f != family) {
                fprintf(stderr,
                        "to_semistructured: mixed-family semistructured conversion is not implemented\n");
                return nullptr;
            }
        }

        if (family == HEX8) {
        const ptrdiff_t n_blocks = static_cast<ptrdiff_t>(mesh->n_blocks());
        std::vector<ptrdiff_t>                    n_e(static_cast<size_t>(n_blocks));
        std::vector<const idx_t *const *>         hex_soa(static_cast<size_t>(n_blocks));
        std::vector<idx_t **>                     ss_soa(static_cast<size_t>(n_blocks));
        std::vector<std::shared_ptr<Mesh::Block>> ss_blocks(static_cast<size_t>(n_blocks));
        std::vector<const element_idx_t *>        hft(static_cast<size_t>(n_blocks));
        std::vector<const block_idx_t *>          hnbb(static_cast<size_t>(n_blocks));
        std::vector<SharedBuffer<element_idx_t>>  hft_keep(static_cast<size_t>(n_blocks));
        std::vector<SharedBuffer<block_idx_t>>    hnbb_keep(static_cast<size_t>(n_blocks));

        const enum ElemType ss_element_type = semistructured_type(HEX8, level);
        const int           nxe             = sshex8_nxe(level);

        for (ptrdiff_t b = 0; b < n_blocks; ++b) {
            auto block = mesh->block(static_cast<size_t>(b));
            if (block->n_nodes_per_element() != 8) {
                fprintf(stderr,
                        "to_semistructured: HEX family block '%s' does not have 8 nodes per element\n",
                        block->name().c_str());
                return nullptr;
            }
            const size_t bi = static_cast<size_t>(b);
            n_e[bi]         = block->n_elements();
            hex_soa[bi]     = block->elements()->data();

            auto ss_elems = create_host_buffer<idx_t>(nxe, static_cast<size_t>(n_e[bi]));
            ss_soa[bi]    = ss_elems->data();

            auto ss_block = std::make_shared<Mesh::Block>();
            ss_block->set_name(block->name());
            ss_block->set_element_type(ss_element_type);
            ss_block->set_elements(ss_elems);
            ss_blocks[bi] = ss_block;

            hft_keep[bi]  = mesh->half_face_table(static_cast<block_idx_t>(b));
            hnbb_keep[bi] = mesh->half_face_neighbor_block(static_cast<block_idx_t>(b));
            hft[bi]       = hft_keep[bi]->data();
            hnbb[bi]      = hnbb_keep[bi]->data();
        }

        ptrdiff_t n_unique_nodes{-1};
        ptrdiff_t interior_start{-1};
        sshex8_generate_elements_blocks(level,
                                        n_blocks,
                                        n_e.data(),
                                        mesh->n_nodes(),
                                        hex_soa.data(),
                                        ss_soa.data(),
                                        hft.data(),
                                        hnbb.data(),
                                        &n_unique_nodes,
                                        &interior_start);

        if (hiearchical_ordering) {
            const int        nlevels = sshex8_hierarchical_n_levels(level);
            std::vector<int> levels(static_cast<size_t>(nlevels));
            sshex8_hierarchical_mesh_levels(level, nlevels, levels.data());
            auto node_mapping = create_host_buffer<idx_t>(n_unique_nodes);
            sshex8_hierarchical_renumbering_blocks(level,
                                                   nlevels,
                                                   levels.data(),
                                                   n_blocks,
                                                   n_e.data(),
                                                   n_unique_nodes,
                                                   ss_soa.data(),
                                                   node_mapping->data(),
                                                   true);
        }

        auto p       = smesh::create_host_buffer<geom_t>(mesh->spatial_dimension(), n_unique_nodes);
        auto macro_p = mesh->points()->data();

        if (use_GLL) {
            const double *qx{nullptr};
            switch (level) {
                case 1: {
                    qx = line_GL_q2_x;
                    break;
                }
                case 2: {
                    qx = line_GL_q3_x;
                    break;
                }
                case 4: {
                    qx = line_GL_q5_x;
                    break;
                }
                case 8: {
                    qx = line_GL_q9_x;
                    break;
                }
                case 16: {
                    qx = line_GL_q17_x;
                    break;
                }
                default: {
                    SMESH_ERROR("Unsupported order %d!", level);
                }
            }

            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                sshex8_fill_points_1D_map(level, n_e[static_cast<size_t>(b)], ss_soa[static_cast<size_t>(b)], macro_p, qx, p->data());
            }
        } else {
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                sshex8_fill_points(level, n_e[static_cast<size_t>(b)], ss_soa[static_cast<size_t>(b)], macro_p, p->data());
            }
        }

        return std::make_shared<Mesh>(mesh->comm(), ss_blocks, p);
        }

        if (family == TET4) {
            if (use_GLL) {
                fprintf(stderr, "to_semistructured: GLL nodes are not implemented for TET SS\n");
                return nullptr;
            }

            const ptrdiff_t n_blocks = static_cast<ptrdiff_t>(mesh->n_blocks());
            std::vector<ptrdiff_t>                    n_e(static_cast<size_t>(n_blocks));
            std::vector<const idx_t *const *>         tet_soa(static_cast<size_t>(n_blocks));
            std::vector<idx_t **>                     ss_soa(static_cast<size_t>(n_blocks));
            std::vector<std::shared_ptr<Mesh::Block>> ss_blocks(static_cast<size_t>(n_blocks));
            std::vector<const element_idx_t *>        hft(static_cast<size_t>(n_blocks));
            std::vector<const block_idx_t *>          hnbb(static_cast<size_t>(n_blocks));
            std::vector<SharedBuffer<element_idx_t>>  hft_keep(static_cast<size_t>(n_blocks));
            std::vector<SharedBuffer<block_idx_t>>    hnbb_keep(static_cast<size_t>(n_blocks));

            const enum ElemType ss_element_type = semistructured_type(TET4, level);
            const int           nxe             = sstet4_nxe(level);

            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                auto block = mesh->block(static_cast<size_t>(b));
                if (block->n_nodes_per_element() != 4) {
                    fprintf(stderr,
                            "to_semistructured: TET family block '%s' does not have 4 nodes per element\n",
                            block->name().c_str());
                    return nullptr;
                }
                const size_t bi = static_cast<size_t>(b);
                n_e[bi]         = block->n_elements();
                tet_soa[bi]     = block->elements()->data();

                auto ss_elems = create_host_buffer<idx_t>(nxe, static_cast<size_t>(n_e[bi]));
                ss_soa[bi]    = ss_elems->data();

                auto ss_block = std::make_shared<Mesh::Block>();
                ss_block->set_name(block->name());
                ss_block->set_element_type(ss_element_type);
                ss_block->set_elements(ss_elems);
                ss_blocks[bi] = ss_block;

                hft_keep[bi]  = mesh->half_face_table(static_cast<block_idx_t>(b));
                hnbb_keep[bi] = mesh->half_face_neighbor_block(static_cast<block_idx_t>(b));
                hft[bi]       = hft_keep[bi]->data();
                hnbb[bi]      = hnbb_keep[bi]->data();
            }

            ptrdiff_t n_unique_nodes{-1};
            ptrdiff_t interior_start{-1};
            sstet4_generate_elements_blocks(level,
                                            n_blocks,
                                            n_e.data(),
                                            mesh->n_nodes(),
                                            tet_soa.data(),
                                            ss_soa.data(),
                                            hft.data(),
                                            hnbb.data(),
                                            &n_unique_nodes,
                                            &interior_start);

            if (hiearchical_ordering) {
                const int        nlevels = sshex8_hierarchical_n_levels(level);
                std::vector<int> levels(static_cast<size_t>(nlevels));
                sshex8_hierarchical_mesh_levels(level, nlevels, levels.data());
                auto node_mapping = create_host_buffer<idx_t>(n_unique_nodes);
                sstet4_hierarchical_renumbering_blocks(level,
                                                       nlevels,
                                                       levels.data(),
                                                       n_blocks,
                                                       n_e.data(),
                                                       n_unique_nodes,
                                                       ss_soa.data(),
                                                       node_mapping->data(),
                                                       true);
            }

            auto p       = smesh::create_host_buffer<geom_t>(mesh->spatial_dimension(), n_unique_nodes);
            auto macro_p = mesh->points()->data();
            for (ptrdiff_t b = 0; b < n_blocks; ++b) {
                sstet4_fill_points(level, n_e[static_cast<size_t>(b)], ss_soa[static_cast<size_t>(b)], macro_p, p->data());
            }

            return std::make_shared<Mesh>(mesh->comm(), ss_blocks, p);
        }

        fprintf(stderr,
                "to_semistructured: SS family %s is not implemented\n",
                family == INVALID ? "INVALID" : type_to_string(family));
        return nullptr;
    }

    void sshex_block_to_hex8_block(const Mesh::Block &block, Mesh::Block &new_block) {
        auto      elements         = block.elements();
        const int level            = semistructured_level(block.element_type());
        ptrdiff_t n_micro_elements = block.n_elements() * level * level * level;
        auto      hex8_elements    = create_host_buffer<idx_t>(8, n_micro_elements);
        sshex8_to_standard_hex8_mesh(level, block.n_elements(), elements->data(), hex8_elements->data());

        new_block.set_name(block.name());
        new_block.set_elements(hex8_elements);
        new_block.set_element_type(HEX8);
    }

    void ssquad_block_to_quad4_block(const Mesh::Block &block, Mesh::Block &new_block) {
        const int level   = proteus_quad_micro_elements_per_dim(block.element_type());
        const int nnxs    = 4;
        const int nexs    = level * level;
        auto      surface = smesh::create_host_buffer<idx_t>(nnxs, block.n_elements() * nexs);

        ssquad4_to_standard_quad4_mesh(level, block.n_elements(), block.elements()->data(), surface->data());

        new_block.set_name(block.name());
        new_block.set_elements(surface);
        new_block.set_element_type(QUAD4);
    }

    std::shared_ptr<Mesh> sshex_to_hex8(const std::shared_ptr<Mesh> &sshex) {
        std::vector<std::shared_ptr<Mesh::Block>> blocks;
        for (auto &block : sshex->blocks()) {
            auto new_block = std::make_shared<Mesh::Block>();
            sshex_block_to_hex8_block(*block, *new_block);
            blocks.push_back(new_block);
        }

        return std::make_shared<Mesh>(sshex->comm(), blocks, sshex->points());
    }

    std::shared_ptr<Mesh> ssquad_to_quad4(const std::shared_ptr<Mesh> &ssquad) {
        std::vector<std::shared_ptr<Mesh::Block>> blocks;
        for (auto &block : ssquad->blocks()) {
            auto new_block = std::make_shared<Mesh::Block>();
            ssquad_block_to_quad4_block(*block, *new_block);
            blocks.push_back(new_block);
        }

        auto ret = std::make_shared<Mesh>(ssquad->comm(), blocks, ssquad->points());
        ret->set_node_mapping(ssquad->node_mapping());
        return ret;
    }

    std::shared_ptr<Mesh> derefine(const std::shared_ptr<Mesh> &mesh, const int to_level) {
        if (mesh->n_blocks() == 0) {
            fprintf(stderr, "derefine: empty mesh\n");
            return nullptr;
        }

        const enum ElemType family = ss_source_family(mesh->element_type(0));
        for (size_t b = 1; b < mesh->n_blocks(); ++b) {
            if (ss_source_family(mesh->element_type(static_cast<block_idx_t>(b))) != family) {
                fprintf(stderr, "derefine: mixed-family meshes are not implemented\n");
                return nullptr;
            }
        }

        if (family == TET4) {
            std::vector<std::shared_ptr<Mesh::Block>> blocks;
            ptrdiff_t                                 n_unique_nodes{-1};
            for (auto &block : mesh->blocks()) {
                if (!is_semistructured_type(block->element_type()) || !is_tet_ss_family(block->element_type())) {
                    fprintf(stderr, "derefine: only TET-family semistructured blocks are implemented\n");
                    return nullptr;
                }

                const int from_level  = semistructured_level(block->element_type());
                if (to_level <= 0 || from_level < to_level || (from_level % to_level) != 0) {
                    fprintf(stderr, "derefine: invalid levels from=%d to=%d\n", from_level, to_level);
                    return nullptr;
                }
                const int step_factor = from_level / to_level;
                const int nxe         = sstet4_nxe(to_level);

                auto elements = block->elements();
                auto view     = std::make_shared<Buffer<idx_t *>>(
                        nxe,
                        block->n_elements(),
                        (idx_t **)SMESH_ALLOC(nxe * sizeof(idx_t *)),
                        [keep_alive = elements](int, void **v) {
                            (void)keep_alive;
                            SMESH_FREE(v);
                        },
                        elements->mem_space());

                for (int z = 0; z <= to_level; ++z) {
                    for (int y = 0; y <= to_level - z; ++y) {
                        for (int x = 0; x <= to_level - z - y; ++x) {
                            const int from_lidx   = sstet4_lidx(from_level, x * step_factor, y * step_factor, z * step_factor);
                            const int to_lidx     = sstet4_lidx(to_level, x, y, z);
                            view->data()[to_lidx] = elements->data()[from_lidx];
                        }
                    }
                }

                auto derefined_block = std::make_shared<Mesh::Block>();
                derefined_block->set_name(block->name());
                derefined_block->set_element_type(semistructured_type(TET4, to_level));
                derefined_block->set_elements(view);
                blocks.push_back(derefined_block);

                auto            vv        = view->data();
                const ptrdiff_t nelements = block->n_elements();
                for (size_t v = 0; v < view->extent(0); v++) {
                    for (ptrdiff_t e = 0; e < nelements; e++) {
                        n_unique_nodes = std::max(static_cast<ptrdiff_t>(vv[v][e]), n_unique_nodes);
                    }
                }
            }

            n_unique_nodes += 1;
            int  sdim   = mesh->spatial_dimension();
            auto points = smesh::view(mesh->points(), 0, sdim, 0, n_unique_nodes);
            return std::make_shared<Mesh>(mesh->comm(), blocks, points);
        }

        std::vector<std::shared_ptr<Mesh::Block>> blocks;
        ptrdiff_t                                 n_unique_nodes{-1};
        for (auto &block : mesh->blocks()) {
            if (!is_semistructured_type(block->element_type()) || !is_hex_ss_family(block->element_type())) {
                fprintf(stderr, "derefine: only HEX-family semistructured blocks are implemented\n");
                return nullptr;
            }

            const int from_level  = semistructured_level(block->element_type());
            const int step_factor = from_level / to_level;
            const int nxe         = (to_level + 1) * (to_level + 1) * (to_level + 1);

            auto elements = block->elements();
            auto view     = std::make_shared<Buffer<idx_t *>>(
                    nxe,
                    block->n_elements(),
                    (idx_t **)SMESH_ALLOC(nxe * sizeof(idx_t *)),
                    [keep_alive = elements](int, void **v) {
                        (void)keep_alive;
                        SMESH_FREE(v);
                    },
                    elements->mem_space());

            for (int zi = 0; zi <= to_level; zi++) {
                for (int yi = 0; yi <= to_level; yi++) {
                    for (int xi = 0; xi <= to_level; xi++) {
                        const int from_lidx   = sshex8_lidx(from_level, xi * step_factor, yi * step_factor, zi * step_factor);
                        const int to_lidx     = sshex8_lidx(to_level, xi, yi, zi);
                        view->data()[to_lidx] = elements->data()[from_lidx];
                    }
                }
            }

            enum ElemType element_type = semistructured_type(HEX8, to_level);

            auto derefined_block = std::make_shared<Mesh::Block>();
            derefined_block->set_name(block->name());
            derefined_block->set_element_type(element_type);
            derefined_block->set_elements(view);
            blocks.push_back(derefined_block);

            // Find max node id
            {
                auto            vv        = view->data();
                const ptrdiff_t nelements = block->n_elements();
                for (size_t v = 0; v < view->extent(0); v++) {
                    for (ptrdiff_t e = 0; e < nelements; e++) {
                        n_unique_nodes = std::max(static_cast<ptrdiff_t>(vv[v][e]), n_unique_nodes);
                    }
                }
            }
        }

        n_unique_nodes += 1;

        int  sdim   = mesh->spatial_dimension();
        auto points = smesh::view(mesh->points(), 0, sdim, 0, n_unique_nodes);

        return std::make_shared<Mesh>(mesh->comm(), blocks, points);
    }

    SharedBuffer<idx_t *> sshex8_device_elements_view(const SharedBuffer<idx_t *> &fine_device_soa,
                                                      const int                    from_level,
                                                      const int                    to_level) {
        if (!fine_device_soa) {
            SMESH_ERROR("sshex8_device_elements_view: fine_device_soa is null");
            return nullptr;
        }

        if (to_level <= 0 || from_level < to_level || (from_level % to_level) != 0) {
            SMESH_ERROR("sshex8_device_elements_view: invalid levels from=%d to=%d", from_level, to_level);
            return nullptr;
        }

        const int       step_factor = from_level / to_level;
        const int       nxe        = (to_level + 1) * (to_level + 1) * (to_level + 1);
        const ptrdiff_t n_elements = fine_device_soa->extent(1);

        auto fill_host_view = [&](idx_t **host_fine_ptrs, idx_t **host_coarse_ptrs) {
            for (int zi = 0; zi <= to_level; zi++) {
                for (int yi = 0; yi <= to_level; yi++) {
                    for (int xi = 0; xi <= to_level; xi++) {
                        const int from_lidx           = sshex8_lidx(from_level, xi * step_factor, yi * step_factor, zi * step_factor);
                        const int to_lidx             = sshex8_lidx(to_level, xi, yi, zi);
                        host_coarse_ptrs[to_lidx] = host_fine_ptrs[from_lidx];
                    }
                }
            }
        };

#ifdef SMESH_ENABLE_CUDA
        if (fine_device_soa->mem_space() == MEMORY_SPACE_DEVICE) {
            const ptrdiff_t        n0_fine = fine_device_soa->extent(0);
            std::vector<idx_t *> fine_host_ptrs(n0_fine);
            device::device_to_host(static_cast<size_t>(n0_fine), fine_device_soa->data(), fine_host_ptrs.data());

            std::vector<idx_t *> coarse_host_ptrs(static_cast<size_t>(nxe));
            fill_host_view(fine_host_ptrs.data(), coarse_host_ptrs.data());

            idx_t **dev_ptrs = device::alloc<idx_t *>(static_cast<size_t>(nxe));
            device::host_to_device(static_cast<size_t>(nxe), coarse_host_ptrs.data(), dev_ptrs);

            return Buffer<idx_t *>::own(
                    static_cast<size_t>(nxe),
                    static_cast<size_t>(n_elements),
                    dev_ptrs,
                    [keep_alive = fine_device_soa](int /*n*/, void **ptr) {
                        (void)keep_alive;
                        device::destroy(ptr);
                    },
                    MEMORY_SPACE_DEVICE);
        }
#endif

        auto view = Buffer<idx_t *>::own(
                static_cast<size_t>(nxe),
                static_cast<size_t>(n_elements),
                (idx_t **)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(idx_t *)),
                [keep_alive = fine_device_soa](int /*n*/, void **v) {
                    (void)keep_alive;
                    SMESH_FREE(v);
                },
                fine_device_soa->mem_space());

        fill_host_view(fine_device_soa->data(), view->data());
        return view;
    }

    SharedBuffer<idx_t *> sshex8_to_hex8_device_elements_view(const SharedBuffer<idx_t *> &fine_device_soa,
                                                              const int                    from_level) {
        if (!fine_device_soa) {
            SMESH_ERROR("sshex8_to_hex8_device_elements_view: fine_device_soa is null");
            return nullptr;
        }

        if (from_level < 1) {
            SMESH_ERROR("sshex8_to_hex8_device_elements_view: invalid from_level=%d", from_level);
            return nullptr;
        }

        // Corner (xi,yi,zi) for HEX8 local nodes — same as sshex8_to_standard_hex8_mesh(level=1).
        static const int hex8_corner_xi[8] = {0, 1, 1, 0, 0, 1, 1, 0};
        static const int hex8_corner_yi[8] = {0, 0, 1, 1, 0, 0, 1, 1};
        static const int hex8_corner_zi[8] = {0, 0, 0, 0, 1, 1, 1, 1};

        constexpr int   nxe        = 8;
        const int       step_factor = from_level;  // view of level-1 corners in the fine SS mesh
        const ptrdiff_t n_elements = fine_device_soa->extent(1);

        auto fill_host_view = [&](idx_t **host_fine_ptrs, idx_t **host_hex8_ptrs) {
            for (int l = 0; l < nxe; l++) {
                const int from_lidx =
                        sshex8_lidx(from_level, hex8_corner_xi[l] * step_factor, hex8_corner_yi[l] * step_factor, hex8_corner_zi[l] * step_factor);
                host_hex8_ptrs[l] = host_fine_ptrs[from_lidx];
            }
        };

#ifdef SMESH_ENABLE_CUDA
        if (fine_device_soa->mem_space() == MEMORY_SPACE_DEVICE) {
            const ptrdiff_t      n0_fine = fine_device_soa->extent(0);
            std::vector<idx_t *> fine_host_ptrs(static_cast<size_t>(n0_fine));
            device::device_to_host(static_cast<size_t>(n0_fine), fine_device_soa->data(), fine_host_ptrs.data());

            std::vector<idx_t *> hex8_host_ptrs(static_cast<size_t>(nxe));
            fill_host_view(fine_host_ptrs.data(), hex8_host_ptrs.data());

            idx_t **dev_ptrs = device::alloc<idx_t *>(static_cast<size_t>(nxe));
            device::host_to_device(static_cast<size_t>(nxe), hex8_host_ptrs.data(), dev_ptrs);

            return Buffer<idx_t *>::own(
                    static_cast<size_t>(nxe),
                    static_cast<size_t>(n_elements),
                    dev_ptrs,
                    [keep_alive = fine_device_soa](int /*n*/, void **ptr) {
                        (void)keep_alive;
                        device::destroy(ptr);
                    },
                    MEMORY_SPACE_DEVICE);
        }
#endif

        auto view = Buffer<idx_t *>::own(
                static_cast<size_t>(nxe),
                static_cast<size_t>(n_elements),
                (idx_t **)SMESH_ALLOC(static_cast<size_t>(nxe) * sizeof(idx_t *)),
                [keep_alive = fine_device_soa](int /*n*/, void **v) {
                    (void)keep_alive;
                    SMESH_FREE(v);
                },
                fine_device_soa->mem_space());

        fill_host_view(fine_device_soa->data(), view->data());
        return view;
    }

}  // namespace smesh
