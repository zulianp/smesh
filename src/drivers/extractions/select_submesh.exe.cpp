#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_extractions.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <stdio.h>
#include <algorithm>
#include <vector>

using namespace smesh;

static std::shared_ptr<Mesh> select_submesh(const Mesh &mesh, const SharedBuffer<element_idx_t> &element_indices) {
    const ptrdiff_t n_blocks            = static_cast<ptrdiff_t>(mesh.n_blocks());
    const ptrdiff_t n_nodes             = mesh.n_nodes();
    const ptrdiff_t n_total_elements    = mesh.n_elements();
    const ptrdiff_t n_selected_elements = static_cast<ptrdiff_t>(element_indices->size());
    const int       dim                 = mesh.spatial_dimension();

    std::vector<ptrdiff_t>     block_offsets(n_blocks + 1, 0);
    std::vector<ptrdiff_t>     selected_per_block(n_blocks, 0);
    std::vector<block_idx_t>   selected_block_ids(n_selected_elements);
    std::vector<element_idx_t> selected_local_ids(n_selected_elements);
    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        block_offsets[b + 1] = block_offsets[b] + mesh.n_elements(static_cast<block_idx_t>(b));
    }

    auto      selected_nodes = create_host_buffer<u8>(n_nodes);
    auto      old_to_new     = create_host_buffer<idx_t>(n_nodes);
    u8       *d_selected     = selected_nodes->data();
    idx_t    *d_old_to_new   = old_to_new->data();
    ptrdiff_t n_selected_nodes{0};

    const element_idx_t *d_selected_elements = element_indices->data();
    for (ptrdiff_t i = 0; i < n_selected_elements; ++i) {
        const ptrdiff_t global_element = static_cast<ptrdiff_t>(d_selected_elements[i]);
        if (global_element < 0 || global_element >= n_total_elements) {
            SMESH_ERROR(
                    "select_submesh: element index %ld out of bounds [0, %ld)\n", (long)global_element, (long)n_total_elements);
            return nullptr;
        }

        const auto      it       = std::upper_bound(block_offsets.begin() + 1, block_offsets.end(), global_element);
        const ptrdiff_t block_id = static_cast<ptrdiff_t>((it - block_offsets.begin()) - 1);
        const ptrdiff_t local_e  = global_element - block_offsets[block_id];
        selected_block_ids[i]    = static_cast<block_idx_t>(block_id);
        selected_local_ids[i]    = static_cast<element_idx_t>(local_e);

        ++selected_per_block[block_id];

        const int           nxe      = mesh.n_nodes_per_element(static_cast<block_idx_t>(block_id));
        const idx_t *const *elements = mesh.elements(static_cast<block_idx_t>(block_id))->data();
        for (int v = 0; v < nxe; ++v) {
            const idx_t node = elements[v][local_e];
            if (!d_selected[node]) {
                d_selected[node] = 1;
                ++n_selected_nodes;
            }
        }
    }

    auto selected_points = create_host_buffer<geom_t>(dim, n_selected_nodes);
    auto node_mapping    = create_host_buffer<idx_t>(n_selected_nodes);

    std::vector<SharedBuffer<idx_t *>>        block_elements(n_blocks);
    std::vector<ptrdiff_t>                    block_write_offset(n_blocks, 0);
    std::vector<std::shared_ptr<Mesh::Block>> out_blocks;
    out_blocks.reserve(n_blocks);

    for (ptrdiff_t b = 0; b < n_blocks; ++b) {
        if (selected_per_block[b] == 0) {
            continue;
        }

        auto elements = create_host_buffer<idx_t>(mesh.n_nodes_per_element(static_cast<block_idx_t>(b)), selected_per_block[b]);
        block_elements[b] = elements;
        out_blocks.push_back(std::make_shared<Mesh::Block>(
                mesh.block(static_cast<size_t>(b))->name(), mesh.element_type(static_cast<block_idx_t>(b)), elements));
    }

    if (out_blocks.empty() && n_blocks > 0) {
        auto elements = create_host_buffer<idx_t>(mesh.n_nodes_per_element(0), 0);
        out_blocks.push_back(std::make_shared<Mesh::Block>(mesh.block(0)->name(), mesh.element_type(0), elements));
    }

    auto      out_points   = selected_points->data();
    auto      src_points   = mesh.points()->data();
    idx_t    *d_mapping    = node_mapping->data();
    ptrdiff_t next_node_id = 0;

    for (ptrdiff_t i = 0; i < n_selected_elements; ++i) {
        const ptrdiff_t block_id = selected_block_ids[i];
        const ptrdiff_t local_e  = selected_local_ids[i];
        const ptrdiff_t out_e    = block_write_offset[block_id]++;

        const int           nxe          = mesh.n_nodes_per_element(static_cast<block_idx_t>(block_id));
        const idx_t *const *src_elements = mesh.elements(static_cast<block_idx_t>(block_id))->data();
        idx_t *const       *dst_elements = block_elements[block_id]->data();

        for (int v = 0; v < nxe; ++v) {
            const idx_t node = src_elements[v][local_e];
            if (d_selected[node] == 1) {
                const idx_t new_node = static_cast<idx_t>(next_node_id++);
                d_selected[node]     = 2;
                d_old_to_new[node]   = new_node;
                d_mapping[new_node]  = node;

                for (int d = 0; d < dim; ++d) {
                    out_points[d][new_node] = src_points[d][node];
                }
            }

            dst_elements[v][out_e] = d_old_to_new[node];
        }
    }

    auto ret = std::make_shared<Mesh>(mesh.comm(), out_blocks, selected_points);
    ret->set_node_mapping(node_mapping);
    return ret;
}

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("select_submesh.exe");
    auto ctx = initialize_serial(argc, argv);

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <mesh_folder> <element_indices.element_idx_t> <output_folder>\n", argv[0]);
        return SMESH_FAILURE;
    }

    auto mesh            = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));
    auto element_indices = Buffer<element_idx_t>::from_file(Path(argv[2]));
    auto output_folder   = Path(argv[3]);

    auto selected_mesh = select_submesh(*mesh, element_indices);
    if (!selected_mesh) {
        return SMESH_FAILURE;
    }

    return selected_mesh->write(output_folder);
}
