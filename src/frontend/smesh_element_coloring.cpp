#include "smesh_element_coloring.hpp"

#include "smesh_buffer.hpp"
#include "smesh_crs_graph.hpp"

#include <algorithm>
#include <vector>

namespace smesh {

    class ElementColoring::Block {
    public:
        std::shared_ptr<Mesh::Block> block;
        SharedBuffer<idx_t>          colors;
        SharedBuffer<ptrdiff_t>      color_ptr;
        int                          n_colors{0};
        ptrdiff_t                    min_per_color{0};
        ptrdiff_t                    max_per_color{0};

        void print(std::ostream &os, const int verbosity) const {
            os << "[ElementColoring] block \"" << (block ? block->name() : std::string("?")) << "\"\n";
            os << "  elements:  " << (block ? block->n_elements() : 0) << "\n";
            os << "  colors:    " << n_colors << "\n";
            os << "  per color: " << min_per_color << " .. " << max_per_color << "\n";
            if (verbosity > 0 && color_ptr) {
                const ptrdiff_t *const ptr = color_ptr->data();
                for (int c = 0; c < n_colors; ++c) {
                    os << "    color " << c << ": [" << ptr[c] << ", " << ptr[c + 1] << ")\n";
                }
            }
        }
    };

    class ElementColoring::Impl {
    public:
        std::shared_ptr<Mesh>                 mesh;
        std::vector<std::shared_ptr<Block>>   blocks;

        void init(const std::shared_ptr<Mesh>                 &in_mesh,
                  const std::vector<std::string>              &block_names,
                  const bool                                   modify_mesh,
                  const std::vector<std::shared_ptr<Sideset>> &sidesets) {
            mesh = in_mesh;
            if (!mesh) {
                return;
            }

            // The conflict relation is "shares a NODE", so the input is the node-to-element graph
            // and not the dual graph: two elements meeting at a single node or along an edge
            // conflict without sharing a face, and a face-adjacency colouring would call them
            // compatible and race. `node_to_element_graph()` is memoized on the mesh and is
            // multiblock-correct, which is why this builds no graph of its own.
            if (!mesh->node_to_element_graph()) {
                return;
            }

            const auto selected = mesh->blocks(block_names);
            blocks.clear();
            blocks.reserve(selected.size());

            for (const auto &mesh_block : selected) {
                if (!mesh_block || !mesh_block->elements()) {
                    continue;
                }
                auto out = std::make_shared<Block>();
                out->block = mesh_block;
                color_one_block(*out, modify_mesh, sidesets);
                blocks.push_back(out);
            }
        }

    private:
        /// The block's position in `Mesh::blocks()`, which is the number the node-to-element
        /// graph tags its entries with and the `block_id` `reorder_elements_from_tags` takes.
        block_idx_t block_number_of(const std::shared_ptr<Mesh::Block> &needle) const {
            const auto &all = mesh->blocks();
            for (size_t i = 0; i < all.size(); ++i) {
                if (all[i] == needle) {
                    return static_cast<block_idx_t>(i);
                }
            }
            return 0;
        }

        void color_one_block(Block                                       &out,
                             const bool                                   modify_mesh,
                             const std::vector<std::shared_ptr<Sideset>> &sidesets) {
            auto graph = mesh->node_to_element_graph();
            const count_t *const       rowptr      = graph->rowptr()->data();
            const element_idx_t *const colidx      = graph->colidx()->data();
            // In a multiblock mesh `colidx` holds BLOCK-LOCAL element indices, and this parallel
            // array says which block each entry belongs to. With a single block it is absent and
            // every entry is this block's.
            auto                     n2e_block   = mesh->node_to_element_block_number();
            const block_idx_t *const d_n2e_block = n2e_block ? n2e_block->data() : nullptr;

            const ptrdiff_t   nelems = out.block->n_elements();
            const int         nxe    = out.block->n_nodes_per_element();
            idx_t *const     *elems  = out.block->elements()->data();
            const block_idx_t bnum   = block_number_of(out.block);

            if (nelems == 0) {
                return;
            }

            // BALANCED GREEDY, VISITING ELEMENTS IN ASCENDING ORDER.
            //
            // The visit order is deliberately NOT degree-sorted. Callers run this after a
            // space-filling reorder, so neighbours are visited close together: the colour count
            // lands near the lower bound, and -- the reason it matters for a sweep -- the elements
            // that end up sharing a colour stay strided through the curve instead of scattered, so
            // a colour's gather keeps most of the locality the element order was built for.
            // Sorting by degree first would give that up for a colour count that is already
            // minimal here.
            //
            // The colour CHOICE is the least loaded feasible one rather than the lowest index.
            // Every colour is a parallel region ending in a barrier, so the classic lowest-index
            // rule -- which leaves the late colours nearly empty -- pays a full barrier for a
            // handful of elements.
            //
            // Serial and deterministic: no randomness, no tie-breaking freedom, and the only
            // input is the element order. Deterministic for a FIXED element order, which is what
            // a reproducibility claim built on this needs; change the ordering and the colouring
            // changes with it, so the colour count is not a property of the mesh.
            std::vector<int>       color(static_cast<size_t>(nelems), -1);
            std::vector<char>      used;
            std::vector<ptrdiff_t> count;
            int                    n_colors = 0;

            for (ptrdiff_t e = 0; e < nelems; ++e) {
                used.assign(static_cast<size_t>(n_colors), 0);
                for (int a = 0; a < nxe; ++a) {
                    const ptrdiff_t node = static_cast<ptrdiff_t>(elems[a][e]);
                    for (count_t j = rowptr[node]; j < rowptr[node + 1]; ++j) {
                        if (d_n2e_block != nullptr && d_n2e_block[j] != bnum) {
                            continue;
                        }
                        const ptrdiff_t other = static_cast<ptrdiff_t>(colidx[j]);
                        if (other < 0 || other >= nelems) {
                            continue;
                        }
                        const int c = color[static_cast<size_t>(other)];
                        if (c >= 0 && c < n_colors) {
                            used[static_cast<size_t>(c)] = 1;
                        }
                    }
                }

                int chosen = -1;
                for (int c = 0; c < n_colors; ++c) {
                    if (used[static_cast<size_t>(c)] != 0) {
                        continue;
                    }
                    if (chosen < 0 || count[static_cast<size_t>(c)] < count[static_cast<size_t>(chosen)]) {
                        chosen = c;
                    }
                }
                if (chosen < 0) {
                    chosen = n_colors++;
                    count.push_back(0);
                }
                color[static_cast<size_t>(e)] = chosen;
                count[static_cast<size_t>(chosen)]++;
            }

            out.n_colors      = n_colors;
            out.min_per_color = count.empty() ? 0 : *std::min_element(count.begin(), count.end());
            out.max_per_color = count.empty() ? 0 : *std::max_element(count.begin(), count.end());

            out.colors = create_host_buffer<idx_t>(static_cast<size_t>(nelems));
            idx_t *const d_colors = out.colors->data();
            for (ptrdiff_t e = 0; e < nelems; ++e) {
                d_colors[e] = static_cast<idx_t>(color[static_cast<size_t>(e)]);
            }

            if (!modify_mesh) {
                return;
            }

            // The renumbering is `Mesh::reorder_elements_from_tags` with the colour as the tag: a
            // counting sort that ALSO remaps registered sidesets and edgesets, which a hand-rolled
            // permutation would silently invalidate. It is stable -- it scans elements ascending
            // and appends -- so the element order inside a colour is the order it came in with,
            // which is what preserves the locality the visit order above was chosen for.
            mesh->reorder_elements_from_tags(bnum, out.colors, sidesets);

            // After that sort colour c is exactly [prefix[c], prefix[c + 1]).
            out.color_ptr = create_host_buffer<ptrdiff_t>(static_cast<size_t>(n_colors) + 1);
            ptrdiff_t *const d_ptr = out.color_ptr->data();
            d_ptr[0]               = 0;
            for (int c = 0; c < n_colors; ++c) {
                d_ptr[c + 1] = d_ptr[c] + count[static_cast<size_t>(c)];
            }
            // `colors` described the OLD numbering; in the new one it is sorted by construction.
            for (int c = 0; c < n_colors; ++c) {
                for (ptrdiff_t e = d_ptr[c]; e < d_ptr[c + 1]; ++e) {
                    d_colors[e] = static_cast<idx_t>(c);
                }
            }
        }
    };

    ElementColoring::ElementColoring() : impl_(std::make_unique<Impl>()) {}
    ElementColoring::~ElementColoring() = default;

    std::shared_ptr<ElementColoring> ElementColoring::create(
            const std::shared_ptr<Mesh>                 &mesh,
            const std::vector<std::string>              &block_names,
            const bool                                   modify_mesh,
            const std::vector<std::shared_ptr<Sideset>> &sidesets) {
        auto coloring = std::make_shared<ElementColoring>();
        coloring->impl_->init(mesh, block_names, modify_mesh, sidesets);
        return coloring;
    }

    std::shared_ptr<Mesh> ElementColoring::mesh() const { return impl_->mesh; }

    ptrdiff_t ElementColoring::n_blocks() const { return static_cast<ptrdiff_t>(impl_->blocks.size()); }

    std::string ElementColoring::block_name(const int block_idx) const {
        return impl_->blocks[static_cast<size_t>(block_idx)]->block->name();
    }

    int ElementColoring::n_colors(const int block_idx) const {
        return impl_->blocks[static_cast<size_t>(block_idx)]->n_colors;
    }

    SharedBuffer<idx_t> ElementColoring::colors(const int block_idx) const {
        return impl_->blocks[static_cast<size_t>(block_idx)]->colors;
    }

    SharedBuffer<ptrdiff_t> ElementColoring::color_ptr(const int block_idx) const {
        return impl_->blocks[static_cast<size_t>(block_idx)]->color_ptr;
    }

    ptrdiff_t ElementColoring::min_elements_per_color(const int block_idx) const {
        return impl_->blocks[static_cast<size_t>(block_idx)]->min_per_color;
    }

    ptrdiff_t ElementColoring::max_elements_per_color(const int block_idx) const {
        return impl_->blocks[static_cast<size_t>(block_idx)]->max_per_color;
    }

    void ElementColoring::print(std::ostream &os, const int verbosity) const {
        for (const auto &block : impl_->blocks) {
            block->print(os, verbosity);
        }
    }

}  // namespace smesh
