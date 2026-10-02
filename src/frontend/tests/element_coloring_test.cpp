#include <algorithm>
#include <cstring>
#include <set>
#include <vector>

#include "smesh_conversion.hpp"
#include "smesh_element_coloring.hpp"
#include "smesh_mesh.hpp"
#include "smesh_mesh_reorder.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

using namespace smesh;

// A hex block and a tet block over the same points, so the two share nodes along the seam. The
// per-block scope of the colouring is only visible on a mesh like this.
static std::shared_ptr<Mesh> create_hex8_tet4_serial(const ptrdiff_t nx, const ptrdiff_t ny, const ptrdiff_t nz) {
    auto            cube       = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
    const ptrdiff_t n_hex_all  = cube->n_elements();
    const ptrdiff_t n_hex_keep = n_hex_all / 2;
    const ptrdiff_t n_hex_conv = n_hex_all - n_hex_keep;
    auto            hex_src    = cube->elements(0)->data();
    auto            hex_keep   = create_host_buffer<idx_t>(8, static_cast<size_t>(n_hex_keep));
    for (int d = 0; d < 8; ++d) {
        std::memcpy(hex_keep->data()[d], hex_src[d], static_cast<size_t>(n_hex_keep) * sizeof(idx_t));
    }
    idx_t *hex_tail[8];
    for (int d = 0; d < 8; ++d) {
        hex_tail[d] = hex_src[d] + n_hex_keep;
    }
    auto tet_buf = create_host_buffer<idx_t>(4, static_cast<size_t>(n_hex_conv * 6));
    mesh_hex8_to_6x_tet4<idx_t>(n_hex_conv, hex_tail, tet_buf->data());
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(std::make_shared<Mesh::Block>("hex", HEX8, hex_keep));
    blocks.push_back(std::make_shared<Mesh::Block>("tet", TET4, tet_buf));
    return std::make_shared<Mesh>(Communicator::self(), blocks, cube->points());
}

// THE CORRECTNESS PROPERTY. No two elements of a colour share a node, so a sweep over one colour
// can write with `+=`. Nothing else the class reports implies this, and a colouring that fails it
// races silently rather than failing, which is why it is checked directly against the
// connectivity rather than against the graph the colouring was built from.
static int colors_are_conflict_free(const std::shared_ptr<Mesh> &mesh,
                                   const block_idx_t            block_id,
                                   const SharedBuffer<idx_t>   &colors) {
    const ptrdiff_t nelems = mesh->n_elements(block_id);
    const int       nxe    = mesh->n_nodes_per_element(block_id);
    idx_t **const   elems  = mesh->elements(block_id)->data();
    const idx_t    *d_col  = colors->data();

    idx_t n_colors = 0;
    for (ptrdiff_t e = 0; e < nelems; ++e) {
        n_colors = std::max(n_colors, static_cast<idx_t>(d_col[e] + 1));
    }

    std::vector<std::set<idx_t>> seen(static_cast<size_t>(n_colors));
    for (ptrdiff_t e = 0; e < nelems; ++e) {
        auto &nodes = seen[static_cast<size_t>(d_col[e])];
        for (int a = 0; a < nxe; ++a) {
            SMESH_TEST_ASSERT(nodes.insert(elems[a][e]).second);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_conflict_free_single_block() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 6, 5, 4);
    auto ec   = ElementColoring::create(mesh);
    SMESH_TEST_ASSERT(ec != nullptr);
    SMESH_TEST_EQ(ec->n_blocks(), static_cast<ptrdiff_t>(1));
    SMESH_TEST_ASSERT(ec->n_colors(0) > 0);

    auto colors = ec->colors(0);
    SMESH_TEST_ASSERT(colors != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(colors->size()), mesh->n_elements(0));
    SMESH_TEST_ASSERT(colors_are_conflict_free(mesh, 0, colors) == SMESH_TEST_SUCCESS);

    // Every element carries exactly one colour in range, and the reported sizes account for all
    // of them -- a colouring that silently left an element at -1 would still be conflict-free.
    std::vector<ptrdiff_t> count(static_cast<size_t>(ec->n_colors(0)), 0);
    const idx_t           *d_col = colors->data();
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        SMESH_TEST_ASSERT(d_col[e] >= 0 && d_col[e] < ec->n_colors(0));
        count[static_cast<size_t>(d_col[e])]++;
    }
    ptrdiff_t total = 0;
    for (const ptrdiff_t c : count) {
        SMESH_TEST_ASSERT(c > 0);
        SMESH_TEST_ASSERT(c >= ec->min_elements_per_color(0));
        SMESH_TEST_ASSERT(c <= ec->max_elements_per_color(0));
        total += c;
    }
    SMESH_TEST_EQ(total, mesh->n_elements(0));

    // Without `modify_mesh` the mesh is untouched, so the colour ranges index `element_order`
    // rather than the block: a device assembly that cannot move its elements reads the colouring
    // in this form.
    auto order     = ec->element_order(0);
    auto color_ptr = ec->color_ptr(0);
    SMESH_TEST_ASSERT(order != nullptr);
    SMESH_TEST_ASSERT(color_ptr != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(order->size()), mesh->n_elements(0));
    SMESH_TEST_EQ(color_ptr->data()[ec->n_colors(0)], mesh->n_elements(0));
    std::vector<unsigned char> seen_elem(static_cast<size_t>(mesh->n_elements(0)), 0);
    for (int c = 0; c < ec->n_colors(0); ++c) {
        for (ptrdiff_t i = color_ptr->data()[c]; i < color_ptr->data()[c + 1]; ++i) {
            const element_idx_t e = order->data()[i];
            SMESH_TEST_ASSERT(e >= 0 && e < mesh->n_elements(0));
            SMESH_TEST_ASSERT(seen_elem[static_cast<size_t>(e)] == 0);
            seen_elem[static_cast<size_t>(e)] = 1;
            SMESH_TEST_EQ(d_col[e], static_cast<idx_t>(c));
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_colors_are_deterministic() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 5, 4, 3);
    auto one  = ElementColoring::create(mesh);
    auto two  = ElementColoring::create(mesh);
    SMESH_TEST_EQ(one->n_colors(0), two->n_colors(0));
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        SMESH_TEST_EQ(one->colors(0)->data()[e], two->colors(0)->data()[e]);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_renumbering_makes_colors_contiguous() {
    auto      mesh   = Mesh::create_hex8_cube(Communicator::self(), 6, 5, 4);
    const int nxe    = mesh->n_nodes_per_element(0);
    const ptrdiff_t nelems = mesh->n_elements(0);

    // The connectivity before the renumbering, as node tuples in element order, so the
    // permutation can be checked element by element rather than only for conserving the set.
    std::vector<std::vector<idx_t>> before(static_cast<size_t>(nelems));
    {
        idx_t **const elems = mesh->elements(0)->data();
        for (ptrdiff_t e = 0; e < nelems; ++e) {
            before[static_cast<size_t>(e)].resize(static_cast<size_t>(nxe));
            for (int a = 0; a < nxe; ++a) {
                before[static_cast<size_t>(e)][static_cast<size_t>(a)] = elems[a][e];
            }
        }
    }

    auto ec = ElementColoring::create(mesh, {}, /*modify_mesh=*/true);
    SMESH_TEST_ASSERT(ec != nullptr);
    SMESH_TEST_EQ(mesh->n_elements(0), nelems);

    auto color_ptr = ec->color_ptr(0);
    SMESH_TEST_ASSERT(color_ptr != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(color_ptr->size()), static_cast<ptrdiff_t>(ec->n_colors(0)) + 1);
    const ptrdiff_t *d_ptr = color_ptr->data();
    SMESH_TEST_EQ(d_ptr[0], static_cast<ptrdiff_t>(0));
    SMESH_TEST_EQ(d_ptr[ec->n_colors(0)], nelems);

    // Colour c is exactly [d_ptr[c], d_ptr[c + 1]) -- the property that lets a sweep run the
    // range with no indirection at all -- and the colour array agrees with it.
    const idx_t *d_col = ec->colors(0)->data();
    for (int c = 0; c < ec->n_colors(0); ++c) {
        SMESH_TEST_ASSERT(d_ptr[c + 1] > d_ptr[c]);
        for (ptrdiff_t e = d_ptr[c]; e < d_ptr[c + 1]; ++e) {
            SMESH_TEST_EQ(d_col[e], static_cast<idx_t>(c));
        }
    }
    SMESH_TEST_ASSERT(colors_are_conflict_free(mesh, 0, ec->colors(0)) == SMESH_TEST_SUCCESS);

    // New element i is old element `element_order[i]`, exactly. This is what makes the two views
    // of the colouring one object rather than two: a caller that cannot move its elements drives
    // a kernel by the order array and sweeps the same elements in the same groups as a caller
    // that let the mesh be renumbered.
    auto order = ec->element_order(0);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(order->size()), nelems);
    idx_t **const after = mesh->elements(0)->data();
    for (ptrdiff_t i = 0; i < nelems; ++i) {
        const element_idx_t src = order->data()[i];
        SMESH_TEST_ASSERT(src >= 0 && src < nelems);
        for (int a = 0; a < nxe; ++a) {
            SMESH_TEST_EQ(after[a][i], before[static_cast<size_t>(src)][static_cast<size_t>(a)]);
        }
    }
    return SMESH_TEST_SUCCESS;
}

// The reason the renumbering goes through `Mesh::reorder_elements_from_tags` rather than a
// hand-rolled permutation: a sideset names elements by number, and a permutation that does not
// remap it leaves it naming different elements with nothing to signal the change.
static int test_sidesets_survive_the_renumbering() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 6, 5, 4);
    auto sides = Sideset::create_from_selector(
            mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-8; });
    SMESH_TEST_ASSERT(!sides.empty());
    auto sideset = sides[0];
    SMESH_TEST_ASSERT(sideset->size() > 0);

    const int     nxe   = mesh->n_nodes_per_element(0);
    idx_t **const elems = mesh->elements(0)->data();

    // What the sideset names NOW, by node tuple rather than by element number, since the number
    // is the thing about to change.
    std::vector<std::pair<std::vector<idx_t>, i16>> named;
    for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
        const element_idx_t  parent = sideset->parent()->data()[i];
        std::vector<idx_t>   nodes(static_cast<size_t>(nxe));
        for (int a = 0; a < nxe; ++a) {
            nodes[static_cast<size_t>(a)] = elems[a][parent];
        }
        named.emplace_back(nodes, sideset->lfi()->data()[i]);
    }

    auto ec = ElementColoring::create(mesh, {}, /*modify_mesh=*/true, {sideset});
    SMESH_TEST_ASSERT(ec != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(named.size()), sideset->size());

    idx_t **const after = mesh->elements(0)->data();
    for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
        const element_idx_t parent = sideset->parent()->data()[i];
        SMESH_TEST_ASSERT(parent >= 0 && parent < mesh->n_elements(0));
        std::vector<idx_t> nodes(static_cast<size_t>(nxe));
        for (int a = 0; a < nxe; ++a) {
            nodes[static_cast<size_t>(a)] = after[a][parent];
        }
        SMESH_TEST_ASSERT(nodes == named[static_cast<size_t>(i)].first);
        SMESH_TEST_EQ(sideset->lfi()->data()[i], named[static_cast<size_t>(i)].second);
    }
    return SMESH_TEST_SUCCESS;
}

// Each block gets its own colouring, and a colour of one block says nothing about the other --
// the precondition the header states, and the reason the conflict check below runs per block.
static int test_multiblock_colors_each_block() {
    auto mesh = create_hex8_tet4_serial(4, 4, 4);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto ec = ElementColoring::create(mesh, {}, /*modify_mesh=*/true);
    SMESH_TEST_ASSERT(ec != nullptr);
    SMESH_TEST_EQ(ec->n_blocks(), static_cast<ptrdiff_t>(2));
    SMESH_TEST_ASSERT(ec->block_name(0) == "hex");
    SMESH_TEST_ASSERT(ec->block_name(1) == "tet");

    for (int b = 0; b < 2; ++b) {
        SMESH_TEST_ASSERT(ec->n_colors(b) > 0);
        SMESH_TEST_ASSERT(colors_are_conflict_free(mesh, static_cast<block_idx_t>(b), ec->colors(b)) ==
                          SMESH_TEST_SUCCESS);
        const ptrdiff_t *d_ptr = ec->color_ptr(b)->data();
        SMESH_TEST_EQ(d_ptr[ec->n_colors(b)], mesh->n_elements(static_cast<block_idx_t>(b)));
    }

    // Naming one block colours that one alone and leaves the other's elements where they are.
    auto only_tet = ElementColoring::create(mesh, {"tet"});
    SMESH_TEST_EQ(only_tet->n_blocks(), static_cast<ptrdiff_t>(1));
    SMESH_TEST_ASSERT(only_tet->block_name(0) == "tet");
    return SMESH_TEST_SUCCESS;
}

// The node-to-element graph is cached on the mesh, indexed by node and VALUED by element, so a
// renumbering leaves it pointing at the wrong elements while it still looks well-formed. Nothing
// about that fails loudly -- it comes back out of `node_to_element_graph()` in place of a correct
// one -- and a second colouring pass over the same mesh is the first thing that reads it.
static int test_cached_graph_follows_the_renumbering() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 5, 4, 3);
    SMESH_TEST_ASSERT(mesh->node_to_element_graph() != nullptr);

    auto first = ElementColoring::create(mesh, {}, /*modify_mesh=*/true);
    SMESH_TEST_ASSERT(first != nullptr);

    // The graph must name, for each node, the elements that contain that node in the numbering
    // the mesh has NOW.
    auto            graph  = mesh->node_to_element_graph();
    const count_t  *rowptr = graph->rowptr()->data();
    const element_idx_t *colidx = graph->colidx()->data();
    const int       nxe    = mesh->n_nodes_per_element(0);
    idx_t **const   elems  = mesh->elements(0)->data();
    const ptrdiff_t nelems = mesh->n_elements(0);
    for (ptrdiff_t node = 0; node < mesh->n_nodes(); ++node) {
        for (count_t j = rowptr[node]; j < rowptr[node + 1]; ++j) {
            const ptrdiff_t e = static_cast<ptrdiff_t>(colidx[j]);
            SMESH_TEST_ASSERT(e >= 0 && e < nelems);
            bool contains = false;
            for (int a = 0; a < nxe; ++a) {
                contains = contains || (elems[a][e] == static_cast<idx_t>(node));
            }
            SMESH_TEST_ASSERT(contains);
        }
    }

    // Which is what makes a second pass meaningful: the mesh is already in colour order, so the
    // colouring is reproduced rather than changed.
    auto second = ElementColoring::create(mesh, {}, /*modify_mesh=*/true);
    SMESH_TEST_EQ(first->n_colors(0), second->n_colors(0));
    for (ptrdiff_t e = 0; e < nelems; ++e) {
        SMESH_TEST_EQ(first->colors(0)->data()[e], second->colors(0)->data()[e]);
    }
    SMESH_TEST_ASSERT(colors_are_conflict_free(mesh, 0, second->colors(0)) == SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

// The colouring is meant to be run AFTER a space-filling reorder, which is the configuration
// production code uses and the one the colour count was tuned for.
static int test_after_sfc_reorder() {
    auto mesh = Mesh::create_hex8_cube(Communicator::self(), 8, 8, 8);
    SFC sfc;
    SMESH_TEST_ASSERT(sfc.reorder(*mesh) == SMESH_SUCCESS);
    auto ec = ElementColoring::create(mesh, {}, /*modify_mesh=*/true);
    SMESH_TEST_ASSERT(ec != nullptr);
    SMESH_TEST_ASSERT(colors_are_conflict_free(mesh, 0, ec->colors(0)) == SMESH_TEST_SUCCESS);
    // A HEX8 element conflicts with the 26 around it, so 8 colours is the lower bound on a
    // regular grid; the balanced greedy must not be paying many more than that.
    SMESH_TEST_ASSERT(ec->n_colors(0) >= 8);
    SMESH_TEST_ASSERT(ec->n_colors(0) <= 16);
    // Balance is what the barrier between colours is paid for; a colour far below the others is
    // a barrier paid for a partly idle sweep.
    SMESH_TEST_ASSERT(ec->max_elements_per_color(0) <= 2 * ec->min_elements_per_color(0));
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_conflict_free_single_block);
    SMESH_RUN_TEST(test_colors_are_deterministic);
    SMESH_RUN_TEST(test_renumbering_makes_colors_contiguous);
    SMESH_RUN_TEST(test_sidesets_survive_the_renumbering);
    SMESH_RUN_TEST(test_multiblock_colors_each_block);
    SMESH_RUN_TEST(test_cached_graph_follows_the_renumbering);
    SMESH_RUN_TEST(test_after_sfc_reorder);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
