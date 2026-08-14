#include <algorithm>
#include <cstring>
#include <cstdint>
#include <vector>

#include "smesh_conversion.hpp"
#include "smesh_mesh.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static std::shared_ptr<Mesh> create_hex8_tet4_serial(const ptrdiff_t nx,
                                                     const ptrdiff_t ny,
                                                     const ptrdiff_t nz) {
    auto            cube      = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
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

static std::shared_ptr<Mesh> split_first_half(const std::shared_ptr<Mesh> &mesh) {
    auto              out     = mesh->clone();
    const ptrdiff_t   n       = out->n_elements(0);
    const ptrdiff_t   n_split = n / 2;
    auto              parents = create_host_buffer<element_idx_t>(static_cast<size_t>(n_split));
    for (ptrdiff_t i = 0; i < n_split; ++i) {
        parents->data()[i] = static_cast<element_idx_t>(i);
    }
    if (out->split_block(parents, "part0", 0) != SMESH_SUCCESS) {
        return nullptr;
    }
    return out;
}

static int node_map_is_permutation(const SharedBuffer<idx_t> &node_map, const ptrdiff_t n_nodes) {
    SMESH_TEST_ASSERT(node_map != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(node_map->size()), n_nodes);
    std::vector<unsigned char> seen(static_cast<size_t>(n_nodes), 0);
    auto                       map = node_map->data();
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        const idx_t dst = map[i];
        SMESH_TEST_ASSERT(dst >= 0);
        SMESH_TEST_ASSERT(dst < n_nodes);
        SMESH_TEST_ASSERT(seen[static_cast<size_t>(dst)] == 0);
        seen[static_cast<size_t>(dst)] = 1;
    }
    return SMESH_TEST_SUCCESS;
}

static int test_checkerboard_to_semistructured() {
    auto comm  = Communicator::self();
    auto multi = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_ASSERT(single != nullptr);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_ASSERT(ss_multi->block(0)->name() == "white");
    SMESH_TEST_ASSERT(ss_multi->block(1)->name() == "black");
    SMESH_TEST_EQ(ss_multi->element_type(0), semistructured_type(HEX8, 2));
    SMESH_TEST_EQ(ss_multi->element_type(1), semistructured_type(HEX8, 2));

    auto hex_multi = sshex_to_hex8(ss_multi);
    SMESH_TEST_ASSERT(hex_multi != nullptr);
    SMESH_TEST_EQ(static_cast<int>(hex_multi->n_blocks()), 2);
    const ptrdiff_t txe = sshex8_txe(2);
    SMESH_TEST_EQ(hex_multi->n_elements(0), multi->n_elements(0) * txe);
    SMESH_TEST_EQ(hex_multi->n_elements(1), multi->n_elements(1) * txe);
    SMESH_TEST_EQ(hex_multi->n_elements(), ss_single->n_elements() * txe);

    return SMESH_TEST_SUCCESS;
}

static int test_checkerboard_to_semistructured_hierarchical() {
    auto comm  = Communicator::self();
    auto multi = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);

    auto ss_multi  = to_semistructured(2, multi, true, false);
    auto ss_single = to_semistructured(2, single, true, false);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    return SMESH_TEST_SUCCESS;
}

static int test_checkerboard_derefine() {
    auto comm  = Communicator::self();
    auto multi = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
    auto single = Mesh::create_hex8_cube(comm, 2, 2, 2);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    auto d_multi   = derefine(ss_multi, 1);
    auto d_single  = derefine(ss_single, 1);
    SMESH_TEST_ASSERT(d_multi != nullptr);
    SMESH_TEST_ASSERT(d_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(d_multi->n_blocks()), 2);
    SMESH_TEST_EQ(d_multi->n_nodes(), d_single->n_nodes());
    SMESH_TEST_EQ(d_multi->element_type(0), semistructured_type(HEX8, 1));
    return SMESH_TEST_SUCCESS;
}

static int test_mixed_to_semistructured_rejected() {
    auto mixed = create_hex8_tet4_serial(2, 2, 2);
    SMESH_TEST_ASSERT(mixed != nullptr);
    SMESH_TEST_ASSERT(to_semistructured(2, mixed) == nullptr);
    return SMESH_TEST_SUCCESS;
}

static int test_tet4_to_semistructured() {
    auto comm  = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);
    SMESH_TEST_ASSERT(single != nullptr);
    SMESH_TEST_ASSERT(multi != nullptr);
    SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

    for (int L : {1, 2, 3, 4}) {
        auto ss_multi  = to_semistructured(L, multi);
        auto ss_single = to_semistructured(L, single);
        SMESH_TEST_ASSERT(ss_multi != nullptr);
        SMESH_TEST_ASSERT(ss_single != nullptr);
        SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
        SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
        SMESH_TEST_EQ(ss_multi->element_type(0), semistructured_type(TET4, L));
        SMESH_TEST_EQ(ss_multi->element_type(1), semistructured_type(TET4, L));
        SMESH_TEST_EQ(static_cast<int>(ss_multi->block(0)->n_nodes_per_element()), sstet4_nxe(L));
        SMESH_TEST_EQ(ss_multi->n_elements(0), multi->n_elements(0));
        SMESH_TEST_EQ(ss_multi->n_elements(1), multi->n_elements(1));
        SMESH_TEST_EQ(ss_multi->n_elements(), single->n_elements());
    }

    return SMESH_TEST_SUCCESS;
}

static int test_tet4_to_semistructured_hierarchical() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);

    auto ss_multi  = to_semistructured(2, multi, true, false);
    auto ss_single = to_semistructured(2, single, true, false);
    SMESH_TEST_ASSERT(ss_multi != nullptr);
    SMESH_TEST_ASSERT(ss_single != nullptr);
    SMESH_TEST_EQ(ss_multi->n_nodes(), ss_single->n_nodes());
    SMESH_TEST_EQ(static_cast<int>(ss_multi->n_blocks()), 2);
    return SMESH_TEST_SUCCESS;
}

static int test_tet4_derefine() {
    auto comm   = Communicator::self();
    auto single = Mesh::create_tet4_cube(comm, 2, 2, 2);
    auto multi  = split_first_half(single);

    auto ss_multi  = to_semistructured(2, multi);
    auto ss_single = to_semistructured(2, single);
    auto d_multi   = derefine(ss_multi, 1);
    auto d_single  = derefine(ss_single, 1);
    SMESH_TEST_ASSERT(d_multi != nullptr);
    SMESH_TEST_ASSERT(d_single != nullptr);
    SMESH_TEST_EQ(static_cast<int>(d_multi->n_blocks()), 2);
    SMESH_TEST_EQ(d_multi->n_nodes(), d_single->n_nodes());
    SMESH_TEST_EQ(d_multi->element_type(0), semistructured_type(TET4, 1));
    SMESH_TEST_EQ(d_multi->n_nodes(), single->n_nodes());
    return SMESH_TEST_SUCCESS;
}

static int test_packed_checkerboard() {
    auto mesh = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto packed = PackedMesh<uint16_t>::create(mesh, {}, false);
    SMESH_TEST_ASSERT(packed != nullptr);
    SMESH_TEST_EQ(packed->n_blocks(), static_cast<ptrdiff_t>(2));
    SMESH_TEST_ASSERT(packed->block_name(0) == "white");
    SMESH_TEST_ASSERT(packed->block_name(1) == "black");
    for (int b = 0; b < 2; ++b) {
        auto elems = packed->elements(b);
        SMESH_TEST_ASSERT(elems != nullptr);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(elems->extent(1)), mesh->n_elements(static_cast<block_idx_t>(b)));
        SMESH_TEST_EQ(static_cast<int>(elems->extent(0)), 8);
    }
    SMESH_TEST_EQ(node_map_is_permutation(packed->node_map(), mesh->n_nodes()), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_packed_hex8_tet4() {
    auto mesh = create_hex8_tet4_serial(2, 2, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto packed = PackedMesh<uint16_t>::create(mesh, {}, false);
    SMESH_TEST_ASSERT(packed != nullptr);
    SMESH_TEST_EQ(packed->n_blocks(), static_cast<ptrdiff_t>(2));
    SMESH_TEST_EQ(static_cast<int>(packed->elements(0)->extent(0)), 8);
    SMESH_TEST_EQ(static_cast<int>(packed->elements(1)->extent(0)), 4);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(packed->elements(0)->extent(1)), mesh->n_elements(0));
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(packed->elements(1)->extent(1)), mesh->n_elements(1));
    SMESH_TEST_EQ(node_map_is_permutation(packed->node_map(), mesh->n_nodes()), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_checkerboard_to_semistructured);
    SMESH_RUN_TEST(test_checkerboard_to_semistructured_hierarchical);
    SMESH_RUN_TEST(test_checkerboard_derefine);
    SMESH_RUN_TEST(test_mixed_to_semistructured_rejected);
    SMESH_RUN_TEST(test_tet4_to_semistructured);
    SMESH_RUN_TEST(test_tet4_to_semistructured_hierarchical);
    SMESH_RUN_TEST(test_tet4_derefine);
    SMESH_RUN_TEST(test_packed_checkerboard);
    SMESH_RUN_TEST(test_packed_hex8_tet4);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

