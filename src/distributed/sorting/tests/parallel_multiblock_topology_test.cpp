#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <limits>
#include <utility>
#include <vector>

#include "smesh_adjacency.hpp"
#include "smesh_conversion.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_dual_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sfc.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static void append_element_edges(const enum ElemType type, idx_t **elems,
                                 const ptrdiff_t e,
                                 std::vector<std::pair<idx_t, idx_t>> &edges) {
  auto push = [&](idx_t a, idx_t b) {
    if (a > b) {
      std::swap(a, b);
    }
    edges.emplace_back(a, b);
  };

  if (type == HEX8) {
    static const int hex8_edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
        {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (int k = 0; k < 12; ++k) {
      push(elems[hex8_edges[k][0]][e], elems[hex8_edges[k][1]][e]);
    }
  } else if (type == TET4) {
    static const int tet4_edges[6][2] = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (int k = 0; k < 6; ++k) {
      push(elems[tet4_edges[k][0]][e], elems[tet4_edges[k][1]][e]);
    }
  } else if (type == QUAD4 || type == QUADSHELL4) {
    static const int quad4_edges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    for (int k = 0; k < 4; ++k) {
      push(elems[quad4_edges[k][0]][e], elems[quad4_edges[k][1]][e]);
    }
  } else if (type == TRI3 || type == TRISHELL3) {
    static const int tri3_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int k = 0; k < 3; ++k) {
      push(elems[tri3_edges[k][0]][e], elems[tri3_edges[k][1]][e]);
    }
  } else if (type == WEDGE6) {
    static const int wedge6_edges[9][2] = {
        {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    for (int k = 0; k < 9; ++k) {
      push(elems[wedge6_edges[k][0]][e], elems[wedge6_edges[k][1]][e]);
    }
  } else if (type == PYRAMID5) {
    static const int pyramid5_edges[8][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    for (int k = 0; k < 8; ++k) {
      push(elems[pyramid5_edges[k][0]][e], elems[pyramid5_edges[k][1]][e]);
    }
  }
}

static ptrdiff_t count_unique_mesh_edges(const Mesh &mesh) {
  std::vector<std::pair<idx_t, idx_t>> edges;
  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    auto block = mesh.block(b);
    auto elems = block->elements()->data();
    const ptrdiff_t n_e = block->n_elements();
    const enum ElemType type = block->element_type();
    for (ptrdiff_t e = 0; e < n_e; ++e) {
      append_element_edges(type, elems, e, edges);
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  return static_cast<ptrdiff_t>(edges.size());
}

static int check_half_face_tables(Mesh &mesh, const int require_inter_block) {
  int n_inter_block = 0;
  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    const block_idx_t bid = static_cast<block_idx_t>(b);
    const enum ElemType type = mesh.element_type(bid);
    const int ns = elem_num_sides(type);
    const ptrdiff_t n_e = mesh.n_elements(bid);
    auto hft = mesh.half_face_table(bid);
    auto hnbb = mesh.half_face_neighbor_block(bid);
    SMESH_TEST_ASSERT(hft != nullptr);
    SMESH_TEST_ASSERT(hnbb != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(hft->size()), n_e * ns);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(hnbb->size()), n_e * ns);

    for (ptrdiff_t e = 0; e < n_e; ++e) {
      for (int s = 0; s < ns; ++s) {
        const element_idx_t neighbor = hft->data()[e * ns + s];
        const block_idx_t nb = hnbb->data()[e * ns + s];
        if (neighbor != invalid_idx<element_idx_t>()) {
          SMESH_TEST_ASSERT(nb >= 0);
          SMESH_TEST_ASSERT(static_cast<size_t>(nb) < mesh.n_blocks());
          SMESH_TEST_ASSERT(neighbor < mesh.n_elements(nb));
          if (nb == bid) {
            SMESH_TEST_ASSERT(neighbor < n_e);
          } else {
            n_inter_block++;
          }
        }
      }
    }
  }
  if (require_inter_block) {
    SMESH_TEST_ASSERT(n_inter_block > 0);
  }
  return SMESH_TEST_SUCCESS;
}

static int meshes_equal(const Mesh &a, const Mesh &b) {
  SMESH_TEST_EQ(a.n_nodes(), b.n_nodes());
  SMESH_TEST_EQ(a.spatial_dimension(), b.spatial_dimension());
  SMESH_TEST_EQ(static_cast<int>(a.n_blocks()), static_cast<int>(b.n_blocks()));

  const int dim = a.spatial_dimension();
  auto ap = a.points()->data();
  auto bp = b.points()->data();
  for (int d = 0; d < dim; ++d) {
    for (ptrdiff_t i = 0; i < a.n_nodes(); ++i) {
      SMESH_TEST_ASSERT(ap[d][i] == bp[d][i]);
    }
  }

  for (size_t bidx = 0; bidx < a.n_blocks(); ++bidx) {
    auto ab = a.block(bidx);
    auto bb = b.block(bidx);
    SMESH_TEST_ASSERT(ab->name() == bb->name());
    SMESH_TEST_EQ(ab->element_type(), bb->element_type());
    SMESH_TEST_EQ(ab->n_elements(), bb->n_elements());
    const int nxe = ab->n_nodes_per_element();
    auto ae = ab->elements()->data();
    auto be = bb->elements()->data();
    for (int d = 0; d < nxe; ++d) {
      for (ptrdiff_t e = 0; e < ab->n_elements(); ++e) {
        SMESH_TEST_EQ(ae[d][e], be[d][e]);
      }
    }
  }
  return SMESH_TEST_SUCCESS;
}

static std::shared_ptr<Mesh>
split_first_half(const std::shared_ptr<Mesh> &mesh) {
  auto out = mesh->clone();
  const ptrdiff_t n = out->n_elements(0);
  const ptrdiff_t n_split = n / 2;
  auto parents = create_host_buffer<element_idx_t>(static_cast<size_t>(n_split));
  for (ptrdiff_t i = 0; i < n_split; ++i) {
    parents->data()[i] = static_cast<element_idx_t>(i);
  }
  if (out->split_block(parents, "part0", 0) != SMESH_SUCCESS) {
    return nullptr;
  }
  return out;
}

static std::shared_ptr<Mesh> create_edge2_line(const ptrdiff_t n_seg, const enum ElemType et = EDGE2) {
  auto elems = create_host_buffer<idx_t>(2, static_cast<size_t>(n_seg));
  auto pts   = create_host_buffer<geom_t>(2, static_cast<size_t>(n_seg + 1));
  for (ptrdiff_t i = 0; i < n_seg; ++i) {
    elems->data()[0][i] = static_cast<idx_t>(i);
    elems->data()[1][i] = static_cast<idx_t>(i + 1);
  }
  for (ptrdiff_t i = 0; i <= n_seg; ++i) {
    pts->data()[0][i] = static_cast<geom_t>(i);
    pts->data()[1][i] = 0;
  }
  return std::make_shared<Mesh>(Communicator::self(), et, elems, pts);
}

static std::shared_ptr<Mesh> create_hex8_tet4_serial(const ptrdiff_t nx,
                                                     const ptrdiff_t ny,
                                                     const ptrdiff_t nz) {
  auto cube = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
  const ptrdiff_t n_hex_all = cube->n_elements();
  const ptrdiff_t n_hex_keep = n_hex_all / 2;
  const ptrdiff_t n_hex_conv = n_hex_all - n_hex_keep;
  auto hex_src = cube->elements(0)->data();
  auto hex_keep = create_host_buffer<idx_t>(8, static_cast<size_t>(n_hex_keep));
  for (int d = 0; d < 8; ++d) {
    std::memcpy(hex_keep->data()[d], hex_src[d],
                static_cast<size_t>(n_hex_keep) * sizeof(idx_t));
  }
  idx_t *hex_tail[8];
  for (int d = 0; d < 8; ++d) {
    hex_tail[d] = hex_src[d] + n_hex_keep;
  }
  auto tet_buf =
      create_host_buffer<idx_t>(4, static_cast<size_t>(n_hex_conv * 6));
  mesh_hex8_to_6x_tet4<idx_t>(n_hex_conv, hex_tail, tet_buf->data());
  std::vector<std::shared_ptr<Mesh::Block>> blocks;
  blocks.push_back(std::make_shared<Mesh::Block>("hex", HEX8, hex_keep));
  blocks.push_back(std::make_shared<Mesh::Block>("tet", TET4, tet_buf));
  return std::make_shared<Mesh>(Communicator::self(), blocks, cube->points());
}

static int test_serial_checkerboard_graphs() {
  auto comm = Communicator::self();
  auto mesh = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
  SMESH_TEST_EQ(check_half_face_tables(*mesh, 1), SMESH_TEST_SUCCESS);

  auto edge = mesh->edge_graph();
  SMESH_TEST_ASSERT(edge != nullptr);
  SMESH_TEST_EQ(edge->nnz(), count_unique_mesh_edges(*mesh));

  auto n2n_hex = mesh->create_node_to_node_graph(HEX8);
  auto n2n = mesh->node_to_node_graph();
  SMESH_TEST_ASSERT(n2n_hex != nullptr);
  SMESH_TEST_ASSERT(n2n != nullptr);
  SMESH_TEST_EQ(n2n_hex->nnz(), n2n->nnz());

  auto n2e = mesh->node_to_element_graph();
  auto n2e_block = mesh->node_to_element_block_number();
  SMESH_TEST_ASSERT(n2e != nullptr);
  SMESH_TEST_ASSERT(n2e_block != nullptr);
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(n2e_block->size()), n2e->nnz());

  auto dual = DualGraph::create(mesh);
  SMESH_TEST_ASSERT(dual != nullptr);
  SMESH_TEST_ASSERT(dual->adj_block() != nullptr);
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(dual->adj_block()->size()),
                static_cast<ptrdiff_t>(dual->adj_idx()->size()));

  return SMESH_TEST_SUCCESS;
}

static int test_serial_hex8_tet4_graphs() {
  auto mesh = create_hex8_tet4_serial(2, 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 2);
  SMESH_TEST_EQ(mesh->element_type(0), HEX8);
  SMESH_TEST_EQ(mesh->element_type(1), TET4);
  SMESH_TEST_EQ(check_half_face_tables(*mesh, 0), SMESH_TEST_SUCCESS);

  auto edge = mesh->edge_graph();
  SMESH_TEST_ASSERT(edge != nullptr);
  SMESH_TEST_EQ(edge->nnz(), count_unique_mesh_edges(*mesh));

  auto n2e = mesh->node_to_element_graph();
  auto n2e_block = mesh->node_to_element_block_number();
  SMESH_TEST_ASSERT(n2e != nullptr);
  SMESH_TEST_ASSERT(n2e_block != nullptr);
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(n2e_block->size()), n2e->nnz());

  auto dual = DualGraph::create(mesh);
  SMESH_TEST_ASSERT(dual != nullptr);
  SMESH_TEST_ASSERT(dual->adj_block() != nullptr);

  return SMESH_TEST_SUCCESS;
}

static int test_serial_transforms() {
  auto comm = Communicator::self();
  auto mesh = Mesh::create_hex8_checkerboard_cube(comm, 2, 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto tet_mesh = convert_to(TET4, mesh);
  SMESH_TEST_ASSERT(tet_mesh != nullptr);
  SMESH_TEST_EQ(static_cast<int>(tet_mesh->n_blocks()), 2);
  SMESH_TEST_EQ(tet_mesh->n_nodes(), mesh->n_nodes());

  auto cloned = mesh->clone();
  SMESH_TEST_ASSERT(cloned != nullptr);
  SMESH_TEST_EQ(meshes_equal(*mesh, *cloned), SMESH_TEST_SUCCESS);

  auto promoted = promote_to(TET10, tet_mesh);
  SMESH_TEST_ASSERT(promoted != nullptr);
  SMESH_TEST_ASSERT(promoted->n_nodes() > tet_mesh->n_nodes());
  SMESH_TEST_EQ(static_cast<int>(promoted->n_blocks()), 2);

  auto refined = refine(tet_mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(static_cast<int>(refined->n_blocks()), 2);

  auto skins = skin_sidesets(mesh);
  SMESH_TEST_EQ(static_cast<int>(skins.size()), 2);
  for (block_idx_t b = 0; b < 2; ++b) {
    SMESH_TEST_EQ(skins[b]->block_id(), b);
    const int ns = elem_num_sides(HEX8);
    for (ptrdiff_t i = 0; i < skins[b]->size(); ++i) {
      const element_idx_t parent = skins[b]->parent()->data()[i];
      const i16 side = skins[b]->lfi()->data()[i];
      const element_idx_t neighbor =
          mesh->half_face_table(b)->data()[parent * ns + side];
      SMESH_TEST_EQ(neighbor, invalid_idx<element_idx_t>());
    }
  }

  return SMESH_TEST_SUCCESS;
}

static int test_serial_promote_refine_vs_single_block() {
  auto single = Mesh::create_tet4_cube(Communicator::self(), 2, 2, 2);
  SMESH_TEST_ASSERT(single != nullptr);
  auto multi = split_first_half(single);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);
  SMESH_TEST_EQ(multi->n_nodes(), single->n_nodes());
  SMESH_TEST_EQ(multi->n_elements(), single->n_elements());

  auto single_p = promote_to(TET10, single);
  auto multi_p = promote_to(TET10, multi);
  SMESH_TEST_ASSERT(single_p != nullptr);
  SMESH_TEST_ASSERT(multi_p != nullptr);
  SMESH_TEST_EQ(multi_p->n_nodes(), single_p->n_nodes());
  SMESH_TEST_EQ(multi_p->n_elements(), single_p->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_p->n_blocks()), 2);

  auto single_r = refine(single, 1);
  auto multi_r = refine(multi, 1);
  SMESH_TEST_ASSERT(single_r != nullptr);
  SMESH_TEST_ASSERT(multi_r != nullptr);
  SMESH_TEST_EQ(multi_r->n_nodes(), single_r->n_nodes());
  SMESH_TEST_EQ(multi_r->n_elements(), single_r->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_r->n_blocks()), 2);

  return SMESH_TEST_SUCCESS;
}

static int test_serial_hex8_multiblock_refine() {
  auto single = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
  auto multi = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 2, 2, 2);
  SMESH_TEST_ASSERT(single != nullptr);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_r = refine(single, 1);
  auto multi_r = refine(multi, 1);
  SMESH_TEST_ASSERT(single_r != nullptr);
  SMESH_TEST_ASSERT(multi_r != nullptr);
  SMESH_TEST_EQ(multi_r->n_nodes(), single_r->n_nodes());
  SMESH_TEST_EQ(multi_r->n_elements(), single_r->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_r->n_blocks()), 2);
  SMESH_TEST_EQ(multi_r->element_type(0), HEX8);
  SMESH_TEST_EQ(multi_r->element_type(1), HEX8);

  return SMESH_TEST_SUCCESS;
}

static int test_serial_quad4_refine() {
  auto mesh = Mesh::create_quad4_square(Communicator::self(), 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  const ptrdiff_t n_e = mesh->n_elements();
  auto            ss  = to_semistructured(2, mesh);
  SMESH_TEST_ASSERT(ss != nullptr);
  auto exploded = ssquad_to_quad4(ss);
  SMESH_TEST_ASSERT(exploded != nullptr);

  auto refined = refine(mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(refined->element_type(0), QUAD4);
  SMESH_TEST_EQ(refined->n_elements(), n_e * 4);
  SMESH_TEST_EQ(refined->n_nodes(), exploded->n_nodes());
  SMESH_TEST_EQ(refined->n_elements(), exploded->n_elements());
  return SMESH_TEST_SUCCESS;
}

static int test_serial_quad4_multiblock_refine() {
  auto single = Mesh::create_quad4_square(Communicator::self(), 4, 4);
  auto multi  = split_first_half(single);
  SMESH_TEST_ASSERT(single != nullptr);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_r = refine(single, 1);
  auto multi_r  = refine(multi, 1);
  SMESH_TEST_ASSERT(single_r != nullptr);
  SMESH_TEST_ASSERT(multi_r != nullptr);
  SMESH_TEST_EQ(multi_r->n_nodes(), single_r->n_nodes());
  SMESH_TEST_EQ(multi_r->n_elements(), single_r->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_r->n_blocks()), 2);
  SMESH_TEST_EQ(multi_r->element_type(0), QUAD4);
  SMESH_TEST_EQ(multi_r->element_type(1), QUAD4);
  return SMESH_TEST_SUCCESS;
}

static int test_serial_quadshell4_refine() {
  auto mesh = Mesh::create_quad4_square(Communicator::self(), 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  mesh->set_element_type(0, QUADSHELL4);
  SMESH_TEST_EQ(mesh->element_type(0), QUADSHELL4);
  auto refined = refine(mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(refined->element_type(0), QUADSHELL4);
  SMESH_TEST_EQ(refined->n_elements(), mesh->n_elements() * 4);
  return SMESH_TEST_SUCCESS;
}

static std::shared_ptr<Mesh> create_wedge6_serial(const ptrdiff_t nx, const ptrdiff_t ny) {
  auto tri = Mesh::create_tri3_square(Communicator::self(), nx, ny);
  return extrude(tri, 1.0, 1);
}

static int test_serial_wedge6_refine() {
  auto mesh = create_wedge6_serial(2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(mesh->element_type(0), WEDGE6);
  const ptrdiff_t n_e = mesh->n_elements();
  auto            ss  = to_semistructured(2, mesh);
  SMESH_TEST_ASSERT(ss != nullptr);
  auto exploded = sswedge_to_wedge6(ss);
  SMESH_TEST_ASSERT(exploded != nullptr);

  auto refined = refine(mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(refined->element_type(0), WEDGE6);
  SMESH_TEST_EQ(refined->n_elements(), n_e * 8);
  SMESH_TEST_EQ(refined->n_nodes(), exploded->n_nodes());
  SMESH_TEST_EQ(refined->n_elements(), exploded->n_elements());
  return SMESH_TEST_SUCCESS;
}

static int test_serial_wedge6_multiblock_refine() {
  auto single = create_wedge6_serial(4, 4);
  auto multi  = split_first_half(single);
  SMESH_TEST_ASSERT(single != nullptr);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_r = refine(single, 1);
  auto multi_r  = refine(multi, 1);
  SMESH_TEST_ASSERT(single_r != nullptr);
  SMESH_TEST_ASSERT(multi_r != nullptr);
  SMESH_TEST_EQ(multi_r->n_nodes(), single_r->n_nodes());
  SMESH_TEST_EQ(multi_r->n_elements(), single_r->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_r->n_blocks()), 2);
  SMESH_TEST_EQ(multi_r->element_type(0), WEDGE6);
  SMESH_TEST_EQ(multi_r->element_type(1), WEDGE6);
  return SMESH_TEST_SUCCESS;
}

static int test_serial_trishell3_refine() {
  auto mesh = Mesh::create_tri3_square(Communicator::self(), 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  const ptrdiff_t n_e = mesh->n_elements();
  mesh->set_element_type(0, TRISHELL3);
  auto refined = refine(mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(refined->element_type(0), TRISHELL3);
  SMESH_TEST_EQ(refined->n_elements(), n_e * 4);
  return SMESH_TEST_SUCCESS;
}

static int test_serial_trishell3_multiblock_refine() {
  auto single = Mesh::create_tri3_square(Communicator::self(), 4, 4);
  SMESH_TEST_ASSERT(single != nullptr);
  single->set_element_type(0, TRISHELL3);
  auto multi = split_first_half(single);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_r = refine(single, 1);
  auto multi_r  = refine(multi, 1);
  SMESH_TEST_ASSERT(single_r != nullptr);
  SMESH_TEST_ASSERT(multi_r != nullptr);
  SMESH_TEST_EQ(multi_r->n_nodes(), single_r->n_nodes());
  SMESH_TEST_EQ(multi_r->n_elements(), single_r->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_r->n_blocks()), 2);
  SMESH_TEST_EQ(multi_r->element_type(0), TRISHELL3);
  SMESH_TEST_EQ(multi_r->element_type(1), TRISHELL3);
  return SMESH_TEST_SUCCESS;
}

static int test_serial_edge2_refine() {
  auto mesh = create_edge2_line(4);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(mesh->element_type(0), EDGE2);
  const ptrdiff_t n_e = mesh->n_elements();
  const ptrdiff_t n_n = mesh->n_nodes();
  auto            refined = refine(mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(refined->element_type(0), EDGE2);
  SMESH_TEST_EQ(refined->n_elements(), n_e * 2);
  SMESH_TEST_EQ(refined->n_nodes(), n_n + n_e);
  return SMESH_TEST_SUCCESS;
}

static int test_serial_edgeshell2_refine() {
  auto mesh = create_edge2_line(4, EDGESHELL2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(mesh->element_type(0), EDGESHELL2);
  auto refined = refine(mesh, 1);
  SMESH_TEST_ASSERT(refined != nullptr);
  SMESH_TEST_EQ(refined->element_type(0), EDGESHELL2);
  SMESH_TEST_EQ(refined->n_elements(), mesh->n_elements() * 2);
  return SMESH_TEST_SUCCESS;
}

static int test_serial_edge2_multiblock_refine() {
  auto single = create_edge2_line(8);
  auto multi  = split_first_half(single);
  SMESH_TEST_ASSERT(single != nullptr);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_r = refine(single, 1);
  auto multi_r  = refine(multi, 1);
  SMESH_TEST_ASSERT(single_r != nullptr);
  SMESH_TEST_ASSERT(multi_r != nullptr);
  SMESH_TEST_EQ(multi_r->n_nodes(), single_r->n_nodes());
  SMESH_TEST_EQ(multi_r->n_elements(), single_r->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_r->n_blocks()), 2);
  SMESH_TEST_EQ(multi_r->element_type(0), EDGE2);
  SMESH_TEST_EQ(multi_r->element_type(1), EDGE2);
  return SMESH_TEST_SUCCESS;
}

static std::shared_ptr<Mesh> create_quad4_square_3d(const ptrdiff_t nx,
                                                    const ptrdiff_t ny) {
  auto q2 = Mesh::create_quad4_square(Communicator::self(), nx, ny);
  auto p3 = create_host_buffer<geom_t>(3, static_cast<size_t>(q2->n_nodes()));
  std::memcpy(p3->data()[0], q2->points()->data()[0],
              static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
  std::memcpy(p3->data()[1], q2->points()->data()[1],
              static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
  std::memset(p3->data()[2], 0, static_cast<size_t>(q2->n_nodes()) * sizeof(geom_t));
  return std::make_shared<Mesh>(Communicator::self(), q2->blocks(), p3);
}

static int test_serial_extrude_vs_single_block() {
  auto single = create_quad4_square_3d(4, 4);
  SMESH_TEST_ASSERT(single != nullptr);
  auto multi = split_first_half(single);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_e = extrude(single, 1.0, 2);
  auto multi_e = extrude(multi, 1.0, 2);
  SMESH_TEST_ASSERT(single_e != nullptr);
  SMESH_TEST_ASSERT(multi_e != nullptr);
  SMESH_TEST_EQ(multi_e->n_nodes(), single_e->n_nodes());
  SMESH_TEST_EQ(multi_e->n_elements(), single_e->n_elements());
  SMESH_TEST_EQ(static_cast<int>(multi_e->n_blocks()), 2);
  SMESH_TEST_EQ(multi_e->element_type(0), HEX8);
  SMESH_TEST_EQ(multi_e->element_type(1), HEX8);

  return SMESH_TEST_SUCCESS;
}

static std::shared_ptr<Mesh> create_n_wedge6(const ptrdiff_t n_e) {
  auto points = create_host_buffer<geom_t>(3, n_e * 6);
  auto elems = create_host_buffer<idx_t>(6, n_e);
  for (ptrdiff_t e = 0; e < n_e; ++e) {
    for (int d = 0; d < 6; ++d) {
      const idx_t node = static_cast<idx_t>(e * 6 + d);
      elems->data()[d][e] = node;
      points->data()[0][node] = static_cast<geom_t>(d % 3);
      points->data()[1][node] = static_cast<geom_t>(d / 3);
      points->data()[2][node] = static_cast<geom_t>(e);
    }
  }
  std::vector<std::shared_ptr<Mesh::Block>> blocks;
  blocks.push_back(std::make_shared<Mesh::Block>("wedge", WEDGE6, elems));
  return std::make_shared<Mesh>(Communicator::self(), blocks, points);
}

static std::shared_ptr<Mesh> create_n_pyramid5(const ptrdiff_t n_e) {
  auto points = create_host_buffer<geom_t>(3, n_e * 5);
  auto elems = create_host_buffer<idx_t>(5, n_e);
  for (ptrdiff_t e = 0; e < n_e; ++e) {
    for (int d = 0; d < 5; ++d) {
      const idx_t node = static_cast<idx_t>(e * 5 + d);
      elems->data()[d][e] = node;
      points->data()[0][node] = static_cast<geom_t>(d);
      points->data()[1][node] = static_cast<geom_t>(e);
      points->data()[2][node] = 0;
    }
  }
  std::vector<std::shared_ptr<Mesh::Block>> blocks;
  blocks.push_back(std::make_shared<Mesh::Block>("pyramid", PYRAMID5, elems));
  return std::make_shared<Mesh>(Communicator::self(), blocks, points);
}

static int test_serial_hex_dominant() {
  auto mesh = Mesh::create_hex_dominant_serial(Communicator::self());
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 4);
  SMESH_TEST_EQ(mesh->n_nodes(), static_cast<ptrdiff_t>(12));
  SMESH_TEST_EQ(mesh->element_type(0), HEX8);
  SMESH_TEST_EQ(mesh->element_type(1), PYRAMID5);
  SMESH_TEST_EQ(mesh->element_type(2), TET4);
  SMESH_TEST_EQ(mesh->element_type(3), WEDGE6);
  SMESH_TEST_EQ(mesh->n_elements(0), static_cast<ptrdiff_t>(1));
  SMESH_TEST_EQ(mesh->n_elements(1), static_cast<ptrdiff_t>(1));
  SMESH_TEST_EQ(mesh->n_elements(2), static_cast<ptrdiff_t>(1));
  SMESH_TEST_EQ(mesh->n_elements(3), static_cast<ptrdiff_t>(1));

  SMESH_TEST_EQ(check_half_face_tables(*mesh, 1), SMESH_TEST_SUCCESS);

  auto hex_hft = mesh->half_face_table(0);
  auto hex_hnb = mesh->half_face_neighbor_block(0);
  SMESH_TEST_EQ(hex_hnb->data()[0], static_cast<block_idx_t>(3));
  SMESH_TEST_EQ(hex_hft->data()[0], static_cast<element_idx_t>(0));
  SMESH_TEST_EQ(hex_hnb->data()[1], static_cast<block_idx_t>(1));
  SMESH_TEST_EQ(hex_hft->data()[1], static_cast<element_idx_t>(0));

  auto pyr_hft = mesh->half_face_table(1);
  auto pyr_hnb = mesh->half_face_neighbor_block(1);
  SMESH_TEST_EQ(pyr_hnb->data()[0], static_cast<block_idx_t>(2));
  SMESH_TEST_EQ(pyr_hft->data()[0], static_cast<element_idx_t>(0));
  SMESH_TEST_EQ(pyr_hnb->data()[4], static_cast<block_idx_t>(0));
  SMESH_TEST_EQ(pyr_hft->data()[4], static_cast<element_idx_t>(0));

  auto edge = mesh->edge_graph();
  SMESH_TEST_ASSERT(edge != nullptr);
  SMESH_TEST_EQ(edge->nnz(), count_unique_mesh_edges(*mesh));

  auto skins = skin_sidesets(mesh);
  SMESH_TEST_ASSERT(!skins.empty());
  int n_ss_block[4] = {0, 0, 0, 0};
  int n_tri_ss = 0;
  int n_quad_ss = 0;
  for (const auto &ss : skins) {
    SMESH_TEST_ASSERT(ss != nullptr);
    const block_idx_t bid = ss->block_id();
    SMESH_TEST_ASSERT(bid >= 0 && bid < 4);
    n_ss_block[bid]++;
    if (ss->size() == 0) {
      continue;
    }
    const enum ElemType et = mesh->element_type(bid);
    const enum ElemType st0 = side_type(et, ss->lfi()->data()[0]);
    for (ptrdiff_t i = 0; i < ss->size(); ++i) {
      SMESH_TEST_EQ(side_type(et, ss->lfi()->data()[i]), st0);
    }
    if (st0 == TRI3) {
      n_tri_ss++;
    } else if (st0 == QUAD4) {
      n_quad_ss++;
    }
  }
  SMESH_TEST_EQ(n_ss_block[0], 1);
  SMESH_TEST_EQ(n_ss_block[2], 1);
  SMESH_TEST_ASSERT(n_ss_block[1] >= 1);
  SMESH_TEST_ASSERT(n_ss_block[3] >= 1);
  SMESH_TEST_ASSERT(n_tri_ss >= 1);
  SMESH_TEST_ASSERT(n_quad_ss >= 1);

  auto wedges = create_n_wedge6(2);
  auto wedge_tets = convert_to(TET4, wedges);
  SMESH_TEST_ASSERT(wedge_tets != nullptr);
  SMESH_TEST_EQ(wedge_tets->n_elements(), static_cast<ptrdiff_t>(6));
  auto wt = wedge_tets->elements(0)->data();
  SMESH_TEST_EQ(wt[0][0], static_cast<idx_t>(0));
  SMESH_TEST_EQ(wt[0][3], static_cast<idx_t>(6));

  auto pyramids = create_n_pyramid5(2);
  auto pyr_tets = convert_to(TET4, pyramids);
  SMESH_TEST_ASSERT(pyr_tets != nullptr);
  SMESH_TEST_EQ(pyr_tets->n_elements(), static_cast<ptrdiff_t>(4));
  auto pt = pyr_tets->elements(0)->data();
  SMESH_TEST_EQ(pt[0][0], static_cast<idx_t>(0));
  SMESH_TEST_EQ(pt[0][2], static_cast<idx_t>(5));

  auto ss = to_semistructured(2, mesh);
  SMESH_TEST_ASSERT(ss != nullptr);
  SMESH_TEST_EQ(static_cast<int>(ss->n_blocks()), 4);
  SMESH_TEST_EQ(ss->element_type(0), semistructured_type(HEX8, 2));
  SMESH_TEST_EQ(ss->element_type(1), semistructured_type(PYRAMID5, 2));
  SMESH_TEST_EQ(ss->element_type(2), semistructured_type(TET4, 2));
  SMESH_TEST_EQ(ss->element_type(3), semistructured_type(WEDGE6, 2));
  SMESH_TEST_ASSERT(ss->n_nodes() > mesh->n_nodes());

  return SMESH_TEST_SUCCESS;
}

#ifdef SMESH_ENABLE_MPI
static int check_sideset_sfc_identity(const Mesh &mesh, const Sideset &sideset,
                                      const element_idx_t *const serial_parent,
                                      const i16 *const serial_lfi,
                                      const ptrdiff_t n_serial,
                                      const ptrdiff_t n_serial_elements,
                                      const int require_not_naive) {
  auto block = mesh.block(sideset.block_id());
  SMESH_TEST_ASSERT(block != nullptr);
  const ptrdiff_t n_owned = block->n_elements_owned();
  auto mapping = block->element_mapping();
  SMESH_TEST_ASSERT(mapping != nullptr);

  std::vector<int> claimed(static_cast<size_t>(n_serial), 0);
  for (ptrdiff_t i = 0; i < sideset.size(); ++i) {
    const element_idx_t local_parent = sideset.parent()->data()[i];
    SMESH_TEST_ASSERT(local_parent >= 0);
    SMESH_TEST_ASSERT(local_parent < n_owned);
    const large_idx_t serial_id = mapping->data()[local_parent];
    SMESH_TEST_ASSERT(serial_id >= 0);
    SMESH_TEST_ASSERT(serial_id < static_cast<large_idx_t>(n_serial_elements));

    const i16 lfi = sideset.lfi()->data()[i];
    bool found = false;
    for (ptrdiff_t j = 0; j < n_serial; ++j) {
      if (claimed[static_cast<size_t>(j)]) {
        continue;
      }
      if (static_cast<large_idx_t>(serial_parent[j]) == serial_id &&
          serial_lfi[j] == lfi) {
        claimed[static_cast<size_t>(j)] = 1;
        found = true;
        break;
      }
    }
    SMESH_TEST_ASSERT(found);
  }

  MPI_Allreduce(MPI_IN_PLACE, claimed.data(), static_cast<int>(n_serial),
                MPI_INT, MPI_SUM, mesh.comm()->get());
  ptrdiff_t n_claimed = 0;
  for (ptrdiff_t j = 0; j < n_serial; ++j) {
    SMESH_TEST_EQ(claimed[static_cast<size_t>(j)], 1);
    n_claimed += claimed[static_cast<size_t>(j)];
  }
  SMESH_TEST_EQ(n_claimed, n_serial);

  ptrdiff_t local_n = sideset.size();
  ptrdiff_t global_n = 0;
  MPI_Allreduce(&local_n, &global_n, 1, mpi_type<ptrdiff_t>(), MPI_SUM,
                mesh.comm()->get());
  SMESH_TEST_EQ(global_n, n_serial);

  int local_not_naive = 0;
  for (ptrdiff_t i = 0; i < n_owned; ++i) {
    const ptrdiff_t serial_id = static_cast<ptrdiff_t>(mapping->data()[i]);
    if (rank_owner(n_serial_elements, serial_id, mesh.comm()->size()) !=
        mesh.comm()->rank()) {
      local_not_naive = 1;
      break;
    }
  }
  int global_not_naive = 0;
  MPI_Allreduce(&local_not_naive, &global_not_naive, 1, MPI_INT, MPI_MAX,
                mesh.comm()->get());
  if (require_not_naive) {
    SMESH_TEST_ASSERT(global_not_naive);
  }

  auto dist = mesh.distributed();
  SMESH_TEST_ASSERT(dist != nullptr);
  auto node_mapping = dist->node_mapping();
  SMESH_TEST_ASSERT(node_mapping != nullptr);

  LocalSideTable lst;
  lst.fill(block->element_type());
  const int nnxs = elem_num_nodes(side_type(block->element_type()));
  auto elems = block->elements()->data();
  for (ptrdiff_t i = 0; i < sideset.size(); ++i) {
    const element_idx_t e = sideset.parent()->data()[i];
    const int s = sideset.lfi()->data()[i];
    for (int ln = 0; ln < nnxs; ++ln) {
      const idx_t local_node = elems[lst(s, ln)][e];
      SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(local_node) >= 0);
      SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(local_node) < mesh.n_nodes());
      SMESH_TEST_ASSERT(node_mapping->data()[local_node] >= 0);
    }
  }

  return SMESH_TEST_SUCCESS;
}

static std::shared_ptr<Sideset> sideset_as_ptr(const Sideset &ss) {
  return Sideset::create(ss.comm(), ss.parent(), ss.lfi(), ss.block_id());
}

static int check_nodeset_and_surface_vs_serial(
    const std::shared_ptr<Mesh> &mesh, const Sideset &sideset,
    const std::vector<u8> &serial_node_mark,
    const std::vector<std::array<idx_t, 4>> &serial_faces) {
  auto dist = mesh->distributed();
  SMESH_TEST_ASSERT(dist != nullptr);
  const ptrdiff_t n_nodes_global = dist->n_nodes_global();
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(serial_node_mark.size()),
                n_nodes_global);
  auto node_mapping = dist->node_mapping();
  SMESH_TEST_ASSERT(node_mapping != nullptr);

  auto ss_ptr = sideset_as_ptr(sideset);
  auto nodeset = create_nodeset_from_sideset(mesh, ss_ptr);
  SMESH_TEST_ASSERT(nodeset != nullptr);

  std::vector<u8> local_mark(static_cast<size_t>(n_nodes_global), 0);
  for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nodeset->size()); ++i) {
    const idx_t local = nodeset->data()[i];
    SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(local) >= 0);
    SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(local) < mesh->n_nodes());
    const large_idx_t gid = node_mapping->data()[local];
    SMESH_TEST_ASSERT(gid >= 0);
    SMESH_TEST_ASSERT(gid < static_cast<large_idx_t>(n_nodes_global));
    local_mark[static_cast<size_t>(gid)] = 1;
    SMESH_TEST_ASSERT(serial_node_mark[static_cast<size_t>(gid)] == 1);
  }

  MPI_Allreduce(MPI_IN_PLACE, local_mark.data(),
                static_cast<int>(n_nodes_global), mpi_type<u8>(), MPI_MAX,
                mesh->comm()->get());
  for (ptrdiff_t i = 0; i < n_nodes_global; ++i) {
    SMESH_TEST_EQ(static_cast<int>(local_mark[static_cast<size_t>(i)]),
                  static_cast<int>(serial_node_mark[static_cast<size_t>(i)]));
  }

  auto [st, surface] = create_surface_from_sideset(mesh, ss_ptr);
  SMESH_TEST_ASSERT(surface != nullptr);
  const int nnxs = elem_num_nodes(st);
  SMESH_TEST_ASSERT(nnxs <= 4);
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(surface->extent(1)), sideset.size());

  std::vector<int> claimed(serial_faces.size(), 0);
  for (ptrdiff_t e = 0; e < static_cast<ptrdiff_t>(surface->extent(1)); ++e) {
    std::array<idx_t, 4> key = {invalid_idx<idx_t>(), invalid_idx<idx_t>(),
                                invalid_idx<idx_t>(), invalid_idx<idx_t>()};
    for (int ln = 0; ln < nnxs; ++ln) {
      const idx_t local = surface->data()[ln][e];
      key[static_cast<size_t>(ln)] =
          static_cast<idx_t>(node_mapping->data()[local]);
    }
    std::sort(key.begin(), key.begin() + nnxs);
    bool found = false;
    for (size_t j = 0; j < serial_faces.size(); ++j) {
      if (claimed[j]) {
        continue;
      }
      if (key == serial_faces[j]) {
        claimed[j] = 1;
        found = true;
        break;
      }
    }
    SMESH_TEST_ASSERT(found);
  }

  if (!claimed.empty()) {
    MPI_Allreduce(MPI_IN_PLACE, claimed.data(), static_cast<int>(claimed.size()),
                  MPI_INT, MPI_SUM, mesh->comm()->get());
    for (size_t j = 0; j < claimed.size(); ++j) {
      SMESH_TEST_EQ(claimed[j], 1);
    }
  }

  return SMESH_TEST_SUCCESS;
}

static void
serial_nodeset_and_faces(const std::shared_ptr<Mesh> &serial,
                         const std::shared_ptr<Sideset> &ss,
                         std::vector<u8> *node_mark,
                         std::vector<std::array<idx_t, 4>> *faces) {
  node_mark->assign(static_cast<size_t>(serial->n_nodes()), 0);
  auto nodeset = create_nodeset_from_sideset(serial, ss);
  for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nodeset->size()); ++i) {
    (*node_mark)[static_cast<size_t>(nodeset->data()[i])] = 1;
  }

  auto [st, surface] = create_surface_from_sideset(serial, ss);
  const int nnxs = elem_num_nodes(st);
  faces->assign(static_cast<size_t>(surface->extent(1)),
                {invalid_idx<idx_t>(), invalid_idx<idx_t>(),
                 invalid_idx<idx_t>(), invalid_idx<idx_t>()});
  for (ptrdiff_t e = 0; e < static_cast<ptrdiff_t>(surface->extent(1)); ++e) {
    for (int ln = 0; ln < nnxs; ++ln) {
      (*faces)[static_cast<size_t>(e)][static_cast<size_t>(ln)] =
          surface->data()[ln][e];
    }
    std::sort((*faces)[static_cast<size_t>(e)].begin(),
              (*faces)[static_cast<size_t>(e)].begin() + nnxs);
  }
}

static int test_mpi_sideset_redistribute_checkerboard() {
  auto comm = Communicator::world();
  if (comm->size() < 2) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr));
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mpi_multiblock_topology_%d_%d", comm->size(), token);
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(4, 2 * comm->size());
  const ptrdiff_t ny = nx;
  const ptrdiff_t nz = 2;

  ptrdiff_t serial_white = 0;
  ptrdiff_t serial_black = 0;
  ptrdiff_t serial_n_nodes = 0;
  ptrdiff_t n_serial_white_ss = 0;
  ptrdiff_t n_serial_black_ss = 0;
  std::vector<element_idx_t> serial_white_parent;
  std::vector<i16> serial_white_lfi;
  std::vector<element_idx_t> serial_black_parent;
  std::vector<i16> serial_black_lfi;
  std::vector<u8> serial_white_nodes;
  std::vector<u8> serial_black_nodes;
  std::vector<std::array<idx_t, 4>> serial_white_faces;
  std::vector<std::array<idx_t, 4>> serial_black_faces;
  ptrdiff_t n_serial_white_faces = 0;
  ptrdiff_t n_serial_black_faces = 0;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial = Mesh::create_hex8_checkerboard_cube(
        Communicator::self(), nx, ny, nz, 0, 0, 0, 1, 1, 1);
    SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
    serial_white = serial->n_elements(0);
    serial_black = serial->n_elements(1);
    serial_n_nodes = serial->n_nodes();

    auto white_ss = Sideset::create_from_plane(
        serial, 1, 0, 0, 0.0, 1e-6, {"white"});
    auto black_ss = Sideset::create_from_plane(
        serial, 0, 0, 1, 1.0, 1e-6, {"black"});
    SMESH_TEST_EQ(static_cast<int>(white_ss.size()), 1);
    SMESH_TEST_EQ(static_cast<int>(black_ss.size()), 1);
    n_serial_white_ss = white_ss[0]->size();
    n_serial_black_ss = black_ss[0]->size();
    serial_white_parent.assign(white_ss[0]->parent()->data(),
                               white_ss[0]->parent()->data() + n_serial_white_ss);
    serial_white_lfi.assign(white_ss[0]->lfi()->data(),
                            white_ss[0]->lfi()->data() + n_serial_white_ss);
    serial_black_parent.assign(black_ss[0]->parent()->data(),
                               black_ss[0]->parent()->data() + n_serial_black_ss);
    serial_black_lfi.assign(black_ss[0]->lfi()->data(),
                            black_ss[0]->lfi()->data() + n_serial_black_ss);

    serial_nodeset_and_faces(serial, white_ss[0], &serial_white_nodes,
                             &serial_white_faces);
    serial_nodeset_and_faces(serial, black_ss[0], &serial_black_nodes,
                             &serial_black_faces);
    n_serial_white_faces = static_cast<ptrdiff_t>(serial_white_faces.size());
    n_serial_black_faces = static_cast<ptrdiff_t>(serial_black_faces.size());

    std::filesystem::create_directories(
        (mesh_path / "sidesets" / "x0_white").to_string());
    std::filesystem::create_directories(
        (mesh_path / "sidesets" / "z1_black").to_string());
    SMESH_TEST_ASSERT(white_ss[0]->write(mesh_path / "sidesets" / "x0_white") ==
                      SMESH_SUCCESS);
    SMESH_TEST_ASSERT(black_ss[0]->write(mesh_path / "sidesets" / "z1_black") ==
                      SMESH_SUCCESS);
  }

  comm->broadcast(&serial_white, 1, 0);
  comm->broadcast(&serial_black, 1, 0);
  comm->broadcast(&serial_n_nodes, 1, 0);
  comm->broadcast(&n_serial_white_ss, 1, 0);
  comm->broadcast(&n_serial_black_ss, 1, 0);
  comm->broadcast(&n_serial_white_faces, 1, 0);
  comm->broadcast(&n_serial_black_faces, 1, 0);
  serial_white_parent.resize(static_cast<size_t>(n_serial_white_ss));
  serial_white_lfi.resize(static_cast<size_t>(n_serial_white_ss));
  serial_black_parent.resize(static_cast<size_t>(n_serial_black_ss));
  serial_black_lfi.resize(static_cast<size_t>(n_serial_black_ss));
  serial_white_nodes.resize(static_cast<size_t>(serial_n_nodes));
  serial_black_nodes.resize(static_cast<size_t>(serial_n_nodes));
  serial_white_faces.resize(static_cast<size_t>(n_serial_white_faces));
  serial_black_faces.resize(static_cast<size_t>(n_serial_black_faces));
  if (n_serial_white_ss > 0) {
    comm->broadcast(serial_white_parent.data(),
                    static_cast<int>(n_serial_white_ss), 0);
    comm->broadcast(serial_white_lfi.data(),
                    static_cast<int>(n_serial_white_ss), 0);
  }
  if (n_serial_black_ss > 0) {
    comm->broadcast(serial_black_parent.data(),
                    static_cast<int>(n_serial_black_ss), 0);
    comm->broadcast(serial_black_lfi.data(),
                    static_cast<int>(n_serial_black_ss), 0);
  }
  if (serial_n_nodes > 0) {
    comm->broadcast(serial_white_nodes.data(),
                    static_cast<int>(serial_n_nodes), 0);
    comm->broadcast(serial_black_nodes.data(),
                    static_cast<int>(serial_n_nodes), 0);
  }
  if (n_serial_white_faces > 0) {
    comm->broadcast(reinterpret_cast<idx_t *>(serial_white_faces.data()),
                    static_cast<int>(n_serial_white_faces * 4), 0);
  }
  if (n_serial_black_faces > 0) {
    comm->broadcast(reinterpret_cast<idx_t *>(serial_black_faces.data()),
                    static_cast<int>(n_serial_black_faces * 4), 0);
  }
  comm->barrier();

  unsetenv("SMESH_REORDER");
  auto mesh = Mesh::create_from_file(comm, mesh_path);
  SMESH_TEST_ASSERT(mesh != nullptr);

  Sideset white_ss;
  SMESH_TEST_ASSERT(white_ss.read_and_redistibute(
                        mesh, mesh_path / "sidesets" / "x0_white", 0) ==
                    SMESH_SUCCESS);
  SMESH_TEST_EQ(check_sideset_sfc_identity(*mesh, white_ss,
                                           serial_white_parent.data(),
                                           serial_white_lfi.data(),
                                           n_serial_white_ss, serial_white, 1),
                SMESH_TEST_SUCCESS);
  SMESH_TEST_EQ(check_nodeset_and_surface_vs_serial(
                    mesh, white_ss, serial_white_nodes, serial_white_faces),
                SMESH_TEST_SUCCESS);

  Sideset black_ss;
  SMESH_TEST_ASSERT(black_ss.read_and_redistibute(
                        mesh, mesh_path / "sidesets" / "z1_black", 1) ==
                    SMESH_SUCCESS);
  SMESH_TEST_EQ(check_sideset_sfc_identity(*mesh, black_ss,
                                           serial_black_parent.data(),
                                           serial_black_lfi.data(),
                                           n_serial_black_ss, serial_black, 1),
                SMESH_TEST_SUCCESS);
  SMESH_TEST_EQ(check_nodeset_and_surface_vs_serial(
                    mesh, black_ss, serial_black_nodes, serial_black_faces),
                SMESH_TEST_SUCCESS);

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  comm->barrier();
  return SMESH_TEST_SUCCESS;
}

static int test_mpi_skin_sidesets_checkerboard() {
  auto comm = Communicator::world();
  if (comm->size() < 2) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr)) + 17;
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mpi_multiblock_skin_%d_%d", comm->size(), token);
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(4, 2 * comm->size());
  const ptrdiff_t ny = nx;
  const ptrdiff_t nz = 2;

  ptrdiff_t serial_n_elements[2] = {0, 0};
  ptrdiff_t n_serial_ss[2] = {0, 0};
  std::vector<element_idx_t> serial_parent[2];
  std::vector<i16> serial_lfi[2];

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial = Mesh::create_hex8_checkerboard_cube(
        Communicator::self(), nx, ny, nz, 0, 0, 0, 1, 1, 1);
    SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
    serial_n_elements[0] = serial->n_elements(0);
    serial_n_elements[1] = serial->n_elements(1);

    auto skins = skin_sidesets(serial);
    SMESH_TEST_EQ(static_cast<int>(skins.size()), 2);
    for (int b = 0; b < 2; ++b) {
      n_serial_ss[b] = skins[b]->size();
      serial_parent[b].assign(skins[b]->parent()->data(),
                              skins[b]->parent()->data() + n_serial_ss[b]);
      serial_lfi[b].assign(skins[b]->lfi()->data(),
                           skins[b]->lfi()->data() + n_serial_ss[b]);
    }
  }

  comm->broadcast(serial_n_elements, 2, 0);
  comm->broadcast(n_serial_ss, 2, 0);
  for (int b = 0; b < 2; ++b) {
    serial_parent[b].resize(static_cast<size_t>(n_serial_ss[b]));
    serial_lfi[b].resize(static_cast<size_t>(n_serial_ss[b]));
    if (n_serial_ss[b] > 0) {
      comm->broadcast(serial_parent[b].data(),
                      static_cast<int>(n_serial_ss[b]), 0);
      comm->broadcast(serial_lfi[b].data(), static_cast<int>(n_serial_ss[b]),
                      0);
    }
  }
  comm->barrier();

  unsetenv("SMESH_REORDER");
  auto mesh = Mesh::create_from_file(comm, mesh_path);
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto skins = skin_sidesets(mesh);
  SMESH_TEST_EQ(static_cast<int>(skins.size()), 2);
  for (int b = 0; b < 2; ++b) {
    SMESH_TEST_ASSERT(skins[b] != nullptr);
    SMESH_TEST_EQ(skins[b]->block_id(), static_cast<block_idx_t>(b));
    SMESH_TEST_EQ(check_sideset_sfc_identity(*mesh, *skins[b],
                                             serial_parent[b].data(),
                                             serial_lfi[b].data(),
                                             n_serial_ss[b],
                                             serial_n_elements[b], 1),
                  SMESH_TEST_SUCCESS);
  }

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  comm->barrier();
  return SMESH_TEST_SUCCESS;
}

static int test_mpi_single_block_sideset_sfc() {
  auto comm = Communicator::world();
  if (comm->size() < 2) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr)) + 23;
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mpi_single_block_ss_%d_%d", comm->size(), token);
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(4, 2 * comm->size());
  const ptrdiff_t ny = nx;
  const ptrdiff_t nz = 2;

  ptrdiff_t serial_n_elements = 0;
  ptrdiff_t n_serial_ss = 0;
  std::vector<element_idx_t> serial_parent;
  std::vector<i16> serial_lfi;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
    SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
    serial_n_elements = serial->n_elements();
    auto ss = Sideset::create_from_plane(serial, 1, 0, 0, 0.0, 1e-6);
    SMESH_TEST_EQ(static_cast<int>(ss.size()), 1);
    n_serial_ss = ss[0]->size();
    serial_parent.assign(ss[0]->parent()->data(),
                         ss[0]->parent()->data() + n_serial_ss);
    serial_lfi.assign(ss[0]->lfi()->data(), ss[0]->lfi()->data() + n_serial_ss);
    std::filesystem::create_directories(
        (mesh_path / "sidesets" / "x0").to_string());
    SMESH_TEST_ASSERT(ss[0]->write(mesh_path / "sidesets" / "x0") ==
                      SMESH_SUCCESS);
  }

  comm->broadcast(&serial_n_elements, 1, 0);
  comm->broadcast(&n_serial_ss, 1, 0);
  serial_parent.resize(static_cast<size_t>(n_serial_ss));
  serial_lfi.resize(static_cast<size_t>(n_serial_ss));
  if (n_serial_ss > 0) {
    comm->broadcast(serial_parent.data(), static_cast<int>(n_serial_ss), 0);
    comm->broadcast(serial_lfi.data(), static_cast<int>(n_serial_ss), 0);
  }
  comm->barrier();

  unsetenv("SMESH_REORDER");
  auto mesh = Mesh::create_from_file(comm, mesh_path);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(static_cast<int>(mesh->n_blocks()), 1);

  Sideset ss;
  SMESH_TEST_ASSERT(
      ss.read_and_redistibute(mesh, mesh_path / "sidesets" / "x0", 0) ==
      SMESH_SUCCESS);
  SMESH_TEST_EQ(check_sideset_sfc_identity(*mesh, ss, serial_parent.data(),
                                           serial_lfi.data(), n_serial_ss,
                                           serial_n_elements, 0),
                SMESH_TEST_SUCCESS);

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  comm->barrier();
  return SMESH_TEST_SUCCESS;
}
#endif

int main(int argc, char **argv) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_serial_checkerboard_graphs);
  SMESH_RUN_TEST(test_serial_hex8_tet4_graphs);
  SMESH_RUN_TEST(test_serial_transforms);
  SMESH_RUN_TEST(test_serial_promote_refine_vs_single_block);
  SMESH_RUN_TEST(test_serial_hex8_multiblock_refine);
  SMESH_RUN_TEST(test_serial_quad4_refine);
  SMESH_RUN_TEST(test_serial_quad4_multiblock_refine);
  SMESH_RUN_TEST(test_serial_quadshell4_refine);
  SMESH_RUN_TEST(test_serial_wedge6_refine);
  SMESH_RUN_TEST(test_serial_wedge6_multiblock_refine);
  SMESH_RUN_TEST(test_serial_trishell3_refine);
  SMESH_RUN_TEST(test_serial_trishell3_multiblock_refine);
  SMESH_RUN_TEST(test_serial_edge2_refine);
  SMESH_RUN_TEST(test_serial_edgeshell2_refine);
  SMESH_RUN_TEST(test_serial_edge2_multiblock_refine);
  SMESH_RUN_TEST(test_serial_extrude_vs_single_block);
  SMESH_RUN_TEST(test_serial_hex_dominant);
#ifdef SMESH_ENABLE_MPI
  SMESH_RUN_TEST(test_mpi_sideset_redistribute_checkerboard);
  SMESH_RUN_TEST(test_mpi_skin_sidesets_checkerboard);
  SMESH_RUN_TEST(test_mpi_single_block_sideset_sfc);
#endif
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
