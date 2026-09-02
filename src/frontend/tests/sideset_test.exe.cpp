#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "smesh_adjacency.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_device_sideset.hpp"
#include "smesh_mask.hpp"
#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sidesets.impl.hpp"
#include "smesh_test.hpp"

using namespace smesh;

namespace {

constexpr ptrdiff_t kNx = 2;
constexpr ptrdiff_t kNy = 3;
constexpr ptrdiff_t kNz = 4;

std::shared_ptr<Mesh> make_test_mesh() {
  return Mesh::create_hex8_cube(Communicator::self(), kNx, kNy, kNz);
}

std::shared_ptr<Sideset>
make_left_boundary_sideset(const std::shared_ptr<Mesh> &mesh) {
  auto sidesets = Sideset::create_from_selector(
      mesh,
      [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });

  if (sidesets.size() != 1) {
    return nullptr;
  }

  return sidesets[0];
}

std::vector<std::pair<element_idx_t, i16>>
sorted_side_keys(const std::shared_ptr<Sideset> &sideset) {
  std::vector<std::pair<element_idx_t, i16>> keys(sideset->size());
  const auto *parent = sideset->parent()->data();
  const auto *lfi = sideset->lfi()->data();

  for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
    keys[static_cast<size_t>(i)] = std::make_pair(parent[i], lfi[i]);
  }

  std::sort(keys.begin(), keys.end());
  return keys;
}

geom_t side_centroid_x(const std::shared_ptr<Mesh> &mesh,
                       const std::shared_ptr<Sideset> &sideset,
                       const ptrdiff_t side) {
  LocalSideTable lst;
  lst.fill(mesh->element_type(0));

  const auto *const *elements = mesh->elements(0)->data();
  const auto *const *points = mesh->points()->data();
  const element_idx_t e = sideset->parent()->data()[side];
  const i16 s = sideset->lfi()->data()[side];

  geom_t x = 0;
  const int nn = lst.nnxs_side[s];
  for (int i = 0; i < nn; ++i) {
    x += points[0][elements[lst(s, i)][e]];
  }

  return x / static_cast<geom_t>(nn);
}

using NodeKey = std::array<long long, 3>;
using FaceKey = std::vector<NodeKey>;

NodeKey quantize_node(const geom_t *const *pts, const idx_t node) {
  return {std::llround(static_cast<double>(pts[0][node]) * 1e9),
          std::llround(static_cast<double>(pts[1][node]) * 1e9),
          std::llround(static_cast<double>(pts[2][node]) * 1e9)};
}

int parse_sideset_meta_file(const Path &folder, ptrdiff_t *size_out, block_idx_t *block_id_out) {
  std::ifstream ifs((folder / "meta.yaml").to_string());
  if (!ifs.good()) {
    return SMESH_FAILURE;
  }
  bool has_size = false;
  bool has_block_id = false;
  std::string line;
  while (std::getline(ifs, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line.resize(hash);
    }
    const auto start = line.find_first_not_of(" \t");
    if (start == std::string::npos) {
      continue;
    }
    line = line.substr(start);
    if (line.compare(0, 5, "size:") == 0) {
      has_size = true;
      *size_out = static_cast<ptrdiff_t>(std::strtoll(line.c_str() + 5, nullptr, 10));
    } else if (line.compare(0, 9, "block_id:") == 0) {
      has_block_id = true;
      *block_id_out = static_cast<block_idx_t>(std::strtol(line.c_str() + 9, nullptr, 10));
    }
  }
  return (has_size && has_block_id) ? SMESH_SUCCESS : SMESH_FAILURE;
}

std::vector<FaceKey> corner_face_keys(const std::shared_ptr<Mesh> &mesh,
                                      const std::shared_ptr<Sideset> &sideset) {
  std::vector<FaceKey> faces;
  auto block = mesh->block(sideset->block_id());
  const enum ElemType et = block->element_type();
  const bool is_ss = is_semistructured_type(et);
  const enum ElemType family = is_ss ? ss_source_family(et) : et;
  int corners[8] = {};
  int n_corners = 0;
  if (is_ss) {
    ss_source_family_corners(family, semistructured_level(et), corners, &n_corners);
  }
  LocalSideTable lst;
  lst.fill(family);
  auto elems = block->elements()->data();
  auto pts = mesh->points()->data();
  const ptrdiff_t n_e = block->n_elements();
  for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
    const element_idx_t e = sideset->parent()->data()[i];
    const i16 s = sideset->lfi()->data()[i];
    if (e < 0 || e >= n_e) {
      return {};
    }
    FaceKey key(static_cast<size_t>(lst.nnxs_side[s]));
    for (int ln = 0; ln < lst.nnxs_side[s]; ++ln) {
      const int soa_row = is_ss ? corners[lst(s, ln)] : lst(s, ln);
      key[static_cast<size_t>(ln)] = quantize_node(pts, elems[soa_row][e]);
    }
    std::sort(key.begin(), key.end());
    faces.push_back(std::move(key));
  }
  std::sort(faces.begin(), faces.end());
  return faces;
}

} // namespace

int test_sideset_creation() {
  auto mesh = make_test_mesh();
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(sideset != nullptr);
  SMESH_TEST_EQ(sideset->block_id(), static_cast<block_idx_t>(0));
  SMESH_TEST_EQ(sideset->size(), kNy * kNz);

  return SMESH_TEST_SUCCESS;
}

int test_sideset_to_nodeset_conversion() {
  auto mesh = make_test_mesh();
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(sideset != nullptr);
  auto nodeset = create_nodeset_from_sideset(mesh, sideset);
  SMESH_TEST_ASSERT(nodeset != nullptr);
  SMESH_TEST_EQ(nodeset->size(), static_cast<size_t>((kNy + 1) * (kNz + 1)));

  auto points = mesh->points()->data();
  for (size_t i = 0; i < nodeset->size(); ++i) {
    const idx_t node = nodeset->data()[i];
    SMESH_TEST_ASSERT(node >= 0);
    SMESH_TEST_ASSERT(node < mesh->n_nodes());
    SMESH_TEST_APPROXEQ(points[0][node], 0.0, 1e-12);
  }

  return SMESH_TEST_SUCCESS;
}

int test_sideset_io_write_read_identity() {
  auto mesh = make_test_mesh();
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(sideset != nullptr);

  const auto token = static_cast<long long>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_sideset_test_%lld", token);
  const Path path(path_buffer);

  std::filesystem::remove_all(path.to_string());

  SMESH_TEST_ASSERT(sideset->write(path) == SMESH_SUCCESS);

  ptrdiff_t meta_size = -1;
  block_idx_t meta_block_id = static_cast<block_idx_t>(-1);
  SMESH_TEST_ASSERT(parse_sideset_meta_file(path, &meta_size, &meta_block_id) == SMESH_SUCCESS);
  SMESH_TEST_EQ(meta_size, sideset->size());
  SMESH_TEST_EQ(meta_block_id, sideset->block_id());

  auto read_back = Sideset::create_from_file(Communicator::self(), path);
  SMESH_TEST_ASSERT(read_back != nullptr);
  SMESH_TEST_EQ(read_back->block_id(), sideset->block_id());
  SMESH_TEST_EQ(read_back->size(), sideset->size());

  auto read_override = std::make_shared<Sideset>();
  SMESH_TEST_ASSERT(read_override->read(Communicator::self(), path, 99) == SMESH_SUCCESS);
  SMESH_TEST_EQ(read_override->block_id(), sideset->block_id());

  for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
    SMESH_TEST_EQ(read_back->parent()->data()[i], sideset->parent()->data()[i]);
    SMESH_TEST_EQ(read_back->lfi()->data()[i], sideset->lfi()->data()[i]);
  }

  std::filesystem::remove_all(path.to_string());
  return SMESH_TEST_SUCCESS;
}

int test_sideset_to_device_preserves_block_and_mapping() {
  auto parent = create_host_buffer<element_idx_t>(2);
  auto lfi = create_host_buffer<i16>(2);
  auto mapping = create_host_buffer<large_idx_t>(4);
  SMESH_TEST_ASSERT(parent != nullptr);
  SMESH_TEST_ASSERT(lfi != nullptr);
  SMESH_TEST_ASSERT(mapping != nullptr);
  parent->data()[0] = 1;
  parent->data()[1] = 3;
  lfi->data()[0] = 0;
  lfi->data()[1] = 2;
  mapping->data()[0] = 10;
  mapping->data()[1] = 20;
  mapping->data()[2] = 30;
  mapping->data()[3] = 40;

  auto host = Sideset::create(Communicator::self(), parent, lfi, 2, mapping);
  SMESH_TEST_ASSERT(host != nullptr);

  auto device_ss = to_device(host);
  SMESH_TEST_ASSERT(device_ss != nullptr);
  SMESH_TEST_EQ(device_ss->block_id(), static_cast<block_idx_t>(2));
  SMESH_TEST_EQ(device_ss->size(), static_cast<ptrdiff_t>(2));
  SMESH_TEST_ASSERT(device_ss->element_mapping() != nullptr);
  SMESH_TEST_ASSERT(device_ss->element_mapping().get() == host->element_mapping().get());
  SMESH_TEST_EQ(device_ss->element_mapping()->data()[1], static_cast<large_idx_t>(20));

  auto host_parent = to_host(device_ss->parent());
  auto host_lfi = to_host(device_ss->lfi());
  SMESH_TEST_ASSERT(host_parent != nullptr);
  SMESH_TEST_ASSERT(host_lfi != nullptr);
  SMESH_TEST_EQ(host_parent->data()[0], static_cast<element_idx_t>(1));
  SMESH_TEST_EQ(host_parent->data()[1], static_cast<element_idx_t>(3));
  SMESH_TEST_EQ(host_lfi->data()[0], static_cast<i16>(0));
  SMESH_TEST_EQ(host_lfi->data()[1], static_cast<i16>(2));

  return SMESH_TEST_SUCCESS;
}

static int check_ss_level_invariance(const std::shared_ptr<Mesh> &mesh,
                                     const std::shared_ptr<Sideset> &sideset,
                                     const int fine_level,
                                     const int coarse_level) {
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  SMESH_TEST_ASSERT(sideset->size() > 0);

  auto ss_fine = to_semistructured(fine_level, mesh);
  auto ss_coarse = to_semistructured(coarse_level, mesh);
  SMESH_TEST_ASSERT(ss_fine != nullptr);
  SMESH_TEST_ASSERT(ss_coarse != nullptr);
  SMESH_TEST_EQ(ss_fine->n_elements(sideset->block_id()), mesh->n_elements(sideset->block_id()));
  SMESH_TEST_EQ(ss_coarse->n_elements(sideset->block_id()), mesh->n_elements(sideset->block_id()));

  auto [fine_st, fine_surf] = create_surface_from_sideset(ss_fine, sideset);
  auto [coarse_st, coarse_surf] = create_surface_from_sideset(ss_coarse, sideset);
  SMESH_TEST_ASSERT(fine_surf != nullptr);
  SMESH_TEST_ASSERT(coarse_surf != nullptr);
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(fine_surf->extent(1)), sideset->size());
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(coarse_surf->extent(1)), sideset->size());
  SMESH_TEST_ASSERT(fine_st != INVALID);
  SMESH_TEST_ASSERT(coarse_st != INVALID);

  const auto coarse_keys = corner_face_keys(ss_coarse, sideset);
  const auto fine_keys = corner_face_keys(ss_fine, sideset);
  SMESH_TEST_EQ(coarse_keys.size(), static_cast<size_t>(sideset->size()));
  SMESH_TEST_EQ(fine_keys.size(), coarse_keys.size());
  for (size_t i = 0; i < coarse_keys.size(); ++i) {
    SMESH_TEST_EQ(fine_keys[i].size(), coarse_keys[i].size());
    for (size_t k = 0; k < coarse_keys[i].size(); ++k) {
      SMESH_TEST_EQ(fine_keys[i][k][0], coarse_keys[i][k][0]);
      SMESH_TEST_EQ(fine_keys[i][k][1], coarse_keys[i][k][1]);
      SMESH_TEST_EQ(fine_keys[i][k][2], coarse_keys[i][k][2]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

int test_sshex_sideset_level_invariance() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  return check_ss_level_invariance(mesh, sideset, 4, 2);
}

int test_sstet_sideset_level_invariance() {
  auto mesh = Mesh::create_tet4_cube(Communicator::self(), 2, 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  auto sidesets = Sideset::create_from_selector(
      mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(sidesets.size(), static_cast<size_t>(1));
  return check_ss_level_invariance(mesh, sidesets[0], 3, 2);
}

int test_sideset_select_propagate_cube_mesh() {
  auto mesh = make_test_mesh();
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto skin = skin_sideset(mesh);
  auto left = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(skin != nullptr);
  SMESH_TEST_ASSERT(left != nullptr);

  const auto expected = sorted_side_keys(left);
  SMESH_TEST_EQ(expected.size(), static_cast<size_t>(kNy * kNz));

  const element_idx_t seed_parent = left->parent()->data()[0];

  ptrdiff_t seed_side = -1;
  for (ptrdiff_t i = 0; i < skin->size(); ++i) {
    if (skin->parent()->data()[i] == seed_parent &&
        skin->lfi()->data()[i] == left->lfi()->data()[0]) {
      seed_side = i;
      break;
    }
  }

  SMESH_TEST_ASSERT(seed_side >= 0);
  SMESH_TEST_APPROXEQ(side_centroid_x(mesh, skin, seed_side), 0.0, 1e-12);

  auto n2e = mesh->node_to_element_graph();
  SMESH_TEST_ASSERT(n2e != nullptr);

  auto selected = create_host_buffer<mask_t>(mask_count(skin->size()));
  SMESH_TEST_ASSERT(selected != nullptr);

  const int err = sideset_select_propagate(
      skin->size(), skin->parent()->data(), skin->lfi()->data(),
      n2e->rowptr()->data(), n2e->colidx()->data(), mesh->element_type(0),
      mesh->n_elements(), mesh->elements(0)->data(),
      static_cast<element_idx_t>(seed_side), selected->data(),
      [&mesh, &skin](const ptrdiff_t, const ptrdiff_t next_side) {
        return std::abs(side_centroid_x(mesh, skin, next_side)) < 1e-12;
      });

  SMESH_TEST_EQ(err, SMESH_SUCCESS);

  std::vector<std::pair<element_idx_t, i16>> actual;
  actual.reserve(expected.size());
  for (ptrdiff_t i = 0; i < skin->size(); ++i) {
    if (!mask_get(i, selected->data())) {
      continue;
    }

    actual.emplace_back(skin->parent()->data()[i], skin->lfi()->data()[i]);
  }

  std::sort(actual.begin(), actual.end());
  SMESH_TEST_EQ(actual.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    SMESH_TEST_EQ(actual[i].first, expected[i].first);
    SMESH_TEST_EQ(actual[i].second, expected[i].second);
  }

  return SMESH_TEST_SUCCESS;
}

int test_hex27_element_contract() {
  SMESH_TEST_ASSERT(HEX27 != PROTEUS_HEX27);
  SMESH_TEST_EQ(type_from_string("HEX27"), HEX27);
  SMESH_TEST_EQ(type_from_string("PROTEUS_HEX27"), PROTEUS_HEX27);
  SMESH_TEST_EQ(elem_num_nodes(HEX27), 27);
  SMESH_TEST_EQ(elem_num_sides(HEX27), 6);
  SMESH_TEST_EQ(side_type(HEX27), QUAD9);
  SMESH_TEST_EQ(shell_type(side_type(HEX27)), QUADSHELL9);
  SMESH_TEST_EQ(side_type(EDGE3), NODE1);
  SMESH_TEST_EQ(elem_num_sides(EDGE3), 2);
  SMESH_TEST_EQ(elem_manifold_dim(EDGE3), 1);
  SMESH_TEST_EQ(type_from_string("EDGESHELL2"), EDGESHELL2);
  SMESH_TEST_EQ(type_from_string("EDGESHELL3"), EDGESHELL3);
  SMESH_TEST_ASSERT(std::strcmp(type_to_string(EDGESHELL3), "EDGESHELL3") == 0);
  SMESH_TEST_EQ(shell_type(EDGE3), EDGESHELL3);
  SMESH_TEST_EQ(elem_num_nodes(EDGESHELL3), 3);
  SMESH_TEST_EQ(elem_num_sides(EDGESHELL3), 2);
  SMESH_TEST_EQ(elem_manifold_dim(EDGESHELL3), 1);
  SMESH_TEST_EQ(shell_type(side_type(QUAD9)), EDGESHELL3);
  SMESH_TEST_EQ(shell_type(side_type(TRI6)), EDGESHELL3);
  SMESH_TEST_ASSERT(!is_semistructured_type(HEX27));
  SMESH_TEST_ASSERT(is_semistructured_type(PROTEUS_HEX27));

  LocalSideTable table;
  table.fill(HEX27);
  SMESH_TEST_EQ(table.nnxs, 9);

  const int expected[6][9] = {
      {0, 1, 5, 4, 8, 17, 12, 16, 20},
      {1, 2, 6, 5, 9, 18, 13, 17, 21},
      {2, 3, 7, 6, 10, 19, 14, 18, 22},
      {3, 0, 4, 7, 11, 16, 15, 19, 23},
      {3, 2, 1, 0, 10, 9, 8, 11, 24},
      {4, 5, 6, 7, 12, 13, 14, 15, 25},
  };
  for (int side = 0; side < 6; ++side) {
    for (int node = 0; node < 9; ++node) {
      SMESH_TEST_EQ(table(side, node), expected[side][node]);
    }
  }

  return SMESH_TEST_SUCCESS;
}

int test_hex27_cube_uses_conventional_ordering() {
  auto mesh = Mesh::create_cube(Communicator::self(), HEX27, 1, 1, 1, 0, 0, 0,
                                1, 1, 1);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_EQ(mesh->element_type(0), HEX27);
  SMESH_TEST_EQ(mesh->n_elements(0), static_cast<ptrdiff_t>(1));
  SMESH_TEST_EQ(mesh->n_nodes(), static_cast<ptrdiff_t>(27));

  const geom_t expected[27][3] = {
      {0, 0, 0},       {1, 0, 0},       {1, 1, 0},       {0, 1, 0},
      {0, 0, 1},       {1, 0, 1},       {1, 1, 1},       {0, 1, 1},
      {0.5, 0, 0},     {1, 0.5, 0},     {0.5, 1, 0},     {0, 0.5, 0},
      {0.5, 0, 1},     {1, 0.5, 1},     {0.5, 1, 1},     {0, 0.5, 1},
      {0, 0, 0.5},     {1, 0, 0.5},     {1, 1, 0.5},     {0, 1, 0.5},
      {0.5, 0, 0.5},   {1, 0.5, 0.5},   {0.5, 1, 0.5},   {0, 0.5, 0.5},
      {0.5, 0.5, 0},   {0.5, 0.5, 1},   {0.5, 0.5, 0.5},
  };
  auto elements = mesh->elements(0)->data();
  auto points = mesh->points()->data();
  for (int node = 0; node < 27; ++node) {
    const idx_t global = elements[node][0];
    for (int component = 0; component < 3; ++component) {
      SMESH_TEST_APPROXEQ(points[component][global],
                          expected[node][component], 1e-12);
    }
  }

  return SMESH_TEST_SUCCESS;
}

int test_upper_triangular_graph_with_unused_node() {
  auto elements = create_host_buffer<idx_t>(3, 1);
  auto points = create_host_buffer<geom_t>(3, 4);
  SMESH_TEST_ASSERT(elements != nullptr);
  SMESH_TEST_ASSERT(points != nullptr);

  auto elems = elements->data();
  elems[0][0] = 0;
  elems[1][0] = 1;
  elems[2][0] = 2;

  auto mesh = std::make_shared<Mesh>(Communicator::self(), TRI3, elements, points);
  auto graph = mesh->node_to_node_graph_upper_triangular();
  SMESH_TEST_ASSERT(graph != nullptr);

  const count_t expected_rowptr[5] = {0, 2, 3, 3, 3};
  const idx_t expected_colidx[3] = {1, 2, 2};

  SMESH_TEST_EQ(graph->rowptr()->size(), static_cast<size_t>(5));
  SMESH_TEST_EQ(graph->colidx()->size(), static_cast<size_t>(3));

  for (size_t i = 0; i < 5; ++i) {
    SMESH_TEST_EQ(graph->rowptr()->data()[i], expected_rowptr[i]);
  }

  for (size_t i = 0; i < 3; ++i) {
    SMESH_TEST_EQ(graph->colidx()->data()[i], expected_colidx[i]);
  }

  return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);

  SMESH_RUN_TEST(test_sideset_creation);
  SMESH_RUN_TEST(test_sideset_to_nodeset_conversion);
  SMESH_RUN_TEST(test_sideset_io_write_read_identity);
  SMESH_RUN_TEST(test_sideset_to_device_preserves_block_and_mapping);
  SMESH_RUN_TEST(test_sshex_sideset_level_invariance);
  SMESH_RUN_TEST(test_sstet_sideset_level_invariance);
  SMESH_RUN_TEST(test_sideset_select_propagate_cube_mesh);
  SMESH_RUN_TEST(test_hex27_element_contract);
  SMESH_RUN_TEST(test_hex27_cube_uses_conventional_ordering);
  SMESH_RUN_TEST(test_upper_triangular_graph_with_unused_node);

  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
