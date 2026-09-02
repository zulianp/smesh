#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "smesh_adjacency.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_device_sideset.hpp"
#include "smesh_mask.hpp"
#include "smesh_mesh.hpp"
#include "smesh_mesh_reorder.hpp"
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

NodeKey quantize_node(const geom_t *const *pts, const idx_t node, const int sdim) {
  return {std::llround(static_cast<double>(pts[0][node]) * 1e9),
          sdim > 1 ? std::llround(static_cast<double>(pts[1][node]) * 1e9) : 0,
          sdim > 2 ? std::llround(static_cast<double>(pts[2][node]) * 1e9) : 0};
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
  const int sdim = mesh->spatial_dimension();
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
      key[static_cast<size_t>(ln)] = quantize_node(pts, elems[soa_row][e], sdim);
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

int test_sideset_remap_on_sfc_reorder() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  const auto keys_before = corner_face_keys(mesh, sideset);
  SMESH_TEST_ASSERT(!keys_before.empty());

  std::vector<std::shared_ptr<Sideset>> ss = {sideset};
  SMESH_TEST_ASSERT(SFC("random3").reorder(*mesh, ss) == SMESH_SUCCESS);

  const auto keys_mapped = corner_face_keys(mesh, sideset);
  auto recreated = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(recreated != nullptr);
  const auto keys_recreated = corner_face_keys(mesh, recreated);
  SMESH_TEST_EQ(keys_mapped.size(), keys_before.size());
  SMESH_TEST_EQ(keys_recreated.size(), keys_before.size());
  for (size_t i = 0; i < keys_before.size(); ++i) {
    SMESH_TEST_EQ(keys_mapped[i].size(), keys_before[i].size());
    SMESH_TEST_EQ(keys_recreated[i].size(), keys_before[i].size());
    for (size_t k = 0; k < keys_before[i].size(); ++k) {
      SMESH_TEST_EQ(keys_mapped[i][k][0], keys_before[i][k][0]);
      SMESH_TEST_EQ(keys_mapped[i][k][1], keys_before[i][k][1]);
      SMESH_TEST_EQ(keys_mapped[i][k][2], keys_before[i][k][2]);
      SMESH_TEST_EQ(keys_recreated[i][k][0], keys_mapped[i][k][0]);
      SMESH_TEST_EQ(keys_recreated[i][k][1], keys_mapped[i][k][1]);
      SMESH_TEST_EQ(keys_recreated[i][k][2], keys_mapped[i][k][2]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

int test_sideset_remap_from_tags() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  const auto keys_before = corner_face_keys(mesh, sideset);

  auto tags = create_host_buffer<idx_t>(mesh->n_elements(0));
  for (ptrdiff_t i = 0; i < mesh->n_elements(0); ++i) {
    tags->data()[i] = static_cast<idx_t>(i & 1);
  }
  std::vector<std::shared_ptr<Sideset>> ss = {sideset};
  mesh->reorder_elements_from_tags(0, tags, ss);

  const auto keys_mapped = corner_face_keys(mesh, sideset);
  SMESH_TEST_EQ(keys_mapped.size(), keys_before.size());
  for (size_t i = 0; i < keys_before.size(); ++i) {
    for (size_t k = 0; k < keys_before[i].size(); ++k) {
      SMESH_TEST_EQ(keys_mapped[i][k][0], keys_before[i][k][0]);
      SMESH_TEST_EQ(keys_mapped[i][k][1], keys_before[i][k][1]);
      SMESH_TEST_EQ(keys_mapped[i][k][2], keys_before[i][k][2]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

static int check_refine_sideset_map(const std::shared_ptr<Mesh> &mesh,
                                    const std::shared_ptr<Sideset> &sideset,
                                    const ptrdiff_t expand) {
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  auto fine = refine(mesh, 1);
  SMESH_TEST_ASSERT(fine != nullptr);
  auto mapped = map_sideset_through_refine(mesh, sideset, fine);
  SMESH_TEST_ASSERT(mapped != nullptr);
  SMESH_TEST_EQ(mapped->size(), sideset->size() * expand);

  auto recreated = Sideset::create_from_selector(
      fine, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(recreated.size(), static_cast<size_t>(1));
  const auto keys_mapped = corner_face_keys(fine, mapped);
  const auto keys_recreated = corner_face_keys(fine, recreated[0]);
  SMESH_TEST_EQ(keys_mapped.size(), keys_recreated.size());
  for (size_t i = 0; i < keys_mapped.size(); ++i) {
    SMESH_TEST_EQ(keys_mapped[i].size(), keys_recreated[i].size());
    for (size_t k = 0; k < keys_mapped[i].size(); ++k) {
      SMESH_TEST_EQ(keys_mapped[i][k][0], keys_recreated[i][k][0]);
      SMESH_TEST_EQ(keys_mapped[i][k][1], keys_recreated[i][k][1]);
      SMESH_TEST_EQ(keys_mapped[i][k][2], keys_recreated[i][k][2]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

int test_sideset_map_hex_refine() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  return check_refine_sideset_map(mesh, sideset, 4);
}

int test_sideset_map_tet_refine() {
  auto mesh = Mesh::create_tet4_cube(Communicator::self(), 2, 2, 2);
  auto sidesets = Sideset::create_from_selector(
      mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(sidesets.size(), static_cast<size_t>(1));
  return check_refine_sideset_map(mesh, sidesets[0], 4);
}

int test_sideset_map_tri_refine() {
  auto mesh = Mesh::create_tri3_square(Communicator::self(), 2, 2);
  auto sidesets = Sideset::create_from_selector(
      mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(sidesets.size(), static_cast<size_t>(1));
  return check_refine_sideset_map(mesh, sidesets[0], 2);
}

int test_sideset_map_refine_unsupported() {
  auto hex = make_test_mesh();
  auto tet = Mesh::create_tet4_cube(Communicator::self(), 1, 1, 1);
  auto sideset = make_left_boundary_sideset(hex);
  SMESH_TEST_ASSERT(sideset != nullptr);
  auto mapped = map_sideset_through_refine(hex, sideset, tet);
  SMESH_TEST_ASSERT(mapped == nullptr);
  return SMESH_TEST_SUCCESS;
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

int test_local_side_table_higher_order_exodus() {
  LocalSideTable tet10;
  SMESH_TEST_EQ(tet10.fill(TET10), SMESH_SUCCESS);
  SMESH_TEST_EQ(tet10.nnxs, 6);
  const int tet10_expected[4][6] = {
      {0, 1, 3, 4, 8, 7},
      {1, 2, 3, 5, 9, 8},
      {0, 3, 2, 7, 9, 6},
      {0, 2, 1, 6, 5, 4},
  };
  for (int s = 0; s < 4; ++s) {
    SMESH_TEST_EQ(tet10.nnxs_side[s], 6);
    for (int n = 0; n < 6; ++n) {
      SMESH_TEST_EQ(tet10(s, n), tet10_expected[s][n]);
    }
  }

  LocalSideTable tri6;
  SMESH_TEST_EQ(tri6.fill(TRI6), SMESH_SUCCESS);
  SMESH_TEST_EQ(tri6.nnxs, 3);
  const int tri6_expected[3][3] = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
  for (int s = 0; s < 3; ++s) {
    for (int n = 0; n < 3; ++n) {
      SMESH_TEST_EQ(tri6(s, n), tri6_expected[s][n]);
    }
  }

  LocalSideTable quad9;
  SMESH_TEST_EQ(quad9.fill(QUAD9), SMESH_SUCCESS);
  SMESH_TEST_EQ(quad9.nnxs, 3);
  const int quad9_expected[4][3] = {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
  for (int s = 0; s < 4; ++s) {
    for (int n = 0; n < 3; ++n) {
      SMESH_TEST_EQ(quad9(s, n), quad9_expected[s][n]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

int test_local_side_table_edge_and_shell_aliases() {
  LocalSideTable edge2;
  SMESH_TEST_EQ(edge2.fill(EDGE2), SMESH_SUCCESS);
  SMESH_TEST_EQ(edge2.nnxs, 1);
  SMESH_TEST_EQ(edge2(0, 0), 0);
  SMESH_TEST_EQ(edge2(1, 0), 1);

  LocalSideTable edgeshell3;
  SMESH_TEST_EQ(edgeshell3.fill(EDGESHELL3), SMESH_SUCCESS);
  SMESH_TEST_EQ(edgeshell3.nnxs, 1);
  SMESH_TEST_EQ(edgeshell3(0, 0), 0);
  SMESH_TEST_EQ(edgeshell3(1, 0), 1);

  LocalSideTable qshell;
  SMESH_TEST_EQ(qshell.fill(QUADSHELL4), SMESH_SUCCESS);
  LocalSideTable quad4;
  SMESH_TEST_EQ(quad4.fill(QUAD4), SMESH_SUCCESS);
  SMESH_TEST_EQ(qshell.nnxs, quad4.nnxs);
  for (int s = 0; s < 4; ++s) {
    SMESH_TEST_EQ(qshell.nnxs_side[s], quad4.nnxs_side[s]);
    for (int n = 0; n < quad4.nnxs_side[s]; ++n) {
      SMESH_TEST_EQ(qshell(s, n), quad4(s, n));
    }
  }

  LocalSideTable tshell;
  SMESH_TEST_EQ(tshell.fill(TRISHELL3), SMESH_SUCCESS);
  LocalSideTable tri3;
  SMESH_TEST_EQ(tri3.fill(TRI3), SMESH_SUCCESS);
  for (int s = 0; s < 3; ++s) {
    for (int n = 0; n < 2; ++n) {
      SMESH_TEST_EQ(tshell(s, n), tri3(s, n));
    }
  }

  LocalSideTable bad;
  SMESH_TEST_EQ(bad.fill(TET20), SMESH_FAILURE);
  SMESH_TEST_EQ(bad.fill(TRI10), SMESH_FAILURE);
  SMESH_TEST_ASSERT(!LocalSideTable::supported(TET20));
  SMESH_TEST_ASSERT(LocalSideTable::supported(HEX8));
  SMESH_TEST_ASSERT(LocalSideTable::supported(EDGE3));
  SMESH_TEST_EQ(shell_type(NODE1), NODE1);
  return SMESH_TEST_SUCCESS;
}

int test_edge2_sideset_create_skin_nodeset() {
  auto elems = create_host_buffer<idx_t>(2, 2);
  auto pts = create_host_buffer<geom_t>(2, 3);
  SMESH_TEST_ASSERT(elems != nullptr);
  SMESH_TEST_ASSERT(pts != nullptr);
  elems->data()[0][0] = 0;
  elems->data()[1][0] = 1;
  elems->data()[0][1] = 1;
  elems->data()[1][1] = 2;
  pts->data()[0][0] = 0;
  pts->data()[1][0] = 0;
  pts->data()[0][1] = 1;
  pts->data()[1][1] = 0;
  pts->data()[0][2] = 2;
  pts->data()[1][2] = 0;
  auto mesh = std::make_shared<Mesh>(Communicator::self(), EDGE2, elems, pts);
  SMESH_TEST_ASSERT(mesh != nullptr);

  auto left = Sideset::create_from_selector(
      mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(left.size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(left[0]->size(), static_cast<ptrdiff_t>(1));
  SMESH_TEST_EQ(left[0]->parent()->data()[0], static_cast<element_idx_t>(0));
  SMESH_TEST_EQ(left[0]->lfi()->data()[0], static_cast<i16>(0));

  auto skin = skin_sideset(mesh);
  SMESH_TEST_ASSERT(skin != nullptr);
  SMESH_TEST_EQ(skin->size(), static_cast<ptrdiff_t>(2));

  auto nodes = create_nodeset_from_sideset(mesh, left[0]);
  SMESH_TEST_ASSERT(nodes != nullptr);
  SMESH_TEST_EQ(nodes->size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(nodes->data()[0], static_cast<idx_t>(0));

  auto [st, surf] = create_surface_from_sideset(mesh, left[0]);
  SMESH_TEST_EQ(st, NODE1);
  SMESH_TEST_ASSERT(surf != nullptr);
  SMESH_TEST_EQ(surf->extent(0), static_cast<size_t>(1));
  SMESH_TEST_EQ(surf->extent(1), static_cast<size_t>(1));
  SMESH_TEST_EQ(surf->data()[0][0], static_cast<idx_t>(0));
  return SMESH_TEST_SUCCESS;
}

int test_quadshell4_and_tet10_tri6_sidesets() {
  auto q4 = Mesh::create_quad4_square(Communicator::self(), 2, 2);
  SMESH_TEST_ASSERT(q4 != nullptr);
  auto q4_ss = Sideset::create_from_selector(
      q4, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(q4_ss.size(), static_cast<size_t>(1));

  q4->set_element_type(0, QUADSHELL4);
  auto qshell_ss = Sideset::create_from_selector(
      q4, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(qshell_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(qshell_ss[0]->size(), q4_ss[0]->size());

  auto tet10 = Mesh::create_cube(Communicator::self(), TET10, 1, 1, 1);
  SMESH_TEST_ASSERT(tet10 != nullptr);
  auto tet_ss = Sideset::create_from_selector(
      tet10, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(tet_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_ASSERT(tet_ss[0]->size() > 0);
  auto tet_ns = create_nodeset_from_sideset(tet10, tet_ss[0]);
  SMESH_TEST_ASSERT(tet_ns != nullptr);
  SMESH_TEST_ASSERT(tet_ns->size() >= 6);

  auto tri6 = Mesh::create_square(Communicator::self(), TRI6, 2, 2);
  SMESH_TEST_ASSERT(tri6 != nullptr);
  auto tri_ss = Sideset::create_from_selector(
      tri6, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(tri_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_ASSERT(tri_ss[0]->size() > 0);

  auto q9_elems = create_host_buffer<idx_t>(9, 1);
  auto q9_pts = create_host_buffer<geom_t>(2, 9);
  SMESH_TEST_ASSERT(q9_elems != nullptr);
  SMESH_TEST_ASSERT(q9_pts != nullptr);
  for (int i = 0; i < 9; ++i) {
    q9_elems->data()[i][0] = static_cast<idx_t>(i);
  }
  const geom_t q9_xy[9][2] = {
      {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0}, {1, 0.5}, {0.5, 1}, {0, 0.5}, {0.5, 0.5},
  };
  for (int i = 0; i < 9; ++i) {
    q9_pts->data()[0][i] = q9_xy[i][0];
    q9_pts->data()[1][i] = q9_xy[i][1];
  }
  auto q9 = std::make_shared<Mesh>(Communicator::self(), QUAD9, q9_elems, q9_pts);
  auto q9_ss = Sideset::create_from_selector(
      q9, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  SMESH_TEST_EQ(q9_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(q9_ss[0]->size(), static_cast<ptrdiff_t>(1));
  SMESH_TEST_EQ(q9_ss[0]->lfi()->data()[0], static_cast<i16>(3));
  auto q9_ns = create_nodeset_from_sideset(q9, q9_ss[0]);
  SMESH_TEST_EQ(q9_ns->size(), static_cast<size_t>(3));
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

static std::shared_ptr<Mesh> split_first_half(const std::shared_ptr<Mesh> &mesh) {
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

static std::set<FaceKey> unique_face_keys(const std::shared_ptr<Mesh> &mesh,
                                          const std::vector<std::shared_ptr<Sideset>> &sidesets) {
  std::set<FaceKey> keys;
  for (const auto &ss : sidesets) {
    if (!ss || ss->size() == 0) {
      continue;
    }
    for (auto &k : corner_face_keys(mesh, ss)) {
      keys.insert(std::move(k));
    }
  }
  return keys;
}

int test_multiblock_sideset_interface_dedup() {
  auto single = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
  auto multi = split_first_half(single);
  SMESH_TEST_ASSERT(single != nullptr);
  SMESH_TEST_ASSERT(multi != nullptr);
  SMESH_TEST_EQ(static_cast<int>(multi->n_blocks()), 2);

  auto single_ss = Sideset::create_from_plane(single, 1, 0, 0, 0.5);
  auto multi_ss = Sideset::create_from_plane(multi, 1, 0, 0, 0.5);
  SMESH_TEST_ASSERT(!single_ss.empty());
  SMESH_TEST_ASSERT(!multi_ss.empty());

  ptrdiff_t n_single = 0;
  ptrdiff_t n_multi = 0;
  for (const auto &ss : single_ss) {
    n_single += ss->size();
  }
  for (const auto &ss : multi_ss) {
    n_multi += ss->size();
  }
  const auto keys_single = unique_face_keys(single, single_ss);
  const auto keys_multi = unique_face_keys(multi, multi_ss);
  SMESH_TEST_EQ(n_single, static_cast<ptrdiff_t>(keys_single.size()));
  SMESH_TEST_EQ(n_multi, static_cast<ptrdiff_t>(keys_multi.size()));
  SMESH_TEST_EQ(keys_single.size(), keys_multi.size());
  SMESH_TEST_ASSERT(keys_single == keys_multi);
  SMESH_TEST_EQ(keys_single.size(), static_cast<size_t>(4));
  return SMESH_TEST_SUCCESS;
}

int test_create_surface_from_sidesets_merge() {
  auto single = Mesh::create_hex8_cube(Communicator::self(), 2, 2, 2);
  auto multi = Mesh::create_hex8_checkerboard_cube(Communicator::self(), 2, 2, 2);
  SMESH_TEST_ASSERT(single != nullptr);
  SMESH_TEST_ASSERT(multi != nullptr);

  auto single_ss = Sideset::create_from_plane(single, 1, 0, 0, 0.0);
  auto multi_ss = Sideset::create_from_plane(multi, 1, 0, 0, 0.0);
  SMESH_TEST_ASSERT(!single_ss.empty());
  SMESH_TEST_ASSERT(multi_ss.size() > 1);

  auto [st_single, surf_single] = create_surface_from_sidesets(single, single_ss);
  auto [st_multi, surf_multi] = create_surface_from_sidesets(multi, multi_ss);
  SMESH_TEST_EQ(st_single, QUADSHELL4);
  SMESH_TEST_EQ(st_multi, QUADSHELL4);
  SMESH_TEST_ASSERT(surf_single != nullptr);
  SMESH_TEST_ASSERT(surf_multi != nullptr);
  SMESH_TEST_EQ(surf_single->extent(1), surf_multi->extent(1));
  SMESH_TEST_EQ(static_cast<ptrdiff_t>(surf_multi->extent(1)),
                static_cast<ptrdiff_t>(unique_face_keys(multi, multi_ss).size()));

  auto skins = skin_sidesets(multi);
  auto [st_skin, surf_skin] = create_surface_from_sidesets(multi, skins);
  auto serial_skin = skin_sideset(single);
  SMESH_TEST_ASSERT(serial_skin != nullptr);
  auto [st_sskin, surf_sskin] = create_surface_from_sideset(single, serial_skin);
  SMESH_TEST_EQ(st_skin, QUADSHELL4);
  SMESH_TEST_EQ(st_sskin, QUADSHELL4);
  SMESH_TEST_ASSERT(surf_skin != nullptr);
  SMESH_TEST_EQ(surf_skin->extent(1), surf_sskin->extent(1));
  return SMESH_TEST_SUCCESS;
}

int test_mesh_sideset_folder_io() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  mesh->add_sideset("left", sideset);
  SMESH_TEST_EQ(mesh->sidesets().size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(mesh->sidesets("left").size(), static_cast<size_t>(1));

  const auto token = static_cast<long long>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mesh_sideset_io_%lld", token);
  const Path path(path_buffer);
  std::filesystem::remove_all(path.to_string());

  SMESH_TEST_ASSERT(mesh->write(path) == SMESH_SUCCESS);
  SMESH_TEST_ASSERT((path / "sidesets" / "left" / "meta.yaml").exists());

  auto loaded = Mesh::create_from_file(Communicator::self(), path);
  SMESH_TEST_ASSERT(loaded != nullptr);
  auto loaded_ss = loaded->sidesets("left");
  SMESH_TEST_EQ(loaded_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(loaded_ss[0]->block_id(), sideset->block_id());
  SMESH_TEST_EQ(loaded_ss[0]->size(), sideset->size());
  for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
    SMESH_TEST_EQ(loaded_ss[0]->parent()->data()[i], sideset->parent()->data()[i]);
    SMESH_TEST_EQ(loaded_ss[0]->lfi()->data()[i], sideset->lfi()->data()[i]);
  }

  auto cloned = mesh->clone();
  SMESH_TEST_ASSERT(cloned != nullptr);
  auto cloned_ss = cloned->sidesets("left");
  SMESH_TEST_EQ(cloned_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_ASSERT(cloned_ss[0].get() != sideset.get());
  SMESH_TEST_EQ(cloned_ss[0]->size(), sideset->size());
  SMESH_TEST_EQ(cloned_ss[0]->parent()->data()[0], sideset->parent()->data()[0]);

  std::filesystem::remove_all(path.to_string());
  return SMESH_TEST_SUCCESS;
}

int test_mesh_multiblock_sideset_folder_io() {
  auto mesh = Mesh::create_hex8_tet4_cube(Communicator::self(), 2, 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);
  auto left = Sideset::create_from_plane(mesh, 1, 0, 0, 0.0);
  SMESH_TEST_ASSERT(left.size() > 1);
  mesh->add_sidesets("left", left);

  const auto token = static_cast<long long>(
      std::chrono::steady_clock::now().time_since_epoch().count() + 1);
  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mesh_mb_sideset_io_%lld", token);
  const Path path(path_buffer);
  std::filesystem::remove_all(path.to_string());

  SMESH_TEST_ASSERT(mesh->write(path) == SMESH_SUCCESS);
  SMESH_TEST_ASSERT((path / "sidesets" / "left").is_dir());
  SMESH_TEST_ASSERT(!(path / "sidesets" / "left" / "meta.yaml").exists());

  auto loaded = Mesh::create_from_file(Communicator::self(), path);
  SMESH_TEST_ASSERT(loaded != nullptr);
  auto loaded_ss = loaded->sidesets("left");
  SMESH_TEST_EQ(loaded_ss.size(), left.size());
  for (size_t i = 0; i < left.size(); ++i) {
    SMESH_TEST_EQ(loaded_ss[i]->block_id(), left[i]->block_id());
    SMESH_TEST_EQ(loaded_ss[i]->size(), left[i]->size());
    for (ptrdiff_t s = 0; s < left[i]->size(); ++s) {
      SMESH_TEST_EQ(loaded_ss[i]->parent()->data()[s], left[i]->parent()->data()[s]);
      SMESH_TEST_EQ(loaded_ss[i]->lfi()->data()[s], left[i]->lfi()->data()[s]);
    }
  }

  std::filesystem::remove_all(path.to_string());
  return SMESH_TEST_SUCCESS;
}

int test_registered_sideset_remap_on_sfc_reorder() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  mesh->add_sideset("left", sideset);
  const auto keys_before = corner_face_keys(mesh, sideset);
  SMESH_TEST_ASSERT(!keys_before.empty());

  SMESH_TEST_ASSERT(SFC("random3").reorder(*mesh) == SMESH_SUCCESS);

  const auto keys_mapped = corner_face_keys(mesh, sideset);
  auto recreated = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(recreated != nullptr);
  const auto keys_recreated = corner_face_keys(mesh, recreated);
  SMESH_TEST_EQ(keys_mapped.size(), keys_before.size());
  SMESH_TEST_EQ(keys_recreated.size(), keys_before.size());
  for (size_t i = 0; i < keys_before.size(); ++i) {
    SMESH_TEST_EQ(keys_mapped[i].size(), keys_before[i].size());
    for (size_t k = 0; k < keys_before[i].size(); ++k) {
      SMESH_TEST_EQ(keys_mapped[i][k][0], keys_recreated[i][k][0]);
      SMESH_TEST_EQ(keys_mapped[i][k][1], keys_recreated[i][k][1]);
      SMESH_TEST_EQ(keys_mapped[i][k][2], keys_recreated[i][k][2]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

int test_registered_sideset_remap_from_tags() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);
  mesh->add_sideset("left", sideset);
  const auto keys_before = corner_face_keys(mesh, sideset);

  auto tags = create_host_buffer<idx_t>(mesh->n_elements(0));
  for (ptrdiff_t i = 0; i < mesh->n_elements(0); ++i) {
    tags->data()[i] = static_cast<idx_t>(i & 1);
  }
  mesh->reorder_elements_from_tags(0, tags);

  const auto keys_mapped = corner_face_keys(mesh, sideset);
  auto recreated = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(recreated != nullptr);
  const auto keys_recreated = corner_face_keys(mesh, recreated);
  SMESH_TEST_EQ(keys_mapped.size(), keys_before.size());
  for (size_t i = 0; i < keys_before.size(); ++i) {
    for (size_t k = 0; k < keys_before[i].size(); ++k) {
      SMESH_TEST_EQ(keys_mapped[i][k][0], keys_recreated[i][k][0]);
      SMESH_TEST_EQ(keys_mapped[i][k][1], keys_recreated[i][k][1]);
      SMESH_TEST_EQ(keys_mapped[i][k][2], keys_recreated[i][k][2]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

int test_split_mixed_arity_wedge_and_pyramid() {
  auto tri = Mesh::create_tri3_square(Communicator::self(), 2, 2);
  auto wedge = extrude(tri, 1.0, 1);
  SMESH_TEST_ASSERT(wedge != nullptr);
  SMESH_TEST_EQ(wedge->element_type(0), WEDGE6);

  auto skin = skin_sideset(wedge);
  SMESH_TEST_ASSERT(skin != nullptr);
  SMESH_TEST_ASSERT(skin->size() > 0);

  auto parts = split_mixed_arity_sideset(wedge, skin);
  SMESH_TEST_EQ(parts.size(), static_cast<size_t>(2));
  SMESH_TEST_EQ(parts[0]->size() + parts[1]->size(), skin->size());
  SMESH_TEST_EQ(parts[0]->block_id(), skin->block_id());
  SMESH_TEST_EQ(parts[1]->block_id(), skin->block_id());

  const i16 *lfi_tri = parts[0]->lfi()->data();
  for (ptrdiff_t i = 0; i < parts[0]->size(); ++i) {
    SMESH_TEST_EQ(side_type(WEDGE6, lfi_tri[i]), TRI3);
  }
  const i16 *lfi_quad = parts[1]->lfi()->data();
  for (ptrdiff_t i = 0; i < parts[1]->size(); ++i) {
    SMESH_TEST_EQ(side_type(WEDGE6, lfi_quad[i]), QUAD4);
  }

  auto surfaces = create_surfaces_from_sidesets(wedge, {skin});
  SMESH_TEST_EQ(surfaces.size(), static_cast<size_t>(2));
  ptrdiff_t n_extracted = 0;
  bool has_tri_shell = false;
  bool has_quad_shell = false;
  for (const auto &kv : surfaces) {
    SMESH_TEST_ASSERT(kv.second != nullptr);
    n_extracted += static_cast<ptrdiff_t>(kv.second->extent(1));
    has_tri_shell = has_tri_shell || (kv.first == TRISHELL3);
    has_quad_shell = has_quad_shell || (kv.first == QUADSHELL4);
  }
  SMESH_TEST_EQ(n_extracted, skin->size());
  SMESH_TEST_ASSERT(has_tri_shell);
  SMESH_TEST_ASSERT(has_quad_shell);

  auto hex = make_test_mesh();
  auto left = make_left_boundary_sideset(hex);
  SMESH_TEST_ASSERT(left != nullptr);
  auto hex_parts = split_mixed_arity_sideset(hex, left);
  SMESH_TEST_EQ(hex_parts.size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(hex_parts[0]->size(), left->size());

  auto hexdom = Mesh::create_hex_dominant_serial(Communicator::self());
  SMESH_TEST_ASSERT(hexdom != nullptr);
  auto pyr_ss = Sideset::create_from_selector(
      hexdom, [](const geom_t, const geom_t, const geom_t) { return true; }, {"pyramid"});
  SMESH_TEST_EQ(pyr_ss.size(), static_cast<size_t>(1));
  SMESH_TEST_EQ(pyr_ss[0]->size(), static_cast<ptrdiff_t>(5));
  auto pyr_parts = split_mixed_arity_sideset(hexdom, pyr_ss[0]);
  SMESH_TEST_EQ(pyr_parts.size(), static_cast<size_t>(2));
  SMESH_TEST_EQ(pyr_parts[0]->size() + pyr_parts[1]->size(), pyr_ss[0]->size());
  SMESH_TEST_EQ(pyr_parts[0]->size(), static_cast<ptrdiff_t>(4));
  SMESH_TEST_EQ(pyr_parts[1]->size(), static_cast<ptrdiff_t>(1));

  auto wedge_ss = Sideset::create_from_selector(
      hexdom, [](const geom_t, const geom_t, const geom_t) { return true; }, {"wedge"});
  SMESH_TEST_EQ(wedge_ss.size(), static_cast<size_t>(1));
  auto wedge_parts = split_mixed_arity_sideset(hexdom, wedge_ss[0]);
  SMESH_TEST_EQ(wedge_parts.size(), static_cast<size_t>(2));
  SMESH_TEST_EQ(wedge_parts[0]->size() + wedge_parts[1]->size(), wedge_ss[0]->size());
  return SMESH_TEST_SUCCESS;
}

static int check_aos_matches_soa(const Path &aos_path, const std::shared_ptr<Mesh> &mesh,
                                 const block_idx_t bid) {
  const int nxe = mesh->n_nodes_per_element(bid);
  const ptrdiff_t ne = mesh->n_elements(bid);
  FILE *fp = fopen(aos_path.c_str(), "rb");
  SMESH_TEST_ASSERT(fp != nullptr);
  auto row = create_host_buffer<idx_t>((size_t)nxe);
  idx_t *const d_row = row->data();
  idx_t *const *elems = mesh->elements(bid)->data();
  for (ptrdiff_t e = 0; e < ne; ++e) {
    SMESH_TEST_ASSERT(fread(d_row, sizeof(idx_t), (size_t)nxe, fp) == (size_t)nxe);
    for (int d = 0; d < nxe; ++d) {
      SMESH_TEST_EQ(d_row[d], elems[d][e]);
    }
  }
  unsigned char extra = 0;
  SMESH_TEST_ASSERT(fread(&extra, 1, 1, fp) == 0);
  fclose(fp);
  return SMESH_TEST_SUCCESS;
}

int test_write_with_xdmf() {
  auto mesh = make_test_mesh();
  auto sideset = make_left_boundary_sideset(mesh);
  SMESH_TEST_ASSERT(mesh != nullptr);
  SMESH_TEST_ASSERT(sideset != nullptr);

  const auto token = static_cast<long long>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer), "/tmp/smesh_write_xdmf_%lld", token);
  const Path path(path_buffer);
  std::filesystem::remove_all(path.to_string());

  mesh->add_sideset("left", sideset);
  SMESH_TEST_ASSERT(mesh->write_with_xdmf(path) == SMESH_SUCCESS);
  SMESH_TEST_ASSERT((path / "mesh.xdmf").exists());
  const Path aos = path / (std::string("connectivity.") + std::string(TypeToString<idx_t>::value()));
  SMESH_TEST_ASSERT(aos.exists());
  SMESH_TEST_EQ(check_aos_matches_soa(aos, mesh, 0), SMESH_TEST_SUCCESS);
  const Path surf = path / "sidesets" / "left" /
                    (std::string("surface.") + std::string(TypeToString<idx_t>::value()));
  SMESH_TEST_ASSERT(surf.exists());

  std::ifstream xdmf((path / "mesh.xdmf").to_string());
  SMESH_TEST_ASSERT(xdmf.good());
  std::string xml((std::istreambuf_iterator<char>(xdmf)), std::istreambuf_iterator<char>());
  SMESH_TEST_ASSERT(xml.find("Hexahedron") != std::string::npos);
  SMESH_TEST_ASSERT(xml.find("Quadrilateral") != std::string::npos);
  SMESH_TEST_ASSERT(xml.find("sidesets/left/surface.") != std::string::npos);

  auto loaded = Mesh::create_from_file(Communicator::self(), path);
  SMESH_TEST_ASSERT(loaded != nullptr);
  SMESH_TEST_EQ(loaded->n_elements(), mesh->n_elements());
  SMESH_TEST_EQ(loaded->n_nodes(), mesh->n_nodes());

  std::filesystem::remove_all(path.to_string());
  return SMESH_TEST_SUCCESS;
}

int test_write_with_xdmf_multiblock() {
  auto mesh = Mesh::create_hex8_tet4_cube(Communicator::self(), 2, 2, 2);
  SMESH_TEST_ASSERT(mesh != nullptr);

  const auto token = static_cast<long long>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer), "/tmp/smesh_write_xdmf_mb_%lld", token);
  const Path path(path_buffer);
  std::filesystem::remove_all(path.to_string());

  SMESH_TEST_ASSERT(mesh->write_with_xdmf(path) == SMESH_SUCCESS);
  SMESH_TEST_ASSERT((path / "mesh.xdmf").exists());
  for (size_t b = 0; b < mesh->n_blocks(); ++b) {
    const Path aos = path / "blocks" / mesh->block(b)->name() /
                     (std::string("connectivity.") + std::string(TypeToString<idx_t>::value()));
    SMESH_TEST_ASSERT(aos.exists());
    SMESH_TEST_EQ(check_aos_matches_soa(aos, mesh, static_cast<block_idx_t>(b)),
                  SMESH_TEST_SUCCESS);
  }

  std::ifstream xdmf((path / "mesh.xdmf").to_string());
  std::string xml((std::istreambuf_iterator<char>(xdmf)), std::istreambuf_iterator<char>());
  SMESH_TEST_ASSERT(xml.find("Collection") != std::string::npos);
  SMESH_TEST_ASSERT(xml.find("blocks/") != std::string::npos);

  std::filesystem::remove_all(path.to_string());
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
  SMESH_RUN_TEST(test_sideset_remap_on_sfc_reorder);
  SMESH_RUN_TEST(test_sideset_remap_from_tags);
  SMESH_RUN_TEST(test_sideset_map_hex_refine);
  SMESH_RUN_TEST(test_sideset_map_tet_refine);
  SMESH_RUN_TEST(test_sideset_map_tri_refine);
  SMESH_RUN_TEST(test_sideset_map_refine_unsupported);
  SMESH_RUN_TEST(test_sideset_select_propagate_cube_mesh);
  SMESH_RUN_TEST(test_hex27_element_contract);
  SMESH_RUN_TEST(test_hex27_cube_uses_conventional_ordering);
  SMESH_RUN_TEST(test_local_side_table_higher_order_exodus);
  SMESH_RUN_TEST(test_local_side_table_edge_and_shell_aliases);
  SMESH_RUN_TEST(test_edge2_sideset_create_skin_nodeset);
  SMESH_RUN_TEST(test_quadshell4_and_tet10_tri6_sidesets);
  SMESH_RUN_TEST(test_upper_triangular_graph_with_unused_node);
  SMESH_RUN_TEST(test_multiblock_sideset_interface_dedup);
  SMESH_RUN_TEST(test_create_surface_from_sidesets_merge);
  SMESH_RUN_TEST(test_mesh_sideset_folder_io);
  SMESH_RUN_TEST(test_mesh_multiblock_sideset_folder_io);
  SMESH_RUN_TEST(test_registered_sideset_remap_on_sfc_reorder);
  SMESH_RUN_TEST(test_registered_sideset_remap_from_tags);
  SMESH_RUN_TEST(test_split_mixed_arity_wedge_and_pyramid);
  SMESH_RUN_TEST(test_write_with_xdmf);
  SMESH_RUN_TEST(test_write_with_xdmf_multiblock);

  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
