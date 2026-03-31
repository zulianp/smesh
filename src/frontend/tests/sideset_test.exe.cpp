#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <utility>
#include <vector>

#include "smesh_mask.hpp"
#include "smesh_mesh.hpp"
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
  for (int i = 0; i < lst.nnxs; ++i) {
    x += points[0][elements[lst(s, i)][e]];
  }

  return x / static_cast<geom_t>(lst.nnxs);
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

  auto read_back = Sideset::create_from_file(Communicator::self(), path);
  SMESH_TEST_ASSERT(read_back != nullptr);
  SMESH_TEST_EQ(read_back->block_id(), sideset->block_id());
  SMESH_TEST_EQ(read_back->size(), sideset->size());

  for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
    SMESH_TEST_EQ(read_back->parent()->data()[i], sideset->parent()->data()[i]);
    SMESH_TEST_EQ(read_back->lfi()->data()[i], sideset->lfi()->data()[i]);
  }

  std::filesystem::remove_all(path.to_string());
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

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);

  SMESH_RUN_TEST(test_sideset_creation);
  SMESH_RUN_TEST(test_sideset_to_nodeset_conversion);
  SMESH_RUN_TEST(test_sideset_io_write_read_identity);
  SMESH_RUN_TEST(test_sideset_select_propagate_cube_mesh);

  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
