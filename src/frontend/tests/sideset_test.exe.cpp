#include <chrono>
#include <cstdio>
#include <filesystem>

#include "smesh_mesh.hpp"
#include "smesh_sideset.hpp"
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

// TODO: create a test for sideset_select_propagate on the cube mesh

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);

  SMESH_RUN_TEST(test_sideset_creation);
  SMESH_RUN_TEST(test_sideset_to_nodeset_conversion);
  SMESH_RUN_TEST(test_sideset_io_write_read_identity);

  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
