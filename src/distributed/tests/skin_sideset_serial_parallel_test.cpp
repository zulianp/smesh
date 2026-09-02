#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "smesh_mesh.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"
#include "smesh_mesh_reorder.hpp"

using namespace smesh;

static std::vector<std::pair<element_idx_t, i16>>
sorted_side_keys(const std::shared_ptr<Sideset> &sideset) {
  std::vector<std::pair<element_idx_t, i16>> keys(sideset->size());
  auto parent = sideset->parent()->data();
  auto lfi = sideset->lfi()->data();

  for (ptrdiff_t i = 0; i < sideset->size(); ++i) {
    keys[(size_t)i] = std::make_pair(parent[i], lfi[i]);
  }

  std::sort(keys.begin(), keys.end());
  return keys;
}

int test_skin_sideset_matches_serial() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();
  if (comm->size() == 1) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr));
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_skin_sideset_serial_parallel_test_%d_%d",
                comm->size(), token);
  const Path root(path_buffer);
  const Path mesh_path = root / "mesh";
  const Path serial_skin_path = root / "serial_skin";
  const Path parallel_skin_path = root / "parallel_skin";

  const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
  const ptrdiff_t ny = 4;
  const ptrdiff_t nz = 3;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
    std::filesystem::create_directories(root.to_string());

    auto serial_mesh = Mesh::create_tet4_cube(Communicator::self(), nx, ny, nz);
    SFC("random3").reorder(*serial_mesh);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SMESH_TEST_ASSERT(serial_mesh->write(mesh_path) == SMESH_SUCCESS);

    auto serial_skin = skin_sideset(serial_mesh);
    SMESH_TEST_ASSERT(serial_skin != nullptr);
    SMESH_TEST_ASSERT(serial_skin->write(serial_skin_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  auto distributed_mesh = std::make_shared<Mesh>(comm);
  SMESH_TEST_ASSERT(distributed_mesh->read(mesh_path) == SMESH_SUCCESS);

  auto parallel_skin = skin_sideset(distributed_mesh);
  SMESH_TEST_ASSERT(parallel_skin != nullptr);
  SMESH_TEST_ASSERT(parallel_skin->write(parallel_skin_path) == SMESH_SUCCESS);

  const i64 local_n = static_cast<i64>(parallel_skin->size());
  const i64 global_n = comm->sum(local_n);

  comm->barrier();

  if (comm->rank() == 0) {
    auto serial_skin =
        Sideset::create_from_file(Communicator::self(), serial_skin_path);
    auto parallel_skin =
        Sideset::create_from_file(Communicator::self(), parallel_skin_path);

    SMESH_TEST_ASSERT(serial_skin != nullptr);
    SMESH_TEST_ASSERT(parallel_skin != nullptr);
    SMESH_TEST_EQ(serial_skin->size(), parallel_skin->size());

    ptrdiff_t meta_size = -1;
    block_idx_t meta_block_id = static_cast<block_idx_t>(-1);
    std::ifstream ifs((parallel_skin_path / "meta.yaml").to_string());
    SMESH_TEST_ASSERT(ifs.good());
    std::string line;
    while (std::getline(ifs, line)) {
      const auto start = line.find_first_not_of(" \t");
      if (start == std::string::npos) {
        continue;
      }
      line = line.substr(start);
      if (line.compare(0, 5, "size:") == 0) {
        meta_size = static_cast<ptrdiff_t>(std::strtoll(line.c_str() + 5, nullptr, 10));
      } else if (line.compare(0, 9, "block_id:") == 0) {
        meta_block_id = static_cast<block_idx_t>(std::strtol(line.c_str() + 9, nullptr, 10));
      }
    }
    SMESH_TEST_EQ(meta_size, serial_skin->size());
    SMESH_TEST_EQ(meta_size, static_cast<ptrdiff_t>(global_n));
    SMESH_TEST_EQ(meta_block_id, serial_skin->block_id());

    const auto serial_keys = sorted_side_keys(serial_skin);
    const auto parallel_keys = sorted_side_keys(parallel_skin);

    SMESH_TEST_EQ(serial_keys.size(), parallel_keys.size());
    for (size_t i = 0; i < serial_keys.size(); ++i) {
      SMESH_TEST_EQ(serial_keys[i].first, parallel_keys[i].first);
      SMESH_TEST_EQ(serial_keys[i].second, parallel_keys[i].second);
    }

    std::filesystem::remove_all(root.to_string());
  }

  comm->barrier();
  return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_skin_sideset_matches_serial);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
