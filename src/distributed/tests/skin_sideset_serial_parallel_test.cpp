#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

static std::shared_ptr<Mesh> pad_to_3d(const std::shared_ptr<Mesh> &m) {
  if (m->spatial_dimension() >= 3) {
    return m;
  }
  const ptrdiff_t n = m->n_nodes();
  auto p3 = create_host_buffer<geom_t>(3, static_cast<size_t>(n));
  geom_t **const d = p3->data();
  geom_t **const s = m->points()->data();
  std::memcpy(d[0], s[0], static_cast<size_t>(n) * sizeof(geom_t));
  if (m->spatial_dimension() > 1) {
    std::memcpy(d[1], s[1], static_cast<size_t>(n) * sizeof(geom_t));
  } else {
    std::memset(d[1], 0, static_cast<size_t>(n) * sizeof(geom_t));
  }
  std::memset(d[2], 0, static_cast<size_t>(n) * sizeof(geom_t));
  return std::make_shared<Mesh>(Communicator::self(), m->blocks(), p3);
}

static std::shared_ptr<Mesh> make_tet4(ptrdiff_t nx, ptrdiff_t ny, ptrdiff_t nz) {
  return Mesh::create_tet4_cube(Communicator::self(), nx, ny, nz);
}

static std::shared_ptr<Mesh> make_hex8(ptrdiff_t nx, ptrdiff_t ny, ptrdiff_t nz) {
  return Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz);
}

static std::shared_ptr<Mesh> make_quad4_3d(ptrdiff_t nx, ptrdiff_t ny, ptrdiff_t) {
  return pad_to_3d(Mesh::create_quad4_square(Communicator::self(), nx, ny));
}

static std::shared_ptr<Mesh> make_tri3_3d(ptrdiff_t nx, ptrdiff_t ny, ptrdiff_t) {
  return pad_to_3d(Mesh::create_tri3_square(Communicator::self(), nx, ny));
}

static std::shared_ptr<Sideset> left_boundary_sideset(const std::shared_ptr<Mesh> &mesh) {
  auto sidesets = Sideset::create_from_selector(
      mesh, [](const geom_t x, const geom_t, const geom_t) { return x < 1e-12; });
  if (sidesets.size() != 1) {
    return nullptr;
  }
  return sidesets[0];
}

using MakeMeshFn = std::shared_ptr<Mesh> (*)(ptrdiff_t, ptrdiff_t, ptrdiff_t);
using MakeSsFn = std::shared_ptr<Sideset> (*)(const std::shared_ptr<Mesh> &);

static int run_serial_parallel_sideset_parity(const char *tag, MakeMeshFn make_mesh, MakeSsFn make_ss) {
#ifndef SMESH_ENABLE_MPI
  SMESH_UNUSED(tag);
  SMESH_UNUSED(make_mesh);
  SMESH_UNUSED(make_ss);
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
                "/tmp/smesh_skin_sideset_serial_parallel_test_%s_%d_%d",
                tag, comm->size(), token);
  const Path root(path_buffer);
  const Path mesh_path = root / "mesh";
  const Path serial_ss_path = root / "serial_ss";
  const Path parallel_ss_path = root / "parallel_ss";

  const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
  const ptrdiff_t ny = 4;
  const ptrdiff_t nz = 3;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
    std::filesystem::create_directories(root.to_string());

    auto serial_mesh = make_mesh(nx, ny, nz);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SFC("random3").reorder(*serial_mesh);
    SMESH_TEST_ASSERT(serial_mesh->write(mesh_path) == SMESH_SUCCESS);

    auto serial_ss = make_ss(serial_mesh);
    SMESH_TEST_ASSERT(serial_ss != nullptr);
    SMESH_TEST_ASSERT(serial_ss->size() > 0);
    SMESH_TEST_ASSERT(serial_ss->write(serial_ss_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  auto distributed_mesh = std::make_shared<Mesh>(comm);
  SMESH_TEST_ASSERT(distributed_mesh->read(mesh_path) == SMESH_SUCCESS);

  auto parallel_ss = make_ss(distributed_mesh);
  SMESH_TEST_ASSERT(parallel_ss != nullptr);
  SMESH_TEST_ASSERT(parallel_ss->write(parallel_ss_path) == SMESH_SUCCESS);

  const i64 local_n = static_cast<i64>(parallel_ss->size());
  const i64 global_n = comm->sum(local_n);

  comm->barrier();

  if (comm->rank() == 0) {
    auto serial_ss = Sideset::create_from_file(Communicator::self(), serial_ss_path);
    auto parallel_loaded = Sideset::create_from_file(Communicator::self(), parallel_ss_path);

    SMESH_TEST_ASSERT(serial_ss != nullptr);
    SMESH_TEST_ASSERT(parallel_loaded != nullptr);
    SMESH_TEST_EQ(serial_ss->size(), parallel_loaded->size());

    ptrdiff_t meta_size = -1;
    block_idx_t meta_block_id = static_cast<block_idx_t>(-1);
    std::ifstream ifs((parallel_ss_path / "meta.yaml").to_string());
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
    SMESH_TEST_EQ(meta_size, serial_ss->size());
    SMESH_TEST_EQ(meta_size, static_cast<ptrdiff_t>(global_n));
    SMESH_TEST_EQ(meta_block_id, serial_ss->block_id());

    const auto serial_keys = sorted_side_keys(serial_ss);
    const auto parallel_keys = sorted_side_keys(parallel_loaded);

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

int test_skin_sideset_matches_serial() {
  return run_serial_parallel_sideset_parity("tet4_skin", make_tet4, skin_sideset);
}

int test_skin_hex8_serial_parallel() {
  return run_serial_parallel_sideset_parity("hex8_skin", make_hex8, skin_sideset);
}

int test_skin_quad4_serial_parallel() {
  return run_serial_parallel_sideset_parity("quad4_skin", make_quad4_3d, skin_sideset);
}

int test_skin_tri3_serial_parallel() {
  return run_serial_parallel_sideset_parity("tri3_skin", make_tri3_3d, skin_sideset);
}

int test_create_quad4_serial_parallel() {
  return run_serial_parallel_sideset_parity("quad4_create", make_quad4_3d, left_boundary_sideset);
}

int test_create_tri3_serial_parallel() {
  return run_serial_parallel_sideset_parity("tri3_create", make_tri3_3d, left_boundary_sideset);
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_skin_sideset_matches_serial);
  SMESH_RUN_TEST(test_skin_hex8_serial_parallel);
  SMESH_RUN_TEST(test_skin_quad4_serial_parallel);
  SMESH_RUN_TEST(test_skin_tri3_serial_parallel);
  SMESH_RUN_TEST(test_create_quad4_serial_parallel);
  SMESH_RUN_TEST(test_create_tri3_serial_parallel);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
