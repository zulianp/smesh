#include <algorithm>
#include <array>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <map>
#include <vector>

#include "smesh_alloc.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_distributed_reorder.hpp"
#include "smesh_mesh.hpp"
#include "smesh_mesh_reorder.hpp"
#include "smesh_test.hpp"

using namespace smesh;

int test_parallel_sfc_sort_and_gather() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr));
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_parallel_sfc_test_%d_%d", comm->size(), token);
  const Path root(path_buffer);
  const Path mesh_path = root / "mesh";

  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
    std::filesystem::create_directories(root.to_string());

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 8);
    auto mesh = Mesh::create_tet4_cube(Communicator::self(), nx, 5, 4);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->write(mesh_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  int nnodesxelem = 0;
  ptrdiff_t n_local_elements = 0;
  ptrdiff_t n_global_elements = 0;
  idx_t **elements = nullptr;
  SMESH_TEST_ASSERT(
      mesh_block_from_folder<idx_t>(comm->get(), mesh_path, &nnodesxelem,
                                    &elements, &n_local_elements,
                                    &n_global_elements) == SMESH_SUCCESS);

  int spatial_dim = 0;
  ptrdiff_t n_local_nodes = 0;
  ptrdiff_t n_global_nodes = 0;
  geom_t **points = nullptr;
  SMESH_TEST_ASSERT(mesh_coordinates_from_folder<geom_t>(
                        comm->get(), mesh_path, &spatial_dim, &points,
                        &n_local_nodes, &n_global_nodes) == SMESH_SUCCESS);
  SMESH_TEST_EQ(nnodesxelem, 4);
  SMESH_TEST_EQ(spatial_dim, 3);

  std::vector<large_idx_t> sorted_ids((size_t)n_local_elements);
  const int reorder_status = distributed_reorder_elements<idx_t, geom_t>(
      comm->get(), nnodesxelem, n_local_elements, n_global_elements, elements,
      n_global_nodes, points, sorted_ids.data());
  SMESH_TEST_ASSERT(reorder_status == SMESH_SUCCESS);

  for (int d = 0; d < nnodesxelem; ++d) {
    SMESH_FREE(elements[d]);
  }
  SMESH_FREE(elements);

  for (int d = 0; d < spatial_dim; ++d) {
    SMESH_FREE(points[d]);
  }
  SMESH_FREE(points);

  comm->barrier();
  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
  }
  comm->barrier();

  SMESH_UNUSED(n_local_nodes);
  SMESH_UNUSED(n_global_nodes);
  return SMESH_TEST_SUCCESS;
#endif
}

int test_parallel_sfc_matches_serial() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr));
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_parallel_sfc_serial_test_%d_%d", comm->size(),
                token);
  const Path root(path_buffer);
  const Path mesh_path = root / "mesh";

  std::vector<large_idx_t> expected_ids;
  ptrdiff_t n_serial_elements = 0;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
    std::filesystem::create_directories(root.to_string());

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 8);
    auto original_mesh = Mesh::create_tet4_cube(Communicator::self(), nx, 5, 4);
    auto serial_mesh = Mesh::create_tet4_cube(Communicator::self(), nx, 5, 4);
    SMESH_TEST_ASSERT(original_mesh != nullptr);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SMESH_TEST_ASSERT(original_mesh->write(mesh_path) == SMESH_SUCCESS);

    const int nxe = original_mesh->n_nodes_per_element(0);
    n_serial_elements = original_mesh->n_elements(0);
    auto original_elements = original_mesh->elements(0)->data();
    auto original_points = original_mesh->points()->data();

    std::map<std::array<geom_t, 3>, large_idx_t> original_id;
    for (ptrdiff_t e = 0; e < n_serial_elements; ++e) {
      std::array<geom_t, 3> key = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
        for (int i = 0; i < nxe; ++i) {
          key[(size_t)d] += original_points[d][original_elements[i][e]];
        }
        key[(size_t)d] /= static_cast<geom_t>(nxe);
      }
      SMESH_TEST_ASSERT(original_id.emplace(key, (large_idx_t)e).second);
    }

    SMESH_TEST_ASSERT(SFC("hilbert3").reorder(*serial_mesh) == SMESH_SUCCESS);

    expected_ids.resize((size_t)n_serial_elements);
    auto serial_elements = serial_mesh->elements(0)->data();
    auto serial_points = serial_mesh->points()->data();
    for (ptrdiff_t e = 0; e < n_serial_elements; ++e) {
      std::array<geom_t, 3> key = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
        for (int i = 0; i < nxe; ++i) {
          key[(size_t)d] += serial_points[d][serial_elements[i][e]];
        }
        key[(size_t)d] /= static_cast<geom_t>(nxe);
      }
      auto iter = original_id.find(key);
      SMESH_TEST_ASSERT(iter != original_id.end());
      expected_ids[(size_t)e] = iter->second;
    }
  }

  comm->barrier();

  int nnodesxelem = 0;
  ptrdiff_t n_local_elements = 0;
  ptrdiff_t n_global_elements = 0;
  idx_t **elements = nullptr;
  SMESH_TEST_ASSERT(
      mesh_block_from_folder<idx_t>(comm->get(), mesh_path, &nnodesxelem,
                                    &elements, &n_local_elements,
                                    &n_global_elements) == SMESH_SUCCESS);

  int spatial_dim = 0;
  ptrdiff_t n_local_nodes = 0;
  ptrdiff_t n_global_nodes = 0;
  geom_t **points = nullptr;
  SMESH_TEST_ASSERT(mesh_coordinates_from_folder<geom_t>(
                        comm->get(), mesh_path, &spatial_dim, &points,
                        &n_local_nodes, &n_global_nodes) == SMESH_SUCCESS);
  SMESH_TEST_EQ(nnodesxelem, 4);
  SMESH_TEST_EQ(spatial_dim, 3);

  std::vector<large_idx_t> parallel_ids((size_t)n_local_elements);
  SMESH_TEST_ASSERT((distributed_reorder_elements<idx_t, geom_t>(
                        comm->get(), nnodesxelem, n_local_elements,
                        n_global_elements, elements, n_global_nodes, points,
                        parallel_ids.data()) == SMESH_SUCCESS));

  std::vector<int> counts((size_t)comm->size());
  std::vector<int> displs((size_t)comm->size());
  for (int r = 0; r < comm->size(); ++r) {
    counts[(size_t)r] = static_cast<int>(rank_split(n_global_elements,
                                                    comm->size(), r));
    displs[(size_t)r] = static_cast<int>(rank_start(n_global_elements,
                                                    comm->size(), r));
  }

  std::vector<large_idx_t> expected_local((size_t)n_local_elements);
  SMESH_MPI_CATCH(MPI_Scatterv(
      comm->rank() == 0 ? expected_ids.data() : nullptr, counts.data(),
      displs.data(), mpi_type<large_idx_t>(), expected_local.data(),
      static_cast<int>(n_local_elements), mpi_type<large_idx_t>(), 0,
      comm->get()));

  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    SMESH_TEST_EQ(parallel_ids[(size_t)e], expected_local[(size_t)e]);
  }

  for (int d = 0; d < nnodesxelem; ++d) {
    SMESH_FREE(elements[d]);
  }
  SMESH_FREE(elements);

  for (int d = 0; d < spatial_dim; ++d) {
    SMESH_FREE(points[d]);
  }
  SMESH_FREE(points);

  comm->barrier();
  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
  }
  comm->barrier();

  SMESH_UNUSED(n_local_nodes);
  SMESH_UNUSED(n_serial_elements);
  return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_parallel_sfc_sort_and_gather);
  SMESH_RUN_TEST(test_parallel_sfc_matches_serial);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
