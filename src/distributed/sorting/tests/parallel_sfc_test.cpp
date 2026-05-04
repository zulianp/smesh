#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <vector>

#include "smesh_alloc.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_distributed_reorder.hpp"
#include "smesh_mesh.hpp"
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

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_parallel_sfc_sort_and_gather);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
