#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <limits>
#include <vector>

#include "smesh_alloc.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_aura.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sfc.hpp"
#include "smesh_test.hpp"

extern "C" {
#include "../../../../external/mpi-sort/include/mpi-sort.h"
}

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

  const ptrdiff_t element_start =
      rank_start(n_global_elements, comm->size(), comm->rank());

  std::vector<geom_t> cx((size_t)n_local_elements);
  std::vector<geom_t> cy((size_t)n_local_elements);
  std::vector<geom_t> cz((size_t)n_local_elements);
  std::vector<idx_t> element_nodes((size_t)n_local_elements *
                                   (size_t)nnodesxelem);
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    for (int d = 0; d < nnodesxelem; ++d) {
      element_nodes[(size_t)e * (size_t)nnodesxelem + (size_t)d] =
          elements[d][e];
    }
  }

  std::vector<geom_t> element_x(element_nodes.size());
  std::vector<geom_t> element_y(element_nodes.size());
  std::vector<geom_t> element_z(element_nodes.size());
  SMESH_TEST_ASSERT(gather_mapped_field(comm->get(),
                                        (ptrdiff_t)element_nodes.size(),
                                        n_global_nodes, element_nodes.data(),
                                        mpi_type<geom_t>(), points[0],
                                        element_x.data()) == SMESH_SUCCESS);
  SMESH_TEST_ASSERT(gather_mapped_field(comm->get(),
                                        (ptrdiff_t)element_nodes.size(),
                                        n_global_nodes, element_nodes.data(),
                                        mpi_type<geom_t>(), points[1],
                                        element_y.data()) == SMESH_SUCCESS);
  SMESH_TEST_ASSERT(gather_mapped_field(comm->get(),
                                        (ptrdiff_t)element_nodes.size(),
                                        n_global_nodes, element_nodes.data(),
                                        mpi_type<geom_t>(), points[2],
                                        element_z.data()) == SMESH_SUCCESS);

  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    geom_t x = 0;
    geom_t y = 0;
    geom_t z = 0;
    for (int d = 0; d < nnodesxelem; ++d) {
      const size_t node = (size_t)e * (size_t)nnodesxelem + (size_t)d;
      x += element_x[node];
      y += element_y[node];
      z += element_z[node];
    }

    cx[(size_t)e] = x / nnodesxelem;
    cy[(size_t)e] = y / nnodesxelem;
    cz[(size_t)e] = z / nnodesxelem;
  }

  geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max()};
  geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest()};
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    const geom_t x = cx[(size_t)e];
    const geom_t y = cy[(size_t)e];
    const geom_t z = cz[(size_t)e];
    local_min[0] = std::min(local_min[0], x);
    local_min[1] = std::min(local_min[1], y);
    local_min[2] = std::min(local_min[2], z);
    local_max[0] = std::max(local_max[0], x);
    local_max[1] = std::max(local_max[1], y);
    local_max[2] = std::max(local_max[2], z);
  }

  geom_t global_min[3];
  geom_t global_max[3];
  SMESH_MPI_CATCH(MPI_Allreduce(local_min, global_min, 3, mpi_type<geom_t>(),
                                MPI_MIN, comm->get()));
  SMESH_MPI_CATCH(MPI_Allreduce(local_max, global_max, 3, mpi_type<geom_t>(),
                                MPI_MAX, comm->get()));

  std::vector<u32> send_keys((size_t)n_local_elements);
  SMESH_TEST_ASSERT(
      encode_hilbert3<geom_t>(n_local_elements, cx.data(), cy.data(), cz.data(),
                              global_min[0], global_max[0], global_min[1],
                              global_max[1], global_min[2], global_max[2],
                              send_keys.data()) == SMESH_SUCCESS);

  std::vector<large_idx_t> send_ids((size_t)n_local_elements);
  for (ptrdiff_t e = 0; e < n_local_elements; ++e) {
    send_ids[(size_t)e] = static_cast<large_idx_t>(element_start + e);
  }

  const ptrdiff_t n_sorted_elements =
      rank_split(n_global_elements, comm->size(), comm->rank());
  std::vector<u32> sorted_keys((size_t)n_sorted_elements);
  std::vector<large_idx_t> sorted_ids((size_t)n_sorted_elements);

  SMESH_TEST_ASSERT(MPI_Sort_bykey(send_keys.data(), send_ids.data(),
                                   static_cast<int>(n_local_elements),
                                   mpi_type<u32>(), mpi_type<large_idx_t>(),
                                   sorted_keys.data(), sorted_ids.data(),
                                   static_cast<int>(n_sorted_elements),
                                   comm->get()) == MPI_SUCCESS);

  u32 local_last = n_sorted_elements ? sorted_keys.back() : 0;
  u32 next_first = 0;
  if (comm->rank() + 1 < comm->size()) {
    SMESH_MPI_CATCH(MPI_Sendrecv(
        &local_last, 1, mpi_type<u32>(), comm->rank() + 1, 0, &next_first, 1,
        mpi_type<u32>(), comm->rank() + 1, 1, comm->get(), MPI_STATUS_IGNORE));
    if (n_sorted_elements > 0) {
      SMESH_TEST_ASSERT(local_last <= next_first);
    }
  }
  if (comm->rank() > 0) {
    u32 prev_last = 0;
    SMESH_MPI_CATCH(MPI_Sendrecv(sorted_keys.data(), 1, mpi_type<u32>(),
                                 comm->rank() - 1, 1, &prev_last, 1,
                                 mpi_type<u32>(), comm->rank() - 1, 0,
                                 comm->get(), MPI_STATUS_IGNORE));
    if (n_sorted_elements > 0) {
      SMESH_TEST_ASSERT(prev_last <= sorted_keys.front());
    }
  }

  std::vector<std::vector<idx_t>> sorted_element_storage(
      (size_t)nnodesxelem, std::vector<idx_t>((size_t)n_sorted_elements));
  for (int d = 0; d < nnodesxelem; ++d) {
    SMESH_TEST_ASSERT(
        gather_mapped_field(comm->get(), n_sorted_elements, n_global_elements,
                            sorted_ids.data(), mpi_type<idx_t>(), elements[d],
                            sorted_element_storage[(size_t)d].data()) ==
        SMESH_SUCCESS);
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
