#include <algorithm>
#include <cstdio>
#include <filesystem>

#ifdef SMESH_ENABLE_MPI
#include <unistd.h>
#endif

#include "smesh_mesh.hpp"
#include "smesh_test.hpp"
#include "smesh_mesh_reorder.hpp"

using namespace smesh;

int test_owned_elements_point_to_owned_nodes() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();
  if (comm->size() == 1) {
    return SMESH_TEST_SUCCESS;
  }

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_owned_elements_owned_nodes_test_%d",
                comm->size());
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(2, comm->size());
  const ptrdiff_t ny = 10;
  const ptrdiff_t nz = 10;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());

    auto serial_mesh = Mesh::create_tet4_cube(Communicator::self(), nx, ny, nz);
    SFC::create_from_env()->reorder(*serial_mesh);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SMESH_TEST_ASSERT(serial_mesh->write(mesh_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  auto mesh = std::make_shared<Mesh>(comm);
  SMESH_TEST_ASSERT(mesh->read(mesh_path) == SMESH_SUCCESS);

  auto dist = mesh->distributed();
  SMESH_TEST_ASSERT(dist != nullptr);

  const ptrdiff_t n_owned_nodes = dist->n_nodes_owned();
  const ptrdiff_t n_owned_not_shared_elements =
      dist->n_elements_owned_not_shared();
  const ptrdiff_t n_owned_elements = dist->n_elements_owned();
  const int nxe = mesh->n_nodes_per_element(0);
  auto elements = mesh->elements(0)->data();

  for (ptrdiff_t e = 0; e < n_owned_not_shared_elements; ++e) {
    for (int d = 0; d < nxe; ++d) {
      SMESH_TEST_ASSERT(elements[d][e] < n_owned_nodes);
    }
  }

  for (ptrdiff_t e = n_owned_not_shared_elements; e < n_owned_elements; ++e) {
    bool touches_non_owned_node = false;
    for (int d = 0; d < nxe; ++d) {
      touches_non_owned_node |= (elements[d][e] >= n_owned_nodes);
    }

    SMESH_TEST_ASSERT(touches_non_owned_node);
  }

  comm->barrier();
  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  comm->barrier();

  return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_owned_elements_point_to_owned_nodes);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
