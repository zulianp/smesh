#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>

#include "smesh_distributed_create.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static int test_hex8_cube_create_distributed_matches_file() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr));
  }
  comm->broadcast(&token, 1, 0);

  const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
  const ptrdiff_t ny = std::max<ptrdiff_t>(comm->size(), 3);
  const ptrdiff_t nz = std::max<ptrdiff_t>(comm->size(), 2);
  const geom_t xmin = 0, ymin = 0, zmin = 0;
  const geom_t xmax = 1, ymax = 1, zmax = 1;

  // Build the reference on disk via the existing serial creator + writer so
  // any breakage in either path is caught here too.
  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_hex8_cube_create_distributed_%d_%d", comm->size(),
                token);
  const Path mesh_path(path_buffer);

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    create_directory(mesh_path);

    auto serial_mesh = Mesh::create_hex8_cube(
        Communicator::self(), nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SMESH_TEST_ASSERT(serial_mesh->write(mesh_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  // Reference distributed read from disk
  int file_nnodesxelem = 0;
  ptrdiff_t file_n_local_e = 0, file_n_global_e = 0;
  idx_t **file_elems = nullptr;
  SMESH_TEST_ASSERT(mesh_block_from_folder<idx_t>(
                        comm->get(), mesh_path, &file_nnodesxelem, &file_elems,
                        &file_n_local_e, &file_n_global_e) == SMESH_SUCCESS);

  int file_spatial_dim = 0;
  ptrdiff_t file_n_local_n = 0, file_n_global_n = 0;
  geom_t **file_points = nullptr;
  SMESH_TEST_ASSERT(mesh_coordinates_from_folder<geom_t>(
                        comm->get(), mesh_path, &file_spatial_dim, &file_points,
                        &file_n_local_n, &file_n_global_n) == SMESH_SUCCESS);

  // In-memory generator
  int gen_nnodesxelem = 0;
  ptrdiff_t gen_n_local_e = 0, gen_n_global_e = 0;
  idx_t **gen_elems = nullptr;

  int gen_spatial_dim = 0;
  ptrdiff_t gen_n_local_n = 0, gen_n_global_n = 0;
  geom_t **gen_points = nullptr;

  SMESH_TEST_ASSERT(
      (hex8_cube_create_distributed<idx_t, geom_t>(
           comm->get(), nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax,
           &gen_nnodesxelem, &gen_n_local_e, &gen_n_global_e, &gen_elems,
           &gen_spatial_dim, &gen_n_local_n, &gen_n_global_n,
           &gen_points) == SMESH_SUCCESS));

  // Sizes must agree (ranks land on the same chunk convention)
  SMESH_TEST_EQ(gen_nnodesxelem, file_nnodesxelem);
  SMESH_TEST_EQ(gen_n_local_e, file_n_local_e);
  SMESH_TEST_EQ(gen_n_global_e, file_n_global_e);
  SMESH_TEST_EQ(gen_spatial_dim, file_spatial_dim);
  SMESH_TEST_EQ(gen_n_local_n, file_n_local_n);
  SMESH_TEST_EQ(gen_n_global_n, file_n_global_n);

  // Connectivity bit-equal
  for (int v = 0; v < gen_nnodesxelem; ++v) {
    for (ptrdiff_t i = 0; i < gen_n_local_e; ++i) {
      SMESH_TEST_EQ(gen_elems[v][i], file_elems[v][i]);
    }
  }

  // Coordinates bit-equal
  for (int d = 0; d < gen_spatial_dim; ++d) {
    for (ptrdiff_t i = 0; i < gen_n_local_n; ++i) {
      SMESH_TEST_EQ(gen_points[d][i], file_points[d][i]);
    }
  }

  // Cleanup
  for (int v = 0; v < file_nnodesxelem; ++v) {
    SMESH_FREE(file_elems[v]);
  }
  SMESH_FREE(file_elems);

  for (int d = 0; d < file_spatial_dim; ++d) {
    SMESH_FREE(file_points[d]);
  }
  SMESH_FREE(file_points);

  for (int v = 0; v < gen_nnodesxelem; ++v) {
    SMESH_FREE(gen_elems[v]);
  }
  SMESH_FREE(gen_elems);

  for (int d = 0; d < gen_spatial_dim; ++d) {
    SMESH_FREE(gen_points[d]);
  }
  SMESH_FREE(gen_points);

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
  SMESH_RUN_TEST(test_hex8_cube_create_distributed_matches_file);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
