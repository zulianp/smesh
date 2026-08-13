#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <limits>
#include <vector>

#include "smesh_distributed_base.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sfc.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static int test_mpi_multiblock_checkerboard_read() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();
  if (comm->size() < 2) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr));
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mpi_multiblock_checkerboard_%d_%d", comm->size(),
                token);
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx =
      std::max<ptrdiff_t>(2 * ((comm->size() + 1) / 2) * 2, 4);
  const ptrdiff_t ny =
      std::max<ptrdiff_t>(2 * ((comm->size() + 1) / 2) * 2, 4);
  const ptrdiff_t nz = 2;

  ptrdiff_t serial_white = 0;
  ptrdiff_t serial_black = 0;
  ptrdiff_t serial_nodes = 0;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial = Mesh::create_hex8_checkerboard_cube(
        Communicator::self(), nx, ny, nz, 0, 0, 0, 1, 1, 1);
    SMESH_TEST_ASSERT(serial != nullptr);
    SMESH_TEST_EQ(static_cast<int>(serial->n_blocks()), 2);
    serial_white = serial->n_elements(0);
    serial_black = serial->n_elements(1);
    serial_nodes = serial->n_nodes();
    SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
  }

  comm->broadcast(&serial_white, 1, 0);
  comm->broadcast(&serial_black, 1, 0);
  comm->broadcast(&serial_nodes, 1, 0);
  comm->barrier();

  // Multi-block MPI defaults to Hilbert SFC when SMESH_REORDER is unset.
  unsetenv("SMESH_REORDER");

  Mesh mesh(comm);
  SMESH_TEST_ASSERT(mesh.read(mesh_path) == SMESH_SUCCESS);
  SMESH_TEST_EQ(static_cast<int>(mesh.n_blocks()), 2);
  SMESH_TEST_ASSERT(mesh.block(0)->name() == "white");
  SMESH_TEST_ASSERT(mesh.block(1)->name() == "black");
  SMESH_TEST_EQ(mesh.block(0)->element_type(), HEX8);
  SMESH_TEST_EQ(mesh.block(1)->element_type(), HEX8);

  auto dist = mesh.distributed();
  SMESH_TEST_ASSERT(dist != nullptr);
  SMESH_TEST_EQ(dist->n_nodes_global(), serial_nodes);
  SMESH_TEST_EQ(dist->n_elements_global(), serial_white + serial_black);
  SMESH_TEST_EQ(mesh.n_nodes(), dist->n_nodes_local());

  const auto element_mapping = dist->element_mapping();
  const ptrdiff_t n_owned = dist->n_elements_owned();
  ptrdiff_t local_owned_white = 0;
  ptrdiff_t local_owned_black = 0;
  for (ptrdiff_t i = 0; i < n_owned; ++i) {
    const large_idx_t cid = element_mapping->data()[i];
    if (cid < static_cast<large_idx_t>(serial_white)) {
      local_owned_white++;
    } else {
      local_owned_black++;
    }
  }

  ptrdiff_t global_owned_white = 0;
  ptrdiff_t global_owned_black = 0;
  MPI_Allreduce(&local_owned_white, &global_owned_white, 1,
                mpi_type<ptrdiff_t>(), MPI_SUM, comm->get());
  MPI_Allreduce(&local_owned_black, &global_owned_black, 1,
                mpi_type<ptrdiff_t>(), MPI_SUM, comm->get());
  SMESH_TEST_EQ(global_owned_white, serial_white);
  SMESH_TEST_EQ(global_owned_black, serial_black);

  // Hilbert keys of owned centroids: max on rank r <= min on rank r+1.
  {
    constexpr int nxe = 8;
    const auto points = mesh.points();
    auto white_elems = mesh.elements(0);
    auto black_elems = mesh.elements(1);

    std::vector<geom_t> cx((size_t)std::max<ptrdiff_t>(n_owned, 1));
    std::vector<geom_t> cy((size_t)std::max<ptrdiff_t>(n_owned, 1));
    std::vector<geom_t> cz((size_t)std::max<ptrdiff_t>(n_owned, 1));
    geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                           std::numeric_limits<geom_t>::max(),
                           std::numeric_limits<geom_t>::max()};
    geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                           std::numeric_limits<geom_t>::lowest(),
                           std::numeric_limits<geom_t>::lowest()};

    ptrdiff_t cursor_w = 0;
    ptrdiff_t cursor_b = 0;
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
      const large_idx_t cid = element_mapping->data()[i];
      idx_t *const *src = nullptr;
      ptrdiff_t local_e = 0;
      if (cid < static_cast<large_idx_t>(serial_white)) {
        src = white_elems->data();
        local_e = cursor_w++;
      } else {
        src = black_elems->data();
        local_e = cursor_b++;
      }

      geom_t x = 0, y = 0, z = 0;
      for (int d = 0; d < nxe; ++d) {
        const idx_t node = src[d][local_e];
        x += points->data()[0][node];
        y += points->data()[1][node];
        z += points->data()[2][node];
      }
      x /= (geom_t)nxe;
      y /= (geom_t)nxe;
      z /= (geom_t)nxe;
      cx[(size_t)i] = x;
      cy[(size_t)i] = y;
      cz[(size_t)i] = z;
      local_min[0] = std::min(local_min[0], x);
      local_min[1] = std::min(local_min[1], y);
      local_min[2] = std::min(local_min[2], z);
      local_max[0] = std::max(local_max[0], x);
      local_max[1] = std::max(local_max[1], y);
      local_max[2] = std::max(local_max[2], z);
    }
    SMESH_TEST_EQ(cursor_w, local_owned_white);
    SMESH_TEST_EQ(cursor_b, local_owned_black);

    geom_t global_min[3];
    geom_t global_max[3];
    MPI_Allreduce(local_min, global_min, 3, mpi_type<geom_t>(), MPI_MIN,
                  comm->get());
    MPI_Allreduce(local_max, global_max, 3, mpi_type<geom_t>(), MPI_MAX,
                  comm->get());

    std::vector<u32> keys((size_t)std::max<ptrdiff_t>(n_owned, 1));
    u32 local_key_min = std::numeric_limits<u32>::max();
    u32 local_key_max = 0;
    if (n_owned > 0) {
      SMESH_TEST_ASSERT(
          encode_hilbert3<geom_t>(n_owned, cx.data(), cy.data(), cz.data(),
                                  global_min[0], global_max[0], global_min[1],
                                  global_max[1], global_min[2], global_max[2],
                                  keys.data()) == SMESH_SUCCESS);
      for (ptrdiff_t e = 0; e < n_owned; ++e) {
        local_key_min = std::min(local_key_min, keys[(size_t)e]);
        local_key_max = std::max(local_key_max, keys[(size_t)e]);
      }
    }

    if (comm->rank() + 1 < comm->size() && n_owned > 0) {
      u32 next_min = 0;
      MPI_Sendrecv(&local_key_max, 1, mpi_type<u32>(), comm->rank() + 1, 0,
                   &next_min, 1, mpi_type<u32>(), comm->rank() + 1, 1,
                   comm->get(), MPI_STATUS_IGNORE);
      SMESH_TEST_ASSERT(local_key_max <= next_min);
    }

    if (comm->rank() > 0 && n_owned > 0) {
      u32 prev_max = 0;
      MPI_Sendrecv(&local_key_min, 1, mpi_type<u32>(), comm->rank() - 1, 1,
                   &prev_max, 1, mpi_type<u32>(), comm->rank() - 1, 0,
                   comm->get(), MPI_STATUS_IGNORE);
      SMESH_TEST_ASSERT(prev_max <= local_key_min);
    }
  }

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_single_block_read_regression() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();
  if (comm->size() < 2) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr)) + 17;
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mpi_single_block_regression_%d_%d", comm->size(),
                token);
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
  const ptrdiff_t ny = 3;
  const ptrdiff_t nz = 2;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz, 0, 0,
                                         0, 1, 1, 1);
    SMESH_TEST_ASSERT(serial != nullptr);
    SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
  }
  comm->barrier();

  unsetenv("SMESH_REORDER");
  Mesh mesh(comm);
  SMESH_TEST_ASSERT(mesh.read(mesh_path) == SMESH_SUCCESS);
  SMESH_TEST_EQ(static_cast<int>(mesh.n_blocks()), 1);
  SMESH_TEST_ASSERT(mesh.block(0)->name() == "default");
  SMESH_TEST_EQ(mesh.block(0)->element_type(), HEX8);
  SMESH_TEST_ASSERT(mesh.distributed() != nullptr);

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_mpi_multiblock_checkerboard_read);
  SMESH_RUN_TEST(test_mpi_single_block_read_regression);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
