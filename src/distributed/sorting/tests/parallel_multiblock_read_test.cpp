#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <limits>
#include <vector>

#include "smesh_conversion.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sfc.hpp"
#include "smesh_test.hpp"

using namespace smesh;

#ifdef SMESH_ENABLE_MPI
static int check_per_block_layout(const Mesh &mesh,
                                  const ptrdiff_t *const serial_n_elements) {
  auto dist = mesh.distributed();
  SMESH_TEST_ASSERT(dist != nullptr);
  const int rank = mesh.comm()->rank();
  auto node_owner_buf = dist->node_owner();
  SMESH_TEST_ASSERT(node_owner_buf != nullptr);
  const int *const node_owner = node_owner_buf->data();
  const ptrdiff_t n_local_nodes = mesh.n_nodes();

  ptrdiff_t sum_owned = 0;
  ptrdiff_t sum_shared = 0;
  ptrdiff_t sum_ghosts = 0;

  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    auto block = mesh.block((block_idx_t)b);
    const ptrdiff_t n_owned = block->n_elements_owned();
    const ptrdiff_t n_shared = block->n_elements_shared();
    const ptrdiff_t n_ghosts = block->n_elements_ghosts();
    const ptrdiff_t n_ons = block->n_elements_owned_not_shared();
    SMESH_TEST_EQ(block->n_elements(), n_owned + n_ghosts);
    SMESH_TEST_EQ(n_owned, n_ons + n_shared);
    SMESH_TEST_ASSERT(n_ons >= 0);

    ptrdiff_t global_owned = 0;
    MPI_Allreduce(&n_owned, &global_owned, 1, mpi_type<ptrdiff_t>(), MPI_SUM,
                  mesh.comm()->get());
    SMESH_TEST_EQ(global_owned, serial_n_elements[b]);

    auto elems = block->elements();
    SMESH_TEST_ASSERT(elems != nullptr);
    const int nxe = block->n_nodes_per_element();
    for (ptrdiff_t e = 0; e < n_ons; ++e) {
      for (int d = 0; d < nxe; ++d) {
        const idx_t node = elems->data()[d][e];
        SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(node) >= 0);
        SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(node) < n_local_nodes);
        SMESH_TEST_EQ(node_owner[node], rank);
      }
    }
    for (ptrdiff_t e = n_ons; e < n_owned; ++e) {
      bool has_remote = false;
      for (int d = 0; d < nxe; ++d) {
        const idx_t node = elems->data()[d][e];
        SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(node) >= 0);
        SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(node) < n_local_nodes);
        has_remote |= (node_owner[node] != rank);
      }
      SMESH_TEST_ASSERT(has_remote);
    }

    auto owned_map = block->element_mapping();
    auto aura_map = block->aura_element_mapping();
    if (n_owned > 0) {
      SMESH_TEST_ASSERT(owned_map != nullptr);
      SMESH_TEST_EQ(static_cast<ptrdiff_t>(owned_map->size()), n_owned);
      for (ptrdiff_t i = 0; i < n_owned; ++i) {
        const large_idx_t id = owned_map->data()[i];
        SMESH_TEST_ASSERT(id >= 0);
        SMESH_TEST_ASSERT(id < static_cast<large_idx_t>(serial_n_elements[b]));
      }
    }
    if (n_ghosts > 0) {
      SMESH_TEST_ASSERT(aura_map != nullptr);
      SMESH_TEST_EQ(static_cast<ptrdiff_t>(aura_map->size()), n_ghosts);
      for (ptrdiff_t i = 0; i < n_ghosts; ++i) {
        const large_idx_t id = aura_map->data()[i];
        SMESH_TEST_ASSERT(id >= 0);
        SMESH_TEST_ASSERT(id < static_cast<large_idx_t>(serial_n_elements[b]));
      }
    }

    sum_owned += n_owned;
    sum_shared += n_shared;
    sum_ghosts += n_ghosts;
  }

  SMESH_TEST_EQ(sum_owned, dist->n_elements_owned());
  SMESH_TEST_EQ(sum_shared, dist->n_elements_shared());
  SMESH_TEST_EQ(sum_ghosts, dist->n_elements_ghosts());
  return SMESH_TEST_SUCCESS;
}

static int check_single_block_matches_distributed(const Mesh &mesh) {
  SMESH_TEST_EQ(static_cast<int>(mesh.n_blocks()), 1);
  auto dist = mesh.distributed();
  SMESH_TEST_ASSERT(dist != nullptr);
  auto block = mesh.block(0);
  SMESH_TEST_EQ(block->n_elements_owned(), dist->n_elements_owned());
  SMESH_TEST_EQ(block->n_elements_shared(), dist->n_elements_shared());
  SMESH_TEST_EQ(block->n_elements_ghosts(), dist->n_elements_ghosts());
  SMESH_TEST_EQ(block->n_elements_owned_not_shared(),
                dist->n_elements_owned_not_shared());
  SMESH_TEST_EQ(block->n_elements(), dist->n_elements_local());

  const ptrdiff_t n_owned = dist->n_elements_owned();
  const ptrdiff_t n_ghosts = dist->n_elements_ghosts();
  if (n_owned > 0) {
    SMESH_TEST_ASSERT(block->element_mapping() != nullptr);
    SMESH_TEST_ASSERT(dist->element_mapping() != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(block->element_mapping()->size()),
                  n_owned);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(dist->element_mapping()->size()),
                  n_owned);
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
      SMESH_TEST_EQ(block->element_mapping()->data()[i],
                    dist->element_mapping()->data()[i]);
    }
  }
  if (n_ghosts > 0) {
    SMESH_TEST_ASSERT(block->aura_element_mapping() != nullptr);
    SMESH_TEST_ASSERT(dist->aura_element_mapping() != nullptr);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(block->aura_element_mapping()->size()),
                  n_ghosts);
    SMESH_TEST_EQ(static_cast<ptrdiff_t>(dist->aura_element_mapping()->size()),
                  n_ghosts);
    for (ptrdiff_t i = 0; i < n_ghosts; ++i) {
      SMESH_TEST_EQ(block->aura_element_mapping()->data()[i],
                    dist->aura_element_mapping()->data()[i]);
    }
  }
  return SMESH_TEST_SUCCESS;
}

static int mpi_write_read_roundtrip(const Mesh &mesh, const Path &out_path,
                                    const ptrdiff_t *const serial_n_elements,
                                    const size_t n_blocks_ref) {
  auto comm = mesh.comm();

  if (comm->rank() == 0) {
    std::filesystem::remove_all(out_path.to_string());
  }
  comm->barrier();

  SMESH_TEST_ASSERT(mesh.write(out_path) == SMESH_SUCCESS);
  comm->barrier();

  if (comm->rank() == 0) {
    Mesh serial_check(Communicator::self());
    SMESH_TEST_ASSERT(serial_check.read(out_path) == SMESH_SUCCESS);
    SMESH_TEST_EQ(serial_check.n_blocks(), n_blocks_ref);
    for (size_t b = 0; b < n_blocks_ref; ++b) {
      SMESH_TEST_EQ(serial_check.n_elements((block_idx_t)b),
                    serial_n_elements[b]);
      SMESH_TEST_ASSERT(serial_check.block((block_idx_t)b)->name() ==
                        mesh.block((block_idx_t)b)->name());
      SMESH_TEST_EQ(serial_check.block((block_idx_t)b)->element_type(),
                    mesh.block((block_idx_t)b)->element_type());
    }
  }
  comm->barrier();

  Mesh roundtrip(comm);
  SMESH_TEST_ASSERT(roundtrip.read(out_path) == SMESH_SUCCESS);
  SMESH_TEST_EQ(roundtrip.n_blocks(), n_blocks_ref);
  for (size_t b = 0; b < n_blocks_ref; ++b) {
    SMESH_TEST_ASSERT(roundtrip.block((block_idx_t)b)->name() ==
                      mesh.block((block_idx_t)b)->name());
    SMESH_TEST_EQ(roundtrip.block((block_idx_t)b)->element_type(),
                  mesh.block((block_idx_t)b)->element_type());
  }
  SMESH_TEST_EQ(check_per_block_layout(roundtrip, serial_n_elements),
                SMESH_TEST_SUCCESS);

  if (comm->rank() == 0) {
    std::filesystem::remove_all(out_path.to_string());
  }
  comm->barrier();
  return SMESH_TEST_SUCCESS;
}
#endif

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
  SMESH_TEST_EQ(mesh.block(0)->n_elements_owned(), local_owned_white);
  SMESH_TEST_EQ(mesh.block(1)->n_elements_owned(), local_owned_black);

  {
    const ptrdiff_t serial_n[2] = {serial_white, serial_black};
    SMESH_TEST_EQ(check_per_block_layout(mesh, serial_n), SMESH_TEST_SUCCESS);
  }

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

  {
    char rt_path_buffer[256];
    std::snprintf(rt_path_buffer, sizeof(rt_path_buffer),
                  "/tmp/smesh_mpi_multiblock_checkerboard_rt_%d_%d",
                  comm->size(), token + 1000);
    const Path rt_path(rt_path_buffer);
    const ptrdiff_t serial_n[2] = {serial_white, serial_black};
    SMESH_TEST_EQ(mpi_write_read_roundtrip(mesh, rt_path, serial_n, 2),
                  SMESH_TEST_SUCCESS);
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

  ptrdiff_t serial_elements = 0;
  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz, 0, 0,
                                         0, 1, 1, 1);
    SMESH_TEST_ASSERT(serial != nullptr);
    serial_elements = serial->n_elements();
    SMESH_TEST_ASSERT(serial->write(mesh_path) == SMESH_SUCCESS);
  }
  comm->broadcast(&serial_elements, 1, 0);
  comm->barrier();

  unsetenv("SMESH_REORDER");
  Mesh mesh(comm);
  SMESH_TEST_ASSERT(mesh.read(mesh_path) == SMESH_SUCCESS);
  SMESH_TEST_EQ(static_cast<int>(mesh.n_blocks()), 1);
  SMESH_TEST_ASSERT(mesh.block(0)->name() == "default");
  SMESH_TEST_EQ(mesh.block(0)->element_type(), HEX8);
  SMESH_TEST_ASSERT(mesh.distributed() != nullptr);
  SMESH_TEST_EQ(check_single_block_matches_distributed(mesh),
                SMESH_TEST_SUCCESS);
  SMESH_TEST_EQ(check_per_block_layout(mesh, &serial_elements),
                SMESH_TEST_SUCCESS);

  {
    char rt_path_buffer[256];
    std::snprintf(rt_path_buffer, sizeof(rt_path_buffer),
                  "/tmp/smesh_mpi_single_block_rt_%d_%d", comm->size(),
                  token + 1000);
    const Path rt_path(rt_path_buffer);
    SMESH_TEST_EQ(mpi_write_read_roundtrip(mesh, rt_path, &serial_elements, 1),
                  SMESH_TEST_SUCCESS);
  }

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
  }
  return SMESH_TEST_SUCCESS;
#endif
}

static int test_mpi_multiblock_hex8_tet4_read() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();
  if (comm->size() < 2) {
    return SMESH_TEST_SUCCESS;
  }

  int token = 0;
  if (comm->rank() == 0) {
    token = static_cast<int>(std::time(nullptr)) + 31;
  }
  comm->broadcast(&token, 1, 0);

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_mpi_multiblock_hex8_tet4_%d_%d", comm->size(),
                token);
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 4);
  const ptrdiff_t ny = 2;
  const ptrdiff_t nz = 2;

  ptrdiff_t serial_hex = 0;
  ptrdiff_t serial_tet = 0;
  ptrdiff_t serial_nodes = 0;

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto cube = Mesh::create_hex8_cube(Communicator::self(), nx, ny, nz, 0, 0,
                                       0, 1, 1, 1);
    SMESH_TEST_ASSERT(cube != nullptr);
    const ptrdiff_t n_hex_all = cube->n_elements();
    SMESH_TEST_ASSERT(n_hex_all >= 2);
    const ptrdiff_t n_hex_keep = n_hex_all / 2;
    const ptrdiff_t n_hex_conv = n_hex_all - n_hex_keep;
    SMESH_TEST_ASSERT(n_hex_keep >= 1);
    SMESH_TEST_ASSERT(n_hex_conv >= 1);

    auto hex_src = cube->elements(0)->data();
    auto hex_keep = create_host_buffer<idx_t>(8, (size_t)n_hex_keep);
    for (int d = 0; d < 8; ++d) {
      std::memcpy(hex_keep->data()[d], hex_src[d],
                  (size_t)n_hex_keep * sizeof(idx_t));
    }

    idx_t *hex_tail[8];
    for (int d = 0; d < 8; ++d) {
      hex_tail[d] = hex_src[d] + n_hex_keep;
    }
    auto tet_buf = create_host_buffer<idx_t>(4, (size_t)(n_hex_conv * 6));
    mesh_hex8_to_6x_tet4<idx_t>(n_hex_conv, hex_tail, tet_buf->data());

    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(std::make_shared<Mesh::Block>("hex", HEX8, hex_keep));
    blocks.push_back(std::make_shared<Mesh::Block>("tet", TET4, tet_buf));
    Mesh serial(Communicator::self(), blocks, cube->points());
    serial_hex = n_hex_keep;
    serial_tet = n_hex_conv * 6;
    serial_nodes = serial.n_nodes();
    SMESH_TEST_ASSERT(serial.write(mesh_path) == SMESH_SUCCESS);
  }

  comm->broadcast(&serial_hex, 1, 0);
  comm->broadcast(&serial_tet, 1, 0);
  comm->broadcast(&serial_nodes, 1, 0);
  comm->barrier();

  unsetenv("SMESH_REORDER");

  Mesh mesh(comm);
  SMESH_TEST_ASSERT(mesh.read(mesh_path) == SMESH_SUCCESS);
  SMESH_TEST_EQ(static_cast<int>(mesh.n_blocks()), 2);
  SMESH_TEST_ASSERT(mesh.block(0)->name() == "hex");
  SMESH_TEST_ASSERT(mesh.block(1)->name() == "tet");
  SMESH_TEST_EQ(mesh.block(0)->element_type(), HEX8);
  SMESH_TEST_EQ(mesh.block(1)->element_type(), TET4);
  SMESH_TEST_EQ(mesh.n_nodes_per_element(0), 8);
  SMESH_TEST_EQ(mesh.n_nodes_per_element(1), 4);

  auto dist = mesh.distributed();
  SMESH_TEST_ASSERT(dist != nullptr);
  SMESH_TEST_EQ(dist->n_nodes_global(), serial_nodes);
  SMESH_TEST_EQ(dist->n_elements_global(), serial_hex + serial_tet);
  SMESH_TEST_EQ(mesh.n_nodes(), dist->n_nodes_local());

  const ptrdiff_t n_local_nodes = mesh.n_nodes();
  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    const int nxe = mesh.n_nodes_per_element((block_idx_t)b);
    auto elems = mesh.elements((block_idx_t)b);
    const ptrdiff_t n_e = mesh.n_elements((block_idx_t)b);
    for (ptrdiff_t e = 0; e < n_e; ++e) {
      for (int d = 0; d < nxe; ++d) {
        const idx_t node = elems->data()[d][e];
        SMESH_TEST_ASSERT(node != invalid_idx<idx_t>());
        SMESH_TEST_ASSERT(static_cast<ptrdiff_t>(node) < n_local_nodes);
      }
    }
  }

  const auto element_mapping = dist->element_mapping();
  const ptrdiff_t n_owned = dist->n_elements_owned();
  ptrdiff_t local_owned_hex = 0;
  ptrdiff_t local_owned_tet = 0;
  for (ptrdiff_t i = 0; i < n_owned; ++i) {
    const large_idx_t cid = element_mapping->data()[i];
    if (cid < static_cast<large_idx_t>(serial_hex)) {
      local_owned_hex++;
    } else {
      local_owned_tet++;
    }
  }

  ptrdiff_t global_owned_hex = 0;
  ptrdiff_t global_owned_tet = 0;
  MPI_Allreduce(&local_owned_hex, &global_owned_hex, 1, mpi_type<ptrdiff_t>(),
                MPI_SUM, comm->get());
  MPI_Allreduce(&local_owned_tet, &global_owned_tet, 1, mpi_type<ptrdiff_t>(),
                MPI_SUM, comm->get());
  SMESH_TEST_EQ(global_owned_hex, serial_hex);
  SMESH_TEST_EQ(global_owned_tet, serial_tet);
  SMESH_TEST_EQ(mesh.block(0)->n_elements_owned(), local_owned_hex);
  SMESH_TEST_EQ(mesh.block(1)->n_elements_owned(), local_owned_tet);

  {
    const ptrdiff_t serial_n[2] = {serial_hex, serial_tet};
    SMESH_TEST_EQ(check_per_block_layout(mesh, serial_n), SMESH_TEST_SUCCESS);
  }

  ptrdiff_t n_owned_nodes = dist->n_nodes_owned();
  ptrdiff_t global_owned_nodes = 0;
  MPI_Allreduce(&n_owned_nodes, &global_owned_nodes, 1, mpi_type<ptrdiff_t>(),
                MPI_SUM, comm->get());
  SMESH_TEST_EQ(global_owned_nodes, serial_nodes);

  {
    const auto points = mesh.points();
    auto hex_elems = mesh.elements(0);
    auto tet_elems = mesh.elements(1);

    std::vector<geom_t> cx((size_t)std::max<ptrdiff_t>(n_owned, 1));
    std::vector<geom_t> cy((size_t)std::max<ptrdiff_t>(n_owned, 1));
    std::vector<geom_t> cz((size_t)std::max<ptrdiff_t>(n_owned, 1));
    geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                           std::numeric_limits<geom_t>::max(),
                           std::numeric_limits<geom_t>::max()};
    geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                           std::numeric_limits<geom_t>::lowest(),
                           std::numeric_limits<geom_t>::lowest()};

    ptrdiff_t cursor_h = 0;
    ptrdiff_t cursor_t = 0;
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
      const large_idx_t cid = element_mapping->data()[i];
      idx_t *const *src = nullptr;
      ptrdiff_t local_e = 0;
      int nxe = 0;
      if (cid < static_cast<large_idx_t>(serial_hex)) {
        src = hex_elems->data();
        local_e = cursor_h++;
        nxe = 8;
      } else {
        src = tet_elems->data();
        local_e = cursor_t++;
        nxe = 4;
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
    SMESH_TEST_EQ(cursor_h, local_owned_hex);
    SMESH_TEST_EQ(cursor_t, local_owned_tet);

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

  {
    char rt_path_buffer[256];
    std::snprintf(rt_path_buffer, sizeof(rt_path_buffer),
                  "/tmp/smesh_mpi_multiblock_hex8_tet4_rt_%d_%d",
                  comm->size(), token + 1000);
    const Path rt_path(rt_path_buffer);
    const ptrdiff_t serial_n[2] = {serial_hex, serial_tet};
    SMESH_TEST_EQ(mpi_write_read_roundtrip(mesh, rt_path, serial_n, 2),
                  SMESH_TEST_SUCCESS);
  }

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
  SMESH_RUN_TEST(test_mpi_multiblock_hex8_tet4_read);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
