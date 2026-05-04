#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <random>
#include <vector>

#include "smesh_distributed_read.hpp"
#include "smesh_distributed_write.hpp"
#include "smesh_decompose.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_read.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static void permute_owned_elements(const ptrdiff_t n_owned_elements,
                                   const int nnodesxelem, idx_t **elements,
                                   large_idx_t *element_mapping,
                                   const unsigned int seed) {
  if (n_owned_elements <= 1) {
    return;
  }

  std::vector<ptrdiff_t> perm((size_t)n_owned_elements);
  for (ptrdiff_t i = 0; i < n_owned_elements; ++i) {
    perm[(size_t)i] = i;
  }

  std::mt19937 gen(seed);
  std::shuffle(perm.begin(), perm.end(), gen);

  std::vector<large_idx_t> mapping_copy((size_t)n_owned_elements);
  for (ptrdiff_t i = 0; i < n_owned_elements; ++i) {
    mapping_copy[(size_t)i] = element_mapping[i];
  }

  for (ptrdiff_t i = 0; i < n_owned_elements; ++i) {
    element_mapping[i] = mapping_copy[(size_t)perm[(size_t)i]];
  }

  for (int d = 0; d < nnodesxelem; ++d) {
    std::vector<idx_t> row_copy((size_t)n_owned_elements);
    for (ptrdiff_t i = 0; i < n_owned_elements; ++i) {
      row_copy[(size_t)i] = elements[d][i];
    }

    for (ptrdiff_t i = 0; i < n_owned_elements; ++i) {
      elements[d][i] = row_copy[(size_t)perm[(size_t)i]];
    }
  }
}

static int compare_mesh_folders_serial(const Path &input_path,
                                       const Path &output_path) {
  int in_nnodesxelem = 0, out_nnodesxelem = 0;
  ptrdiff_t in_nelems = 0, out_nelems = 0;
  idx_t **in_elements = nullptr, **out_elements = nullptr;

  int in_spatial_dim = 0, out_spatial_dim = 0;
  ptrdiff_t in_nnodes = 0, out_nnodes = 0;
  geom_t **in_points = nullptr, **out_points = nullptr;

  const int in_read =
      mesh_from_folder(input_path, &in_nnodesxelem, &in_nelems, &in_elements,
                       &in_spatial_dim, &in_nnodes, &in_points);
  const int out_read = mesh_from_folder(
      output_path, &out_nnodesxelem, &out_nelems, &out_elements,
      &out_spatial_dim, &out_nnodes, &out_points);

  SMESH_TEST_ASSERT(in_read == SMESH_SUCCESS);
  SMESH_TEST_ASSERT(out_read == SMESH_SUCCESS);

  SMESH_TEST_EQ(in_nnodesxelem, out_nnodesxelem);
  SMESH_TEST_EQ(in_nelems, out_nelems);
  SMESH_TEST_EQ(in_spatial_dim, out_spatial_dim);
  SMESH_TEST_EQ(in_nnodes, out_nnodes);

  for (int d = 0; d < in_nnodesxelem; ++d) {
    for (ptrdiff_t i = 0; i < in_nelems; ++i) {
      SMESH_TEST_EQ(in_elements[d][i], out_elements[d][i]);
    }
  }

  for (int d = 0; d < in_spatial_dim; ++d) {
    for (ptrdiff_t i = 0; i < in_nnodes; ++i) {
      SMESH_TEST_EQ(in_points[d][i], out_points[d][i]);
    }
  }

  for (int d = 0; d < in_nnodesxelem; ++d) {
    SMESH_FREE(in_elements[d]);
    SMESH_FREE(out_elements[d]);
  }
  SMESH_FREE(in_elements);
  SMESH_FREE(out_elements);

  for (int d = 0; d < in_spatial_dim; ++d) {
    SMESH_FREE(in_points[d]);
    SMESH_FREE(out_points[d]);
  }
  SMESH_FREE(in_points);
  SMESH_FREE(out_points);

  return SMESH_TEST_SUCCESS;
}

int test_random_distribution_roundtrip() {
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
                "/tmp/smesh_random_distribution_roundtrip_%d_%d", comm->size(),
                token);

  const Path root(path_buffer);
  const Path input_path = root / "input_mesh";
  const Path output_path = root / "output_mesh";

  if (comm->rank() == 0) {
    std::filesystem::remove_all(root.to_string());
    create_directory(root);
    create_directory(input_path);
    create_directory(output_path);

    const ptrdiff_t nx = std::max<ptrdiff_t>(2 * comm->size(), 8);
    auto serial_mesh = Mesh::create_tet4_cube(Communicator::self(), nx, 4, 3);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SMESH_TEST_ASSERT(serial_mesh->write(input_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  int nnodesxelem = 0;
  ptrdiff_t n_local_elements = 0;
  ptrdiff_t n_global_elements = 0;
  idx_t **raw_elements = nullptr;
  SMESH_TEST_ASSERT(mesh_block_from_folder<idx_t>(
                        comm->get(), input_path, &nnodesxelem, &raw_elements,
                        &n_local_elements, &n_global_elements) == SMESH_SUCCESS);

  int spatial_dim = 0;
  ptrdiff_t n_local2global = 0;
  ptrdiff_t n_global_nodes = 0;
  geom_t **raw_points = nullptr;
  SMESH_TEST_ASSERT(mesh_coordinates_from_folder<geom_t>(
                        comm->get(), input_path, &spatial_dim, &raw_points,
                        &n_local2global, &n_global_nodes) == SMESH_SUCCESS);

  auto *input_element_mapping =
      (large_idx_t *)SMESH_ALLOC(n_local_elements * sizeof(large_idx_t));
  const ptrdiff_t element_start =
      rank_start(n_global_elements, comm->size(), comm->rank());
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    input_element_mapping[i] = element_start + i;
  }

  const unsigned int local_seed =
      static_cast<unsigned int>(token + 7919 * comm->rank());
  permute_owned_elements(n_local_elements, nnodesxelem, raw_elements,
                         input_element_mapping, local_seed);

  ptrdiff_t n_owned_elements = 0;
  ptrdiff_t n_shared_elements = 0;
  ptrdiff_t n_ghost_elements = 0;
  large_idx_t *element_mapping = nullptr;
  large_idx_t *aura_element_mapping = nullptr;
  idx_t **elements = nullptr;

  ptrdiff_t n_owned_nodes = 0;
  ptrdiff_t n_shared_nodes = 0;
  ptrdiff_t n_ghost_nodes = 0;
  ptrdiff_t n_aura_nodes = 0;
  large_idx_t *node_mapping = nullptr;
  geom_t **points = nullptr;

  int *node_owner = nullptr;
  ptrdiff_t *node_offsets = nullptr;
  idx_t *ghosts = nullptr;

  SMESH_TEST_ASSERT((mesh_create_parallel<idx_t, geom_t, large_idx_t>(
                        comm->get(), comm->size(), comm->rank(), nnodesxelem,
                        raw_elements, n_local_elements, n_global_elements,
                        spatial_dim, raw_points, n_local2global, n_global_nodes,
                        input_element_mapping, &nnodesxelem, &n_global_elements,
                        &n_owned_elements, &n_shared_elements, &n_ghost_elements,
                        &element_mapping, &aura_element_mapping, &elements,
                        &spatial_dim, &n_global_nodes, &n_owned_nodes,
                        &n_shared_nodes, &n_ghost_nodes, &n_aura_nodes,
                        &node_mapping, &points, &node_owner, &node_offsets,
                        &ghosts) == SMESH_SUCCESS));
  SMESH_FREE(input_element_mapping);

  const enum ElemType element_type = static_cast<enum ElemType>(nnodesxelem);
  SMESH_TEST_ASSERT(write_distributed_mesh_topology(
                        comm->get(), output_path, element_type, spatial_dim,
                        n_global_elements, n_owned_elements, element_mapping,
                        nnodesxelem, elements, n_global_nodes, n_owned_nodes,
                        node_mapping, points) == SMESH_SUCCESS);

  comm->barrier();

  if (comm->rank() == 0) {
    SMESH_TEST_ASSERT(compare_mesh_folders_serial(input_path, output_path) ==
                      SMESH_TEST_SUCCESS);
    std::filesystem::remove_all(root.to_string());
  }

  comm->barrier();

  for (int d = 0; d < nnodesxelem; ++d) {
    SMESH_FREE(elements[d]);
  }
  SMESH_FREE(elements);

  for (int d = 0; d < spatial_dim; ++d) {
    SMESH_FREE(points[d]);
  }
  SMESH_FREE(points);

  SMESH_FREE(element_mapping);
  SMESH_FREE(aura_element_mapping);
  SMESH_FREE(node_mapping);
  SMESH_FREE(node_owner);
  SMESH_FREE(node_offsets);
  SMESH_FREE(ghosts);

  SMESH_UNUSED(n_shared_elements);
  SMESH_UNUSED(n_ghost_elements);
  SMESH_UNUSED(n_shared_nodes);
  SMESH_UNUSED(n_ghost_nodes);
  SMESH_UNUSED(n_aura_nodes);

  return SMESH_TEST_SUCCESS;
#endif
}

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
  SMESH_RUN_TEST(test_random_distribution_roundtrip);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
