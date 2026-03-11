#include "smesh_context.hpp"
#include "smesh_exchange.hpp"
#include "smesh_output.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <limits>
#include <stdio.h>

#if defined(SMESH_ENABLE_MPI)
#include "smesh_distributed_base.hpp"
#include <mpi.h>
#endif

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("exchange_test.exe");

  auto ctx = smesh::initialize(argc, argv);

  if (argc != 2) {
    if (!ctx->communicator()->rank()) {
      fprintf(stderr, "Usage: %s <input_folder>\n", argv[0]);
    }
    return SMESH_FAILURE;
  }

  {
    auto comm = ctx->communicator();

    if (comm->size() == 1) {
      fprintf(stderr, "%s: This test is only for distributed meshes\n",
              argv[0]);
      return SMESH_FAILURE;
    }

    auto mesh = Mesh::create_from_file(comm, Path(argv[1]));
    auto ghost_exchange =
        Exchange::create_nodal(mesh, Exchange::ExchangeScope::GhostsOnly);
    auto full_exchange =
        Exchange::create_nodal(mesh, Exchange::ExchangeScope::GhostsAndAura);

    auto dist = mesh->distributed();
    const ptrdiff_t n_nodes_local = dist->n_nodes_local();
    const ptrdiff_t n_nodes_owned = dist->n_nodes_owned();
    const ptrdiff_t n_nodes_ghosts = dist->n_nodes_ghosts();
    const ptrdiff_t n_nodes_aura = dist->n_nodes_aura();
    const ptrdiff_t ghost_begin = n_nodes_owned;
    const ptrdiff_t ghost_end = n_nodes_owned + n_nodes_ghosts;
    const ptrdiff_t aura_begin = ghost_end;
    auto nodal = create_host_buffer<geom_t>(n_nodes_local);
    auto b_nodal = nodal->data();
    auto b_points = mesh->points()->data();

    auto init_for_gather = [&]() {
      for (ptrdiff_t i = 0; i < n_nodes_owned; i++) {
        b_nodal[i] = b_points[0][i];
      }
      for (ptrdiff_t i = n_nodes_owned; i < n_nodes_local; i++) {
        b_nodal[i] = std::numeric_limits<geom_t>::max();
      }
    };

    init_for_gather();
    ghost_exchange->gather(b_nodal);

    for (ptrdiff_t i = 0; i < ghost_end; i++) {
      if (b_nodal[i] != b_points[0][i]) {
        SMESH_ERROR("ghost gather mismatch at node %ld: %f != %f", i,
                    b_nodal[i], b_points[0][i]);
        return SMESH_FAILURE;
      }
    }
    for (ptrdiff_t i = aura_begin; i < n_nodes_local; i++) {
      if (b_nodal[i] != std::numeric_limits<geom_t>::max()) {
        SMESH_ERROR("ghost gather touched aura node %ld: %f", i, b_nodal[i]);
        return SMESH_FAILURE;
      }
    }

    init_for_gather();
    full_exchange->gather(b_nodal);

    for (ptrdiff_t i = 0; i < n_nodes_local; i++) {
      if (b_nodal[i] != b_points[0][i]) {
        SMESH_ERROR("full gather mismatch at node %ld: %f != %f", i, b_nodal[i],
                    b_points[0][i]);
        return SMESH_FAILURE;
      }
    }

    auto check_scatter = [&](const char *label, Exchange &exchange,
                             const bool with_aura) {
      for (ptrdiff_t i = 0; i < n_nodes_owned; i++) {
        b_nodal[i] = 0;
      }
      for (ptrdiff_t i = ghost_begin; i < ghost_end; i++) {
        b_nodal[i] = 1;
      }
      for (ptrdiff_t i = aura_begin; i < n_nodes_local; i++) {
        b_nodal[i] = with_aura ? 1 : 0;
      }

      exchange.scatter_add(b_nodal);

      double local_owned_sum = 0.0;
      double local_import_count =
          (double)n_nodes_ghosts + (with_aura ? (double)n_nodes_aura : 0.0);
      for (ptrdiff_t i = 0; i < n_nodes_owned; i++) {
        local_owned_sum += b_nodal[i];
      }

      double global_owned_sum = local_owned_sum;
      double global_import_count = local_import_count;
#if defined(SMESH_ENABLE_MPI)
      MPI_Allreduce(&local_owned_sum, &global_owned_sum, 1, MPI_DOUBLE, MPI_SUM,
                    comm->get());
      MPI_Allreduce(&local_import_count, &global_import_count, 1, MPI_DOUBLE,
                    MPI_SUM, comm->get());
#endif

      if (global_owned_sum != global_import_count) {
        SMESH_ERROR("%s scatter mismatch: %f != %f", label, global_owned_sum,
                    global_import_count);
        return SMESH_FAILURE;
      }

      // comm->print_callback([&](std::ostream &os) {
      // const ptrdiff_t n_nodes_owned_not_shared =
      // dist->n_nodes_owned_not_shared();
      //   os << label << " nodal data (shared):\n";
      //   for (ptrdiff_t i = n_nodes_owned_not_shared; i < n_nodes_owned; i++)
      //   {
      //     os << b_nodal[i] << "\n";
      //   }
      //   os << "\n";
      // });

      return SMESH_SUCCESS;
    };

    if (check_scatter("ghost-only", *ghost_exchange, false) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    if (check_scatter("ghost+aura", *full_exchange, true) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    constexpr ptrdiff_t block_size = 3;
    auto nodal_vec = create_host_buffer<geom_t>(n_nodes_local * block_size);
    auto b_nodal_vec = nodal_vec->data();

    auto init_vec_for_gather = [&]() {
      for (ptrdiff_t i = 0; i < n_nodes_owned; ++i) {
        for (ptrdiff_t c = 0; c < block_size; ++c) {
          b_nodal_vec[i * block_size + c] = b_points[0][i] + (geom_t)c;
        }
      }
      for (ptrdiff_t i = n_nodes_owned; i < n_nodes_local; ++i) {
        for (ptrdiff_t c = 0; c < block_size; ++c) {
          b_nodal_vec[i * block_size + c] = std::numeric_limits<geom_t>::max();
        }
      }
    };

    init_vec_for_gather();
    ghost_exchange->gather(b_nodal_vec, block_size);
    for (ptrdiff_t i = 0; i < ghost_end; ++i) {
      for (ptrdiff_t c = 0; c < block_size; ++c) {
        const geom_t expected = b_points[0][i] + (geom_t)c;
        const geom_t actual = b_nodal_vec[i * block_size + c];
        if (actual != expected) {
          SMESH_ERROR("ghost gather(vector) mismatch at node %ld comp %ld: %f != %f",
                      i, c, actual, expected);
          return SMESH_FAILURE;
        }
      }
    }
    for (ptrdiff_t i = aura_begin; i < n_nodes_local; ++i) {
      for (ptrdiff_t c = 0; c < block_size; ++c) {
        if (b_nodal_vec[i * block_size + c] !=
            std::numeric_limits<geom_t>::max()) {
          SMESH_ERROR("ghost gather(vector) touched aura node %ld comp %ld: %f",
                      i, c, b_nodal_vec[i * block_size + c]);
          return SMESH_FAILURE;
        }
      }
    }

    init_vec_for_gather();
    full_exchange->gather(b_nodal_vec, block_size);
    for (ptrdiff_t i = 0; i < n_nodes_local; ++i) {
      for (ptrdiff_t c = 0; c < block_size; ++c) {
        const geom_t expected = b_points[0][i] + (geom_t)c;
        const geom_t actual = b_nodal_vec[i * block_size + c];
        if (actual != expected) {
          SMESH_ERROR("full gather(vector) mismatch at node %ld comp %ld: %f != %f",
                      i, c, actual, expected);
          return SMESH_FAILURE;
        }
      }
    }

    auto check_scatter_vec = [&](const char *label, Exchange &exchange,
                                 const bool with_aura) {
      for (ptrdiff_t i = 0; i < n_nodes_owned; ++i) {
        for (ptrdiff_t c = 0; c < block_size; ++c) {
          b_nodal_vec[i * block_size + c] = 0;
        }
      }
      for (ptrdiff_t i = ghost_begin; i < ghost_end; ++i) {
        for (ptrdiff_t c = 0; c < block_size; ++c) {
          b_nodal_vec[i * block_size + c] = 1;
        }
      }
      for (ptrdiff_t i = aura_begin; i < n_nodes_local; ++i) {
        for (ptrdiff_t c = 0; c < block_size; ++c) {
          b_nodal_vec[i * block_size + c] = with_aura ? 1 : 0;
        }
      }

      exchange.scatter_add(b_nodal_vec, block_size);

      double local_owned_sum = 0.0;
      double local_import_count =
          ((double)n_nodes_ghosts + (with_aura ? (double)n_nodes_aura : 0.0)) *
          (double)block_size;
      for (ptrdiff_t i = 0; i < n_nodes_owned * block_size; ++i) {
        local_owned_sum += b_nodal_vec[i];
      }

      double global_owned_sum = local_owned_sum;
      double global_import_count = local_import_count;
#if defined(SMESH_ENABLE_MPI)
      MPI_Allreduce(&local_owned_sum, &global_owned_sum, 1, MPI_DOUBLE, MPI_SUM,
                    comm->get());
      MPI_Allreduce(&local_import_count, &global_import_count, 1, MPI_DOUBLE,
                    MPI_SUM, comm->get());
#endif
      if (global_owned_sum != global_import_count) {
        SMESH_ERROR("%s scatter(vector) mismatch: %f != %f", label,
                    global_owned_sum, global_import_count);
        return SMESH_FAILURE;
      }
      return SMESH_SUCCESS;
    };

    if (check_scatter_vec("ghost-only", *ghost_exchange, false) !=
        SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    if (check_scatter_vec("ghost+aura", *full_exchange, true) !=
        SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    constexpr ptrdiff_t wide_block_size = 5;
    auto nodal_wide = create_host_buffer<geom_t>(n_nodes_local * wide_block_size);
    auto b_nodal_wide = nodal_wide->data();

    for (ptrdiff_t i = 0; i < n_nodes_owned; ++i) {
      for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
        b_nodal_wide[i * wide_block_size + c] =
            (geom_t)(10 * b_points[0][i] + (geom_t)(c + 1));
      }
    }
    for (ptrdiff_t i = n_nodes_owned; i < n_nodes_local; ++i) {
      for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
        b_nodal_wide[i * wide_block_size + c] =
            -std::numeric_limits<geom_t>::max();
      }
    }

    full_exchange->gather(b_nodal_wide, wide_block_size);
    for (ptrdiff_t i = 0; i < n_nodes_local; ++i) {
      for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
        const geom_t expected =
            (geom_t)(10 * b_points[0][i] + (geom_t)(c + 1));
        const geom_t actual = b_nodal_wide[i * wide_block_size + c];
        if (actual != expected) {
          SMESH_ERROR(
              "full gather(wide-vector) mismatch at node %ld comp %ld: %f != %f",
              i, c, actual, expected);
          return SMESH_FAILURE;
        }
      }
    }

    for (ptrdiff_t i = 0; i < n_nodes_owned; ++i) {
      for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
        b_nodal_wide[i * wide_block_size + c] = 0;
      }
    }
    for (ptrdiff_t i = ghost_begin; i < ghost_end; ++i) {
      for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
        b_nodal_wide[i * wide_block_size + c] = (geom_t)(c + 1);
      }
    }
    for (ptrdiff_t i = aura_begin; i < n_nodes_local; ++i) {
      for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
        b_nodal_wide[i * wide_block_size + c] = (geom_t)(2 * (c + 1));
      }
    }

    full_exchange->scatter_add(b_nodal_wide, wide_block_size);

    for (ptrdiff_t c = 0; c < wide_block_size; ++c) {
      double local_owned_sum = 0.0;
      for (ptrdiff_t i = 0; i < n_nodes_owned; ++i) {
        local_owned_sum += b_nodal_wide[i * wide_block_size + c];
      }

      double local_expected =
          (double)(c + 1) * (double)n_nodes_ghosts +
          (double)(2 * (c + 1)) * (double)n_nodes_aura;
      double global_owned_sum = local_owned_sum;
      double global_expected = local_expected;
#if defined(SMESH_ENABLE_MPI)
      MPI_Allreduce(&local_owned_sum, &global_owned_sum, 1, MPI_DOUBLE, MPI_SUM,
                    comm->get());
      MPI_Allreduce(&local_expected, &global_expected, 1, MPI_DOUBLE, MPI_SUM,
                    comm->get());
#endif

      if (global_owned_sum != global_expected) {
        SMESH_ERROR(
            "ghost+aura scatter(wide-vector) mismatch at comp %ld: %f != %f",
            c, global_owned_sum, global_expected);
        return SMESH_FAILURE;
      }
    }
  }

  return SMESH_SUCCESS;
}
