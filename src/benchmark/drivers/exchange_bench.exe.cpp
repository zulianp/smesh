// Exchange::gather of a nodal vector on a mesh read from disk, as a distributed operator
// synchronises its input before every application (ghosts and aura, block_size components per
// node). Times the mean gather over `repeats`, once with the ranks arriving together (a barrier
// before each, so what remains is the exchange itself) and once back to back; checks the
// gathered values; reports the sizes. Usage: exchange_bench <mesh_folder> [repeats] [block_size]

#include "smesh_context.hpp"
#include "smesh_exchange.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <set>
#include <vector>

using namespace smesh;

int main(int argc, char **argv) {
  auto ctx = smesh::initialize(argc, argv);
  int ret = SMESH_SUCCESS;
  {
    auto comm = ctx->communicator();
    if (argc < 2) {
      if (!comm->rank())
        fprintf(stderr, "usage: %s <mesh_folder> [repeats] [block_size]\n",
                argv[0]);
      return SMESH_FAILURE;
    }
    const int repeats = argc > 2 ? atoi(argv[2]) : 2000;
    const ptrdiff_t bs = argc > 3 ? atol(argv[3]) : 3;

    auto mesh = Mesh::create_from_file(comm, Path(argv[1]));
    auto exchange = Exchange::create_nodal(
        mesh, Exchange::ExchangeScope::GhostsAndAura);
    auto dist = mesh->distributed();
    const ptrdiff_t n_owned = dist->n_nodes_owned();
    const ptrdiff_t n_local =
        n_owned + dist->n_nodes_ghosts() + dist->n_nodes_aura();

    // Each node's x coordinate in every component: what a gather must deliver.
    const geom_t *const x = mesh->points()->data()[0];
    std::vector<double> v(n_local * bs);
    auto reset = [&]() {
      for (ptrdiff_t i = 0; i < n_local; i++)
        for (ptrdiff_t c = 0; c < bs; c++)
          v[i * bs + c] = i < n_owned ? double(x[i]) + c : -1;
    };
    reset();
    exchange->gather(v.data(), bs);
    long wrong = 0;
    for (ptrdiff_t i = 0; i < n_local; i++)
      for (ptrdiff_t c = 0; c < bs; c++)
        wrong += v[i * bs + c] != double(x[i]) + c;

    for (int k = 0; k < 50; k++) exchange->gather(v.data(), bs);

    double synced = 0;
    for (int k = 0; k < repeats; k++) {
      MPI_Barrier(comm->get());
      const double t0 = MPI_Wtime();
      exchange->gather(v.data(), bs);
      synced += MPI_Wtime() - t0;
    }
    MPI_Barrier(comm->get());
    const double t0 = MPI_Wtime();
    for (int k = 0; k < repeats; k++) exchange->gather(v.data(), bs);
    const double back_to_back = MPI_Wtime() - t0;

    // The ranks this one receives ghost and aura values from.
    const int *const owner = dist->node_owner()->data();
    std::set<int> sources;
    for (ptrdiff_t i = n_owned; i < n_local; i++) sources.insert(owner[i]);

    double mine[5] = {1e6 * synced / repeats, 1e6 * back_to_back / repeats,
                      double(n_local - n_owned), double(sources.size()),
                      double(n_owned)};
    double most[5], sum[5];
    MPI_Reduce(mine, most, 5, MPI_DOUBLE, MPI_MAX, 0, comm->get());
    MPI_Reduce(mine, sum, 5, MPI_DOUBLE, MPI_SUM, 0, comm->get());
    long all_wrong = 0;
    MPI_Reduce(&wrong, &all_wrong, 1, MPI_LONG, MPI_SUM, 0, comm->get());
    if (!comm->rank()) {
      const int p = comm->size();
      printf("exchange_bench: %d ranks, %ld nodes, block %ld, %d repeats\n", p,
             (long)sum[4], (long)bs, repeats);
      printf("  imported nodes per rank: mean %.0f, max %.0f; source ranks: "
             "mean %.1f, max %.0f\n",
             sum[2] / p, most[2], sum[3] / p, most[3]);
      printf("  gather, ranks together: mean %.1f us, slowest rank %.1f us\n",
             sum[0] / p, most[0]);
      printf("  gather, back to back:   mean %.1f us, slowest rank %.1f us\n",
             sum[1] / p, most[1]);
      printf("  wrong values: %ld\n", all_wrong);
      if (all_wrong) ret = SMESH_FAILURE;
    }
  }
  return ret;
}
