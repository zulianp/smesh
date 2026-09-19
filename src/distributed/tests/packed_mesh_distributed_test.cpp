// Packing a distributed mesh must leave the ownership metadata describing the mesh.
//
// PackedMesh::create(..., modify_mesh = true) renumbers the nodes in place, and until this
// test existed nothing checked what that did to Distributed. It did the wrong thing, and the
// wrongness was not visible from the rank doing it:
//
//   node_mapping (local -> global id) and node_owner are indexed by local node, so a
//   renumbering that did not permute them left every local slot claiming whatever belonged
//   to the node that used to sit there.
//
//   ghosts_and_aura is worse, because it is not this rank's to repair. Its entries name
//   nodes by GLOBAL owned index -- exchange_create recovers a local index by subtracting
//   node_offsets[rank] -- so they record the position a node holds inside its OWNER's owned
//   block, and the ranks holding those entries are the neighbours. Renumbering a node that a
//   neighbour references silently redirects that neighbour's gather.
//
// The bound that makes packing safe is owned-not-shared. A node is shared iff one of its
// incident elements belongs to another rank, so an owned-not-shared node has every incident
// element here and can appear in no neighbour's ghost list. That prefix may be permuted
// freely; nothing else may move.
//
// ORDER MATTERS HERE. The checks that need no node_map come first, deliberately, so that
// this test still reaches the cross-rank assertion when run against a build whose packing
// releases the map. Otherwise it would stop at a null pointer and report nothing about the
// corruption it exists to catch.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <vector>

#include "smesh_communicator.hpp"
#include "smesh_exchange.hpp"
#include "smesh_mesh.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_test.hpp"

using namespace smesh;

namespace {

using pack_idx_t = uint16_t;

// A value derived from the node's COORDINATES, and the choice is load-bearing rather than
// cosmetic.
//
// The obvious witness is a function of the global id read from node_mapping, and it was
// tried first. It does not work, and measuring why is what produced this comment: packing
// without repairing the metadata leaves node_mapping, node_owner and ghosts_and_aura ALL
// stale together, so filling owned slots from node_mapping and checking ghosts against
// node_mapping uses one consistent stale convention on both sides. The gather round-trips
// perfectly and the corruption goes unseen.
//
// The damage is not that the metadata disagrees with itself; it is that the metadata
// describes where the nodes USED to be while points and connectivity describe where they
// are now. Only a quantity tied to the mesh itself can see that -- so the witness is built
// from the node's position, which no renumbering can alter.
double witness(const geom_t *const *const pts, const ptrdiff_t node, const int dim) {
  double v = 0.0;
  for (int d = 0; d < dim; ++d) {
    v = 1.7 * v + static_cast<double>(pts[d][node]);
  }
  return v;
}

}  // namespace

int test_packed_mesh_preserves_distributed_metadata() {
#ifndef SMESH_ENABLE_MPI
  return SMESH_TEST_SUCCESS;
#else
  auto comm = Communicator::world();
  if (comm->size() == 1) {
    // Nothing to corrupt: one rank has no shared nodes, no ghosts, and no neighbour
    // holding indices into it.
    return SMESH_TEST_SUCCESS;
  }

  char path_buffer[256];
  std::snprintf(path_buffer, sizeof(path_buffer),
                "/tmp/smesh_packed_mesh_distributed_test_%d", comm->size());
  const Path mesh_path(path_buffer);

  const ptrdiff_t nx = std::max<ptrdiff_t>(2, comm->size());

  if (comm->rank() == 0) {
    std::filesystem::remove_all(mesh_path.to_string());
    auto serial_mesh = Mesh::create_tet4_cube(Communicator::self(), nx, 6, 6);
    SMESH_TEST_ASSERT(serial_mesh != nullptr);
    SMESH_TEST_ASSERT(serial_mesh->write(mesh_path) == SMESH_SUCCESS);
  }

  comm->barrier();

  auto mesh = std::make_shared<Mesh>(comm);
  SMESH_TEST_ASSERT(mesh->read(mesh_path) == SMESH_SUCCESS);

  auto dist = mesh->distributed();
  SMESH_TEST_ASSERT(dist != nullptr);

  const ptrdiff_t n_nodes = mesh->n_nodes();
  const ptrdiff_t n_owned = dist->n_nodes_owned();
  const ptrdiff_t n_ghosts = dist->n_nodes_ghosts();
  const ptrdiff_t n_fixed = dist->n_nodes_owned_not_shared();
  const int dim = mesh->spatial_dimension();

  // Reported rather than asserted. Packing restricts itself to the owned-not-shared
  // ELEMENT segment, and if a block carries no distributed element counts that restriction
  // silently falls back to the whole block. Printing the numbers is what distinguishes
  // "the restriction held" from "the restriction was never engaged", which two runs with
  // the same verdict otherwise cannot tell apart.
  {
    auto block = mesh->block(0);
    std::printf(
        "[rank %d] nodes: local %td owned %td owned-not-shared %td ghosts %td | "
        "elements: local %td owned %td owned-not-shared %td ghosts %td\n",
        comm->rank(), n_nodes, n_owned, n_fixed, n_ghosts, block->n_elements(),
        block->n_elements_owned(), block->n_elements_owned_not_shared(),
        block->n_elements_ghosts());
    std::fflush(stdout);
  }

  // Snapshot before packing: the global id and the coordinates sitting in each local slot.
  std::vector<large_idx_t> g_before(static_cast<size_t>(n_nodes));
  {
    auto d_map = dist->node_mapping()->data();
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
      g_before[static_cast<size_t>(i)] = d_map[i];
    }
  }
  std::vector<geom_t> x_before(static_cast<size_t>(n_nodes * dim));
  {
    auto pts = mesh->points()->data();
    for (int d = 0; d < dim; ++d) {
      for (ptrdiff_t i = 0; i < n_nodes; ++i) {
        x_before[static_cast<size_t>(i * dim + d)] = pts[d][i];
      }
    }
  }

  // A small pack size so the mesh really is split over many packs; with a single pack the
  // permutation is near the identity and proves much less.
  auto packed = PackedMesh<pack_idx_t>::create(mesh, {}, true, 8);
  SMESH_TEST_ASSERT(packed != nullptr);

  // A. Shared, ghost and aura nodes did not move.
  //
  // Checked against the coordinate snapshot, not against node_map, so this runs even where
  // the map has been released. Exact comparison is right: a permutation copies coordinates,
  // it does not compute with them.
  for (ptrdiff_t i = n_fixed; i < n_nodes; ++i) {
    auto pts = mesh->points()->data();
    for (int d = 0; d < dim; ++d) {
      SMESH_TEST_ASSERT(pts[d][i] == x_before[static_cast<size_t>(i * dim + d)]);
    }
  }

  // B. The cross-rank assertion: a gather must still deliver each ghost the value belonging
  //    to the node that is actually sitting there.
  //
  // Each owned slot is filled from its own coordinates. The owning rank produces the ghost
  // values out of that array, indexed through ghosts_and_aura, so the value arriving in a
  // ghost slot is whatever the owner holds at the position its exchange indices name. The
  // receiver then compares it against the coordinates of the node in that slot. The two
  // agree only if both ranks still agree which physical node lives where.
  //
  // Both sides read the same float bits for the same node, and both apply the same
  // arithmetic, so agreement is exact; the tolerance guards the comparison, it does not
  // absorb a discrepancy.
  {
    auto pts = mesh->points()->data();
    std::vector<double> f(static_cast<size_t>(n_nodes), 0.0);
    for (ptrdiff_t i = 0; i < n_owned; ++i) {
      f[static_cast<size_t>(i)] = witness(pts, i, dim);
    }

    auto exchange =
        Exchange::create_nodal(mesh, Exchange::ExchangeScope::GhostsOnly);
    SMESH_TEST_ASSERT(exchange != nullptr);
    SMESH_TEST_ASSERT(exchange->gather<double>(f.data()) == SMESH_SUCCESS);

    for (ptrdiff_t i = n_owned; i < n_owned + n_ghosts; ++i) {
      const double expected = witness(pts, i, dim);
      SMESH_TEST_ASSERT(std::abs(f[static_cast<size_t>(i)] - expected) < 1e-12);
    }
  }

  // C. node_map outlives the renumbering, is a permutation, and describes where each node
  //    went. Releasing it is what left map_to_packed and map_to_unpacked with no working
  //    path in any configuration: both refuse unless synched_with_mesh, and that flag was
  //    set only on the branch that had just freed the map they dereference.
  auto node_map = packed->node_map();
  SMESH_TEST_ASSERT(node_map != nullptr);
  auto d_node_map = node_map->data();

  {
    std::vector<int> seen(static_cast<size_t>(n_nodes), 0);
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
      const ptrdiff_t to = static_cast<ptrdiff_t>(d_node_map[i]);
      SMESH_TEST_ASSERT(to >= 0);
      SMESH_TEST_ASSERT(to < n_nodes);
      SMESH_TEST_ASSERT(seen[static_cast<size_t>(to)] == 0);
      seen[static_cast<size_t>(to)] = 1;
    }
  }

  for (ptrdiff_t i = n_fixed; i < n_nodes; ++i) {
    SMESH_TEST_ASSERT(d_node_map[i] == static_cast<idx_t>(i));
  }

  {
    auto d_map_after = dist->node_mapping()->data();
    auto pts = mesh->points()->data();
    for (ptrdiff_t i = 0; i < n_nodes; ++i) {
      const ptrdiff_t to = static_cast<ptrdiff_t>(d_node_map[i]);
      SMESH_TEST_ASSERT(d_map_after[to] == g_before[static_cast<size_t>(i)]);
      for (int d = 0; d < dim; ++d) {
        SMESH_TEST_ASSERT(pts[d][to] == x_before[static_cast<size_t>(i * dim + d)]);
      }
    }
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
  SMESH_RUN_TEST(test_packed_mesh_preserves_distributed_metadata);
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
