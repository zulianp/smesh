// Build a 2D mesh by hand, with no input files and no built-in generators.
//
// Domain: [0, 2] x [0, 1], made of two conforming blocks that share the nodes
// on x = 1:
//
//   y=1 +-----+-----+-----+-----+
//       |     |     | \   | \   |
//       |quad |quad |  \  |  \  |     block 0 "quads": QUAD4 on [0, 1] x [0, 1]
//       +-----+-----+-----+-----+     block 1 "tris":  TRI3  on [1, 2] x [0, 1]
//       |     |     | \   | \   |
//       |     |     |  \  |  \  |
//   y=0 +-----+-----+-----+-----+
//      x=0         x=1         x=2
//
// Sidesets:
//   "left"   (x = 0) quads only, built explicitly as (parent, lfi) pairs
//   "right"  (x = 2) tris only,  built explicitly as (parent, lfi) pairs
//   "top"    (y = 1) both blocks, built with Sideset::create_from_selector
//   "bottom" (y = 0) both blocks, built with Sideset::create_from_selector
//
// Local side (lfi) convention for QUAD4 and TRI3: side s connects local node s
// to local node (s + 1) % n_corners.

#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_sideset.hpp"
#include "smesh_tracer.hpp"

#include <cmath>
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("create_mesh_from_scratch.exe");

  auto ctx = smesh::initialize(argc, argv);

  if (argc != 4) {
    fprintf(stderr, "Usage: %s <nx_per_block> <ny> <output_folder>\n", argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto comm = ctx->communicator();

    if (comm->size() != 1) {
      fprintf(stderr, "%s: this demo builds a serial mesh, run with 1 rank\n",
              argv[0]);
      return SMESH_FAILURE;
    }

    const ptrdiff_t nx = atol(argv[1]);
    const ptrdiff_t ny = atol(argv[2]);
    const Path output_folder(argv[3]);

    const geom_t xmin = 0, xmid = 1, xmax = 2;
    const geom_t ymin = 0, ymax = 1;

    // ------------------------------------------------------------------------
    // 1) Points: one shared (2 x n_nodes) SoA array for all blocks
    // ------------------------------------------------------------------------
    const ptrdiff_t nxn = 2 * nx + 1;
    const ptrdiff_t nyn = ny + 1;
    const ptrdiff_t n_nodes = nxn * nyn;

    auto points = create_host_buffer<geom_t>(2, n_nodes);
    {
      auto xy = points->data();
      const geom_t hx = (xmax - xmin) / (2 * nx);
      const geom_t hy = (ymax - ymin) / ny;
      for (ptrdiff_t j = 0; j < nyn; j++) {
        for (ptrdiff_t i = 0; i < nxn; i++) {
          const ptrdiff_t node = j * nxn + i;
          xy[0][node] = xmin + i * hx;
          xy[1][node] = ymin + j * hy;
        }
      }
    }

    auto node_id = [nxn](ptrdiff_t i, ptrdiff_t j) -> idx_t {
      return static_cast<idx_t>(j * nxn + i);
    };

    // ------------------------------------------------------------------------
    // 2) Block 0: QUAD4 cells (counter-clockwise), columns [0, nx)
    //    Connectivity is SoA: quads[d][e] is local node d of element e.
    // ------------------------------------------------------------------------
    const ptrdiff_t n_quads = nx * ny;
    auto quads = create_host_buffer<idx_t>(4, n_quads);
    {
      auto q = quads->data();
      for (ptrdiff_t j = 0; j < ny; j++) {
        for (ptrdiff_t i = 0; i < nx; i++) {
          const ptrdiff_t e = j * nx + i;
          q[0][e] = node_id(i, j);
          q[1][e] = node_id(i + 1, j);
          q[2][e] = node_id(i + 1, j + 1);
          q[3][e] = node_id(i, j + 1);
        }
      }
    }

    // ------------------------------------------------------------------------
    // 3) Block 1: TRI3 cells, columns [nx, 2 nx), each cell split along the
    //    diagonal (i, j) -> (i + 1, j + 1) into two counter-clockwise triangles.
    // ------------------------------------------------------------------------
    const ptrdiff_t n_tris = 2 * nx * ny;
    auto tris = create_host_buffer<idx_t>(3, n_tris);
    {
      auto t = tris->data();
      for (ptrdiff_t j = 0; j < ny; j++) {
        for (ptrdiff_t ii = 0; ii < nx; ii++) {
          const ptrdiff_t i = nx + ii;
          const ptrdiff_t e = 2 * (j * nx + ii);

          // Lower-right triangle: side 0 bottom, side 1 right, side 2 diagonal
          t[0][e] = node_id(i, j);
          t[1][e] = node_id(i + 1, j);
          t[2][e] = node_id(i + 1, j + 1);

          // Upper-left triangle: side 0 diagonal, side 1 top, side 2 left
          t[0][e + 1] = node_id(i, j);
          t[1][e + 1] = node_id(i + 1, j + 1);
          t[2][e + 1] = node_id(i, j + 1);
        }
      }
    }

    // ------------------------------------------------------------------------
    // 4) Assemble the mesh from its blocks and points
    // ------------------------------------------------------------------------
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(std::make_shared<Mesh::Block>("quads", QUAD4, quads));
    blocks.back()->set_geom_map(AXIS_ALIGNED);
    blocks.push_back(std::make_shared<Mesh::Block>("tris", TRI3, tris));
    blocks.back()->set_geom_map(AFFINE);

    auto mesh = std::make_shared<Mesh>(comm, blocks, points);

    // ------------------------------------------------------------------------
    // 5a) Explicit sidesets: a side is (parent element, local face index) in
    //     a given block. Parent indices are local to that block.
    // ------------------------------------------------------------------------
    {
      const block_idx_t quad_block = 0;
      auto parent = create_host_buffer<element_idx_t>(ny);
      auto lfi = create_host_buffer<i16>(ny);
      for (ptrdiff_t j = 0; j < ny; j++) {
        parent->data()[j] = static_cast<element_idx_t>(j * nx + 0);
        lfi->data()[j] = 3; // QUAD4 side 3: node 3 -> node 0 (x = xmin)
      }
      mesh->add_sideset("left",
                        Sideset::create(comm, parent, lfi, quad_block));
    }

    {
      const block_idx_t tri_block = 1;
      auto parent = create_host_buffer<element_idx_t>(ny);
      auto lfi = create_host_buffer<i16>(ny);
      for (ptrdiff_t j = 0; j < ny; j++) {
        // Lower-right triangle of the last column
        parent->data()[j] = static_cast<element_idx_t>(2 * (j * nx + nx - 1));
        lfi->data()[j] = 1; // TRI3 side 1: node 1 -> node 2 (x = xmax)
      }
      mesh->add_sideset("right", Sideset::create(comm, parent, lfi, tri_block));
    }

    // ------------------------------------------------------------------------
    // 5b) Geometric sidesets: selector is evaluated at each side's barycenter.
    //     Returns one Sideset per block that has matching sides.
    // ------------------------------------------------------------------------
    const geom_t tol = 1e-8;
    mesh->add_sidesets(
        "top", Sideset::create_from_selector(
                   mesh, [=](const geom_t, const geom_t y, const geom_t) {
                     return std::abs(y - ymax) < tol;
                   }));
    mesh->add_sidesets(
        "bottom", Sideset::create_from_selector(
                      mesh, [=](const geom_t, const geom_t y, const geom_t) {
                        return std::abs(y - ymin) < tol;
                      }));

    // ------------------------------------------------------------------------
    // 6) Report and write
    // ------------------------------------------------------------------------
    printf("#nodes %ld, spatial dimension %d, x-interface at %g\n",
           (long)mesh->n_nodes(), mesh->spatial_dimension(), (double)xmid);
    for (size_t b = 0; b < mesh->n_blocks(); b++) {
      auto block = mesh->block(b);
      printf("block %zu \"%s\": %s x %ld\n", b, block->name().c_str(),
             type_to_string(block->element_type()), (long)block->n_elements());
    }
    for (auto &named : mesh->sidesets()) {
      printf("sideset \"%s\": block %d, %ld sides\n", named.first.c_str(),
             (int)named.second->block_id(), (long)named.second->size());
    }

    if (mesh->write(output_folder) != SMESH_SUCCESS) {
      SMESH_ERROR("Failed to write mesh to %s", output_folder.c_str());
      return SMESH_FAILURE;
    }

    // Round trip: blocks and sidesets are restored from disk
    auto reloaded = Mesh::create_from_file(comm, output_folder);
    printf("reloaded: %zu blocks, %zu sidesets\n", reloaded->n_blocks(),
           reloaded->sidesets().size());
  }

  return SMESH_SUCCESS;
}
