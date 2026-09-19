# MULTIBLOCK.md

## Goal
Enable all existing mesh features to work when `n_blocks() > 1` in both serial and MPI paths, including mixed-topology meshes, while keeping current single-block behavior and output formats unchanged.

This plan is functional-first: it is organized by capability, not by abstraction layers.

## Scope
- No breaking API changes.
- Public method signatures remain unchanged.
- Single-block behavior remains byte-compatible unless explicitly documented.
- Mixed-topology block support is enabled where technically feasible.
- Unsupported mixed-topology/MPI combinations must fail with explicit diagnostics.

---

## Execution plan (feature workstreams)

### 1) Make core transformations and extractors block-complete
Implement block-aware behavior for all operations below, remove single-block hard stops and run per block with deterministic remapping.

- `src/frontend/smesh_mesh.cpp`
  - `Mesh::convert_to`
  - `Mesh::promote_to`
  - `Mesh::refine`
  - `Mesh::split_block`
  - `Mesh::split_boundary_layer`
  - `Mesh::renumber_nodes` and `renumber_nodes(const std::vector<node_idx_t>&)`
  - `Mesh::skin`
  - `Mesh::skin_sideset`
  - `Mesh::extrude`
  - `Mesh::extract_sharp_edges`
  - `Mesh::extract_disconnected_faces`

Required behavior for each:
- remove `n_blocks() != 1` assertions/errors
- iterate blocks that participate in the request
- preserve deterministic ID order across blocks
- build outputs with global ID remapping:
  - node ids
  - edge ids
  - face ids
  - element ids
  - sideset ids / names
- for operations dependent on element types, dispatch per block using each block’s `element_type` and element type-specific topology, never rely on `element_type(0)`.

Common cross-cutting rules:
- use explicit ownership checks on shared interfaces before deduplicating faces/edges across block boundaries.
- avoid implicit reordering by block iteration order unless documented and deterministic.
- ensure all newly created global ids are monotonic across all blocks.

---

### 2) Core graph / topology APIs must be multiblock-complete
- `Mesh::half_face_table` in `src/frontend/smesh_mesh.cpp`
  - currently blocked on `n_blocks() > 1` and must be replaced by deterministic multi-block construction.
  - derive dual relationships using existing multiblock graph helpers in `src/graph/smesh_multiblock_graph.*`.
  - ensure interface face matching uses canonical `(block0, side, block1, side)` ordering to keep stable ids and face orientation.
- `Mesh::create_node_to_node_graph` in `src/frontend/smesh_mesh.cpp`
  - remove single-block guard and add cross-block node sharing handling.
  - support both scalar and per-element-type dispatch paths.
- `src/graph/smesh_graph.*`, `src/graph/smesh_adjacency.*`
  - add/adjust tests-level parity checks to ensure shared-boundary nodes/edges are represented once and consistently.
- `src/frontend/smesh_dual_graph.cpp`
  - validate all dual graph construction paths can accept multi-block dual from updated half-face table.

Required behavior:
- no duplicate adjacency entries for shared boundary objects
- shared boundary orientation must be consistent and stable
- every adjacency query result deterministic independent of block read order

---

### 3) Sidesets behave correctly over multiple blocks
- `src/frontend/smesh_sideset.cpp`
- `src/frontend/smesh_sideset.hpp`
- `src/frontend/smesh_extractions.cpp`

Targets:
- `Sideset::create_from_selector`
- selector-based extraction/conversion flows
- `Mesh::skin_sideset`
- `Mesh::create_surface_from_sidesets` and `Mesh::create_surface_from_sideset`

Required behavior:
- create/select/update sidesets across all blocks that satisfy selectors.
- deduplicate shared boundary faces/edges once, with deterministic owner assignment.
- enforce global numbering remap for all parent references (`Face`, `Node` ids).
- resolve multiple same-named sidesets from different blocks deterministically to one global sideset id per resolved canonical name.
- emit explicit error with operation/block context when deterministic merge is impossible.

---

### 4) Semistructured and extraction paths
- `src/frontend/smesh_semistructured.cpp`
  - `to_semistructured`, `derefine`
  - remove one-block guard
  - process all blocks and merge artifacts globally

- `src/frontend/smesh_extractions.cpp`
  - ensure extracted faces/edges/sets are merged across blocks
  - avoid duplicates on shared interfaces

- `src/frontend/smesh_extract_shape_features.impl.hpp`
  - validate feature extraction paths for merged results on multi-block input.

Required behavior:
- mixed-topology blocks are preserved when output format allows mixed blocks
- extracted geometry is globally unique and reuses shared entities.

---

### 5) Reorder/layout-dependent operations
- `src/frontend/smesh_mesh_reorder.cpp`
  - `SFC::reorder`

Requirements:
- remove multi-block guard
- maintain one global ordering strategy that works with:
  - multiple blocks
  - mixed topologies
  - shared-node boundaries
- determinism guarantee:
  - same mesh + same block definitions + same run => same ids
  - no order dependence from partitioning or internal map iteration

---

### 6) MPI distributed read/write with multi-block support
- `src/distributed/smesh_distributed_read.*`
  - discover and load multi-block layout from metadata/folder structure.
  - construct per-block topologies and global remap under partitioning.
  - preserve single-block distributed path behavior for compatibility.
- `src/distributed/smesh_distributed_write.*`
  - emit multi-block metadata and per-block payloads consistently.
  - keep single-block emitted layout identical to current behavior.
- `src/frontend/smesh_mesh.cpp`
  - `Mesh::read` / `Mesh::write`
  - remove all explicit multi-block hard-stop branches in MPI modes.
- `src/io/smesh_read.impl.hpp`
- `src/io/smesh_write.cpp` / `src/io/smesh_write.hpp`
  - verify format negotiation and metadata contracts across ranks.

Decision policy:
- if a particular MPI mode cannot support mixed-topology multi-block today, emit explicit error message:
  - operation name
  - block id(s)
  - reason
  - known workaround

---

### 7) Validation, diagnostics, and unsupported cases
Across all touched paths, replace generic rejections by explicit diagnostics with:
- operation name
- block index and/or block set
- concrete reason
- actionable workaround

If mixed-topology support is fundamentally unsafe for a specific operation+mode, it remains unsupported and must fail with a single clear message; do not hide behind assertions.

---

## Deterministic multi-block remapping rules
Use a consistent order for all cross-block entity synthesis:
1. block index ascending
2. local index ascending inside each block
3. stable tie-breaker on local/global parent ids
4. stable ordering of sideset names

This ordering rule must be documented in implementation comments near remap code.

---

## File-by-file worklist
- `src/frontend/smesh_mesh.cpp` (primary feature lifting)
- `src/frontend/smesh_sideset.cpp`
- `src/frontend/smesh_sideset.hpp`
- `src/frontend/smesh_extractions.cpp`
- `src/frontend/smesh_semistructured.cpp`
- `src/frontend/smesh_mesh_reorder.cpp`
- `src/frontend/smesh_dual_graph.cpp`
- `src/graph/smesh_multiblock_graph.hpp`
- `src/graph/smesh_multiblock_graph.impl.hpp`
- `src/graph/smesh_adjacency.*`
- `src/graph/smesh_graph.*`
- `src/distributed/smesh_distributed_read.*`
- `src/distributed/smesh_distributed_write.*`
- `src/io/smesh_read.impl.hpp`
- `src/io/smesh_write.*`
- `src/shape/smesh_extract_shape_features.*`

---

## Required test matrix
- Regression for every previously blocked method now succeeding with `n_blocks() > 1`.
- Serial read/write round trip for:
  - homogenous topology multi-block
  - mixed topologies
- MPI multi-block read/write round trip.
- Cross-block sideset extraction/creation correctness:
  - duplicate interface entity dedup
  - global id/name resolution
- `half_face_table` parity checks on:
  - single-block
  - two-block interface mesh
  - mixed-topology two-block mesh
- dual graph and node/edge adjacency consistency on shared boundaries.
- deterministic reorder checks with mixed block order.

---

## Risks and explicit unsupporteds to track
- Distributed refinement/promote/extrude paths with cross-rank shared-boundary synchronization.
- Fully conforming mixed-topology write/read semantics for external formats that do not preserve multi-topology natively.
- Any operation that cannot preserve uniqueness of shared entities under current container assumptions.

When blocked, issue explicit error and keep single-block path unchanged.
