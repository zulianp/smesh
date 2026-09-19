# Concepts

## Mesh

A `Mesh` is one shared node pool plus one or more named **blocks**. Each block has a single element type and a connectivity SoA. Serial and MPI use the same types; distributed meshes add ownership, ghosts, and aura.

## Sets

Sidesets, edgesets, and nodesets are topology, not geometry.

- **Sideset:** `(block_id, parent element, local face index)`. `lfi` is 0-based (Exodus side minus one).
- **Edgeset:** `(block_id, parent, local edge index)`.
- **Nodeset:** node ids (local; MPI write stores serial GIDs).

`Mesh::read` / `Mesh::write` round-trip named sets when the registry is non-empty (`sidesets/<name>/`, `nodesets/<name>/`).

## Refine, promote, semistructured

| Operation | What changes | Typical use |
|-----------|----------------|-------------|
| `refine(mesh, levels)` | More elements, same linear type | h-refinement |
| `promote_to(TET10, mesh)` | Same elements, extra nodes | p-refinement |
| `to_semistructured` / `derefine` | Lattice nodes on the same macros | `PROTEUS_*` types |

`refine` refuses higher-order input. Converters (`db_to_raw` / `raw_to_db`) only translate file formats.

Element-type coverage is in [features.md](features.md). There are no public C++ Python bindings yet.
