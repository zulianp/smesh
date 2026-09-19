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

## What works today

| | Serial | MPI |
|--|:------:|:---:|
| Folder IO | yes | yes |
| `refine` (linear HEX/TET/TRI/QUAD/WEDGE/PYRAMID/EDGE, including shells) | yes | yes |
| `promote_to` (TET10/TET15, TRI6, QUAD9, HEX27, matching shells) | yes | yes |
| Skin / selector sidesets | yes | yes |
| Python VTK/Exodus converters | yes | n/a (serial scripts) |

There are no public C++ Python bindings yet. `refine` of TET10/HEX27/`PROTEUS_*` is out of scope.
