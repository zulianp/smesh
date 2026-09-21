# Features and element types

Serial and MPI unless noted. Shells (`TRISHELL*`, `QUADSHELL*`, `EDGESHELL*`) follow the matching TRI / QUAD / EDGE row unless a cell says otherwise.

Marks: **yes** supported · **part** partial · **—** not supported or not applicable.

## Linear families

| | HEX8 | TET4 | WEDGE6 | PYRAMID5 | QUAD4 | TRI3 | EDGE2 |
|--|:----:|:----:|:------:|:--------:|:-----:|:----:|:-----:|
| Folder IO | yes | yes | yes | yes | yes | yes | yes |
| `write_with_xdmf` | yes | yes | yes | yes | yes | yes | yes |
| `create_cube` / `create_square` | cube | cube | — | — | square | square | — |
| `refine` | yes | yes | yes | yes | yes | yes | yes |
| `adapt_refine` | — | yes | — | — | yes | yes | — |
| `improve` / `remesh` | — | yes | — | — | yes | yes | — |
| `to_semistructured` / `derefine` | yes | yes | yes | yes | yes | — | — |
| Sidesets (skin / selector) | yes | yes | yes | yes | yes | yes | yes |
| Sidesets through `refine` | yes | yes | yes | part | yes | yes | yes |
| Edgesets through `refine` | yes | yes | yes | — | yes | yes | yes |
| Nodesets through `refine` / `promote` | yes | yes | yes | yes | yes | yes | yes |
| `extrude` | — | — | — | — | → HEX8 | → WEDGE6 | — |
| Jacobians / FFF | yes | yes | — | — | planar | planar | — |
| SS restrict (host) | yes | yes | — | — | yes | — | — |
| SS restrict (CUDA) | yes | — | — | — | — | — | — |
| SS prolong | yes | yes | — | — | yes | — | — |
| SFC reorder | yes | yes | yes | yes | yes | yes | yes |
| Graphs (n2n / n2e / dual) | yes | yes | yes | yes | yes | yes | yes |
| `split_block` | yes | yes | yes | yes | yes | yes | yes |

`refine` is linear (P1) only. PYRAMID `refine` emits a PYRAMID5 block plus a TET4 block; sidesets through refine cover the quad base (`lfi == 4`), not the triangular sides. Nodesets through refine/promote add mid-edge nodes when both coarse endpoints are already in the set (not face or body mids). Planar QUAD4 / TRI3 use 2-component points; shells stay in 3D.

`adapt_refine` is a separate serial OpenMP API (not uniform `refine()`). Single-block `TET4` and `TRI3`/`TRISHELL3`/`QUAD4`/`QUADSHELL4` only. Size field from Meyer curvature (sharp creases skipped), 2:1 gradation, conforming newest-vertex / longest-edge bisection (TRI/TET) or quad 4-split with 2:1 closure. Optional Jacobi smooth with corner pin / crease slide, then parametrization `apply`. HEX/WEDGE/PYRAMID/mixed/MPI are rejected.

`improve` / `remesh` is a separate serial OpenMP quality remesher (not `refine()` / `adapt_refine()`). Same types as `adapt_refine`. Mean-ratio split / collapse / swap (QUAD: 4-split + 2:1 only), feature locks from sharp edges (corners pinned, creases slide), surface nodes clipped to a displacement band about the input (`0.02 * bbox_diag` default). `improve` mutates in place; `remesh` clones then `improve`. HEX/WEDGE/PYRAMID/mixed/MPI are rejected.

Mixed-volume `refine` and SS (HEX / TET / WEDGE / PYRAMID, including hex-dominant) are supported. Mixed HEX+QUAD is not. Mixed QUAD4+QUADSHELL4 is supported. GLL nodes on SS are HEX-only.

## Promote (`promote_to`)

Same element count, extra nodes. Same-type blocks only. Sidesets keep `(parent, lfi)`.

| From | To |
|------|-----|
| HEX8 | HEX27 |
| TET4 | TET10, TET15 |
| QUAD4 | QUAD9 |
| QUADSHELL4 | QUADSHELL9 |
| TRI3 | TRI6 |
| TRISHELL3 | TRISHELL6 |
| EDGE2 | — |

## Convert (`convert_to`)

| From | To |
|------|-----|
| HEX8 | TET4 (6 tets) |
| WEDGE6 | TET4 (3 tets) |
| PYRAMID5 | TET4 (2 tets) |
| QUAD4 | TRI3 (2 tris) |
| TET15 | HEX8 (4 hexes) |
| QUADSHELL4 | — (no TRISHELL3) |
| PROTEUS_HEX* | HEX8 (not `PROTEUS_HEX4913`) |
| PROTEUS_TET* | TET4 |
| PROTEUS_QUAD* / QUADSHELL* | linear QUAD (not `*289`) |
| PROTEUS_WEDGE* | WEDGE6 |
| PROTEUS_PYRAMID* | — (use `ss_to_linear` / `derefine`) |

## Higher-order and lattice types

| | TET10 / TET15 | TRI6 | QUAD9 | HEX27 | PROTEUS_* |
|--|:-------------:|:----:|:-----:|:-----:|:---------:|
| Folder IO | yes | yes | yes | yes | yes |
| `refine` | — | — | — | — | — |
| `promote_to` source | — | — | — | — | — |
| Jacobians / FFF | TET10 | — | — | — | HEX macros |
| Unstructured restrict | TET10 / MACRO_TET4 → TET4 | TRI6 / MACRO_TRI3 → TRI3 | — | — | — |
| SS restrict / prolong | — | — | — | — | HEX / TET / QUAD (see linear table) |

Device SS restrict (`Restrict` with `EXECUTION_SPACE_DEVICE`) is HEX-family only. Unstructured multi-block restrict is not implemented.

## CUDA

Off by default (`-DSMESH_ENABLE_CUDA=ON`). Device `Buffer`s (`to_device` / `to_host`), `Mesh::device_points_*` / `device_elements_*`, and kernels for:

| Kernel | Coverage |
|--------|----------|
| SSHEX restrict / prolong | HEX SS |
| SSQUAD restrict / prolong | QUAD SS (kernel; mesh `Restrict` DEVICE path is HEX-only) |
| TET4 ↔ MACRO_TET4 prolong / restrict | unstructured TET |
| Unstructured TET10 / TRI6 restrict / prolong | see higher-order table |

WEDGE / PYRAMID have no CUDA restrict or prolong.

## Factories

| Factory | Types |
|---------|--------|
| `create_cube` | HEX8, TET4, HEX27, TET10, `PROTEUS_HEX*` |
| `create_square` | QUAD4, TRI3, TRI6, `PROTEUS_QUAD*` |
| `create_semistructured_hex_cube` | `PROTEUS_HEX*` |
| `create_semistructured_quad_square` | `PROTEUS_QUAD*` |
| `create_quad4_ring` | QUAD4 |
| `create_half_sphere` | HEX8, TET4 |
| `create_hex8_nozzle` / `create_hex8_lshape` | HEX8 |
| `create_hex8_checkerboard_cube` / `create_hex8_bidomain_cube` | HEX8 (two blocks) |
| `create_hex8_tet4_cube` | HEX8 + TET4 |
| `create_hex_dominant_serial` / `create_hex_dominant_cylinder` | HEX8 + WEDGE6 + PYRAMID5 + TET4 |
| `create_wall_mounted_hump` | HEX8, TET4 |

No WEDGE6/PYRAMID5 cube, no QUADSHELL/TRISHELL square, no mixed QUAD+TRI factory.

## Also

- **GeomMap** (per block): `Affine`, `AxisAligned` (HEX/QUAD families, including shells and Proteus), `IsoParametric` (default). `meta.yaml` `geom_map:` round-trips; `detect_geom_map` classifies from coordinates.
- **Parametrization:** circle (including planar `sdim == 2`), sphere, polynomial surface; nodeset-based. Not written by `Mesh::read` / `write`.
- **SFC:** serial and MPI, including `sdim == 2` and multiblock. SS hierarchical numbering stays single-block.
- **Other mesh ops:** `clone`, `concatenate`, `skin` / `mesh_from_sideset`, `split_boundary_layer`, `renumber_nodes`.
- **Python converters:** VTK / Exodus ↔ folder mesh; translation only.
- No public C++ Python bindings yet.
