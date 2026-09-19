# SIDESET.md

## Goal

Make sidesets first-class, mesh-lifecycle-safe objects for all supported mesh types under serial and MPI execution, including refinement/hierarchies, reordering, and IO.

This plan is functional-first. Multiblock sideset merge/naming overlaps [MULTIBLOCK.md](MULTIBLOCK.md) §3; this document owns sideset lifecycle correctness beyond multiblock alone.

## Scope

- No breaking public API signatures unless explicitly versioned.
- Single-block serial behavior remains the compatibility baseline.
- Unsupported combinations must fail with explicit diagnostics (operation, block, reason, workaround).
- Keep representation: `(block_id, parent[], lfi[])` plus optional `element_mapping` for MPI serial GIDs.

---

## Current state (summary)

### Representation

`Sideset` (`src/frontend/smesh_sideset.hpp/.cpp`) stores:

| Field | Meaning |
|-------|---------|
| `parent` | Local parent element indices in `block_id` |
| `lfi` | Local face/side index (0-based in memory) |
| `block_id` | Owning block |
| `element_mapping` | Optional serial/global element GIDs (MPI write / redistribute) |

Connectivity is derived on demand via `LocalSideTable` + element SoA. **Mesh optionally owns named sidesets** (empty registry by default). `Mesh::read` / `Mesh::write` round-trip `sidesets/<name>/` when the registry is non-empty.

### What already works

| Capability | Status |
|------------|--------|
| Create from selector / plane | HEX/TET/WEDGE/PYRAMID/QUAD/TRI + SS families (corners + family LST) |
| Skin (`skin_sideset` / `skin_sidesets`) | Single-block + per-block; MPI ownership via aura face masks |
| Surface / nodeset extract | Unstructured + SSHEX/SSTET/SSQUAD/SSWEDGE/SSPYRAMID |
| MPI create ownership filter | `filter_sidesets_to_owned` |
| MPI redistribute | `read_and_redistibute` / `redistribute` via `element_mapping` |
| Parallel write | Parents remapped to serial GIDs |
| SS level change | Same `(parent,lfi)` remains valid across SS refine/derefine levels |
| Nodeset MPI sync | `synchronize_nodeset_gids` |

### Core files

- `src/frontend/smesh_sideset.*`, `smesh_device_sideset.*`
- `src/mesh/sets/smesh_sidesets.*`
- `src/graph/smesh_adjacency.*`, `smesh_volume_to_surface.*`
- `src/frontend/smesh_mesh.cpp` (`skin_sideset*`, `mesh_from_sideset*`)
- SS extractors: `sshex8` / `sstet4` / `ssquad4` / `sswedge` / `sspyramid` graph impls
- Drivers: `create_sideset`, `skin`, `surface_from_sideset`
- Python Exodus: `python/smesh/exodusII_to_raw.py`, `raw_to_exodusII.py`

---

## Gap matrix

### Mesh types

| Family | Create | Skin | Surface/nodeset | SS extract | Gaps |
|--------|--------|------|-----------------|------------|------|
| HEX8 / HEX27 | Yes | Yes | Yes | SSHEX | — |
| TET4 / TET10 | Yes | Yes | Yes | SSTET | — |
| WEDGE6 | Yes | Yes (TRI/QUAD split) | `create_surfaces_from_sidesets` auto-splits | SSWEDGE | — |
| PYRAMID5 | Yes | Yes (TRI/QUAD split) | Same | SSPYRAMID | — |
| QUAD / TRI | Yes | Yes | Yes | SSQUAD | Thin MPI parity tests |
| Mixed / hex-dominant | Per-block | `skin_sidesets` | Per-block family extract; merge by shell type | Per-block SS extract | — |
| EDGE2 / EDGE3 / EDGESHELL* | Yes | Yes | NODE1 surface / nodeset | N/A | — |
| QUADSHELL* / TRISHELL3/6 | Yes (alias of QUAD/TRI) | Yes | Yes | SSQUAD | — |
| Higher-order PROTEUS_* | Create via family corners | Skin via family | SS extractors | SS* | `LocalSideTable` stays on the linear family; higher PROTEUS uses SS paths |
| TET20 / TRI10 / NODE1 / BEAM2 | No | No | No | N/A | `fill` returns FAILURE; create/skin/extract print one diagnostic |

Hard stops today:

- `skin_sideset`: more than one non-empty per-block skin still errors (use `skin_sidesets`)
- Skin face mask assumes `elem_num_sides <= 8`

### Parallel / MPI

| Item | Status | Gap |
|------|--------|-----|
| Owned-only create/skin | Done | Aura sides intentionally dropped |
| Redistribute from serial GIDs | Done | No live remap after repartition without full redistribute |
| `Sideset::read` with `comm_size>1` | Rejected | Only `read_and_redistibute` |
| `meta.yaml` `size:` under MPI | Done | Global `MPI_SUM` size + `block_id` |
| `to_device(Sideset)` | Done | Keeps `block_id`; host `element_mapping` |
| Cross-block interface dedup | Done | Unique face by sorted corner nodes; owner `(block, parent, lfi)` |
| Named sideset registry on Mesh | Done | `add_sideset(s)` / `sidesets()`; empty by default |

### Refinement / hierarchies

| Path | Sideset behavior | Gap |
|------|------------------|-----|
| SS `to_semistructured` / `derefine` | `(parent,lfi)` stable | Document + test as contract |
| Unstructured `refine` / `promote` | Mesh registry remapped via `map_sideset_through_refine` | HEX/TET/TRI/QUAD/WEDGE/EDGE; promote is identity; unsupported pair fails `refine()` if the registry is non-empty |
| Field prolong/restrict on SS surface | Works in hierarchies test | Not a topology mapper |
| Multiblock refine + sidesets | Not integrated | MULTIBLOCK.md workstream 1 |

### Reordering

| Path | Sideset behavior | Gap |
|------|------------------|-----|
| `SFC::reorder` | Optional `vector<Sideset>&` remap **and** Mesh registry | Skip pointer-identity duplicates |
| `reorder_elements_from_tags` | Optional remap **and** Mesh registry | Same |
| `distributed_reorder_elements` | **No live update** | Mesh IO loads sidesets *after* partition (redistribute). Live in-place distributed reorder still needs a redistribute pass |
| Tests that reorder then skin | Create skin **after** reorder | Live remap covered by serial SFC/tag tests (argument + registry) |

### IO

| Format | Read | Write | Parallel | Gap |
|--------|------|-------|----------|-----|
| Custom folder (`parent.*`, `lfi.*`, `meta.yaml`) | Serial / redistribute | Yes | Write yes; global size + `block_id` | — |
| Mesh folder | Yes | Yes | MPI read uses `read_and_redistibute`; write uses `Sideset::write` | 1 sideset → `sidesets/<name>/`; N>1 → `sidesets/<name>/<block_id>/` |
| Exodus II | Python only | Python only | — | Optional native C++ later |
| VTK | No | No | — | Optional |

---

## Work plan

### Phase 0 — Correctness fixes (done)

#### 0.1 MPI `Sideset::write` metadata

- `meta.yaml` writes **global** side count (`MPI_SUM` of local sizes) and `block_id`.
- `Sideset::read` prefers meta `block_id` when present (argument is fallback for old folders) and checks meta `size` against the loaded array.

#### 0.2 `to_device(Sideset)`

- Preserves `block_id`. `element_mapping` stays a host buffer (MPI IO / redistribute only).

#### 0.3 SS sideset invariance

- `(parent, lfi)` stay valid across SS level changes (macro element + local face index unchanged). Unstructured refine/promote does **not** remap sidesets.

---

### Phase 1 — Lifecycle: reorder + refine (done)

#### 1.1 Remap on element reorder

- `Sideset::remap_parents` / `apply_element_permutation` / `remap_sidesets`.
- `SFC::reorder(Mesh&, const vector<Sideset>& = {})` and `Mesh::reorder_elements_from_tags(..., const vector<Sideset>& = {})`. Empty vector skips remap.
- Mesh still does not own sidesets; callers pass the vector. Distributed in-place SFC is unchanged (no Mesh-owned registry yet).

#### 1.2 Unstructured refine map

- `map_sideset_through_refine(coarse_mesh, coarse_ss, fine_mesh)` for HEX8 (face → L² children), QUAD4/QUADSHELL4 (edge → L children), WEDGE6 (TRI or QUAD face → L² children), TET4 (face → 4 children per level), TRI3/TRISHELL3 (edge → 2 children per level), EDGE2/EDGESHELL2 (NODE1 endpoint → 1 child).
- Factor 1 with the same side count is identity (promote).
- `refine()` remaps the Mesh sideset registry with this mapper (empty stays empty; non-empty + unsupported pair fails).

#### 1.3 Unsupported diagnostics

- SS meshes, non-cube HEX factors, non-8^k TET / non-4^k TRI, and other type pairs print a diagnostic and return `nullptr` instead of leaving stale parents.

---

### Phase 2 — Multiblock + surface ergonomics (done)

Align with [MULTIBLOCK.md](MULTIBLOCK.md) §3.

#### 2.1 Multiblock create / skin / name resolution

- Create/select already walks all matching blocks.
- Shared faces (same sorted corner node ids) are kept once; owner is `(block_id, parent, lfi)` lexicographic.
- Names remain folder/registry (Phase 3). Callers merge same-named sets by passing the vector to extract APIs.
- `skin_sideset`: if exactly one non-empty per-block skin, return it; otherwise error pointing to `skin_sidesets()`.

#### 2.2 Multi-sideset / multi-block surface extract

- `create_surface_from_sidesets` concatenates same-type surfaces (volume node ids, orientation preserved).
- Mixed TRI/QUAD: `split_mixed_arity_sideset` / `create_surfaces_from_sidesets` (one buffer per shell type). Singular API prints and returns INVALID when multiple shell types remain.

#### 2.3 Hex-dominant / SSMIXED extractors

- Per-block family extract (`create_surface_from_sideset` / SS HEX/TET/QUAD/WEDGE/PYRAMID). No mixed-block SSMIXED extractor; `create_surfaces_from_sidesets` groups by shell type across blocks.

---

### Phase 3 — Mesh-integrated IO and registry (done)

#### 3.1 Optional Mesh sideset registry

- Mesh holds named `vector` of `(name, Sideset)` (non-breaking: empty by default).
- `SFC::reorder` / `reorder_elements_from_tags` remap registered sidesets (skip pointers also passed in the argument vector).
- Names are first-class (`add_sideset` / `sidesets(name)`). One name may own several sidesets (one per block).

#### 3.2 Mesh folder IO

- `Mesh::write`: emit `sidesets/<name>/` (or `sidesets/<name>/<block_id>/` when a name has multiple members) beside mesh data when the registry is non-empty.
- `Mesh::read` / distributed read: discover and load; MPI uses `read_and_redistibute`.
- Respects Phase 0.1 meta (`size`, `block_id`). `Mesh::clone` deep-copies the registry.

#### 3.3 Exodus / VTK (optional, lower priority)

- IO in SFEM is always using raw plus meta yaml files. No need for extra formats in the code
- External converter is used instead
  - Keep Python raw↔Exodus as supported path (includes blocks and sidesets from raw). 
  - raw_to_db should also support (when the target format supports it the inclusion of sideset and block information)

---

### Phase 4 — Element-type completeness (done)

#### 4.1 `LocalSideTable` coverage

- `LocalSideTable::fill` returns SUCCESS/FAILURE (no abort). Unsupported types print one diagnostic from create/skin/extract.
- EDGE2 / EDGE3 / EDGESHELL2 / EDGESHELL3: sides are the two endpoints (`NODE1`).
- Shell aliases: QUADSHELL4→QUAD4, TRISHELL3→TRI3, TRISHELL6→TRI6 (QUADSHELL9 already with QUAD9).
- PROTEUS higher-order stays on SS family extractors; `PROTEUS_HEX8` keeps its own corner layout.

#### 4.2 Higher-order face node layout

- HEX27 / TET10 / TRI6 / QUAD9 match Exodus/PATRAN (same as Python `raw_to_exodusII` HEX27). QUAD9 sides are EDGE3 (face center is not a side node).
- Contract tests in `sideset_test`.

---

### Phase 5 — Test and validation matrix

Add/extend tests so every Phase 0–4 item has a serial and (where applicable) MPI check.

| Test focus | Suggested location |
|------------|-------------------|
| MPI meta global size + block_id | `sideset_test` / new MPI IO test |
| Device round-trip | unit next to `smesh_device_sideset` |
| Reorder + live sideset remap | extend `skin_sideset_serial_parallel_test` |
| Unstructured refine map | new `sideset_refine_test` |
| SS invariance across levels | `parallel_ss_sideset_test`, hierarchies |
| Multiblock skin merge / names | `parallel_multiblock_topology_test` |
| Multi-sideset surface merge | frontend sideset test |
| WEDGE/PYRAMID auto-split extract | hex-dominant / wedge unit |
| Mesh read/write sidesets folder | `sideset_test` (serial + multiblock); `parallel_transforms_test` MPI redistribute |
| QUAD/TRI MPI skin parity | distributed skin test |

Parity rule for MPI tests: compare **sorted global face keys** (or global parent+lfi), not local indices.

---

## Execution order (recommended)

```
Phase 0 (meta, device, SS docs)
    → Phase 1.1 reorder remap
    → Phase 1.2 refine map
    → Phase 2 multiblock + multi-surface (with MULTIBLOCK.md §3)
    → Phase 3 registry + Mesh IO (done)
    → Phase 4 type coverage (done)
    → Phase 5 fill remaining coverage gaps continuously
```

Ship each phase with tests green under `build` (serial) and `build_mpi_debug` (MPI).

---

## Non-goals (for now)

- Changing the core `(parent, lfi)` representation to store explicit face connectivity.
- Keeping aura/ghost sideset entries after create/skin (owned-only remains the default; exchange APIs only if a consumer requires them).
- Full native Exodus/VTK C++ stack unless a concrete in-process requirement appears.

---

## Definition of done

Sidesets are fully functioning when:

1. **Types:** create, skin, surface, and nodeset work for every type listed as supported in `LocalSideTable`, including SS families and documented mixed/hex-dominant policy.
2. **MPI:** create/skin/write/redistribute preserve global uniqueness and correct ownership; meta size and device conversion are trustworthy.
3. **Refine:** unstructured refine/promote can map sidesets; SS level changes keep sidesets valid by contract + tests.
4. **Reorder:** element reorder updates sideset parents (via API and/or Mesh registry hooks).
5. **IO:** Mesh folder round-trips named sidesets in serial and MPI; custom sideset folder remains valid; Exodus via Python remains documented.

---

## Cross-references

- Multiblock create/skin/name merge: [MULTIBLOCK.md](MULTIBLOCK.md) §3, §5 (reorder), §6 (MPI IO)
- Deterministic cross-block ordering: MULTIBLOCK.md “Deterministic multi-block remapping rules”
