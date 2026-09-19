# REFINE.md

## Goal

Make unstructured `refine()` a first-class, mesh-lifecycle-safe h-refinement for every linear family that already has a closed same-type split (or an SS lattice that explodes to one), under serial and MPI, including multiblock, sidesets/edgesets/nodesets, and IO.

This plan is functional-first. SS lattice conversion (`to_semistructured` / `derefine`) is a related but distinct path; HEX `refine()` already uses it as an implementation detail. Sideset remap through refine overlaps [SIDESET.md](SIDESET.md) Phase 1.2; this document owns element-type completeness and the refine mesh transform itself.

## Scope

- No breaking public API signatures unless explicitly versioned.
- Single-block serial HEX/TET/TRI behavior remains the compatibility baseline (child count and type).
- Unsupported combinations must fail with explicit diagnostics (operation, block, type, reason, workaround).
- `refine()` stays **same-type h-refinement**: more elements of the same `ElemType`, midpoint geometry.
- `promote_to` stays **p-refinement** (same elements, extra nodes). Do not fold it into `refine()`.
- `to_semistructured` / `derefine` stay **SS lattice** (macro element count unchanged). Do not make `refine()` return PROTEUS_* types.

---

## Current state (summary)

Three operations look like “refinement” and must not be confused:

| API | What it does | Element count | Type out |
|-----|----------------|---------------|----------|
| `refine(mesh, levels)` | Unstructured h-refine | × 8^ℓ (HEX/TET/WEDGE) or × 4^ℓ (TRI/QUAD) or × 2^ℓ (EDGE) | same linear type |
| `to_semistructured(L, mesh)` | Embed each element in an SS lattice | unchanged macros | PROTEUS_* |
| `derefine(ss, to_level)` | Drop SS lattice nodes | unchanged macros | coarser PROTEUS_* |
| `promote_to(TET10/TRI6/TET15)` | p-refine | unchanged | higher-order |

### `refine()` dispatch

Entry: `refine()` in `src/frontend/smesh_mesh.cpp`. Distributed meshes go to `MeshTransformsDistributed::refine` (`src/frontend/smesh_mesh_transforms_distributed.cpp`).

Hard stops today:

- Mixed types other than HEX8+WEDGE6 and QUAD4+QUADSHELL4 are rejected (HEX+TET is Bey vs Kuhn; HEX+QUAD has no shared SS lattice; PYRAMID is SS-only).
- Type must be `HEX8`, `TET4`, `TRI3`, `TRISHELL3`, `QUAD4`, `QUADSHELL4`, `WEDGE6`, `EDGE2`, or `EDGESHELL2`, or mixed HEX+WEDGE / QUAD4+QUADSHELL4.
- Without MPI, `comm_size > 1` is rejected.

| Family | Serial | MPI | Mechanism | Children | Same-type multiblock |
|--------|--------|-----|-----------|----------|----------------------|
| HEX8 | Yes | Yes | `to_semistructured(2^levels)` then `sshex_to_hex8` (one shot) | L³ with L = 2^levels | Yes; mixed HEX+WEDGE via hex-dominant SS + `ss_to_linear` |
| TET4 | Yes | Yes | Edge-midpoint Bey/octahedron 8-split, iterated `levels` times | 8 per level | Yes (same-type). Mixed HEX+TET rejected (Bey ≠ Kuhn) |
| TRI3 / TRISHELL3 | Yes | Yes | Edge-midpoint 4-split, iterated `levels` times | 4 per level | Yes (same-type) |
| QUAD4 / QUADSHELL4 | Yes | Yes | `to_semistructured(2^levels)` then `ssquad_to_quad4` (one shot) | L² with L = 2^levels | Yes; mixed QUAD4+QUADSHELL4 allowed |
| WEDGE6 | Yes | Yes | `to_semistructured(2^levels)` then `sswedge_to_wedge6` (one shot) | L³ with L = 2^levels | Yes; mixed HEX+WEDGE as above |
| EDGE2 / EDGESHELL2 | Yes | Yes | Edge-midpoint 2-split, iterated `levels` times | 2 per level | Yes |

Kernel: `mesh_refine` (`src/mesh/refinement/smesh_refine.impl.hpp`) implements **TET4, TRI3/TRISHELL3, and EDGE2/EDGESHELL2**. HEX/QUAD/WEDGE never enter that kernel. MPI TET/TRI/EDGE uses the same child-pattern tables in `refine_edges_once` (global edge unique + midpoint).

HEX child layout is the SSHEX lexicographic microhex order (`sshex8_to_standard_hex8_mesh`). QUAD child layout is `ssquad4_to_standard_quad4_mesh` (`le = yi*L + xi`). WEDGE child layout is `sswedge_to_standard_wedge6_mesh` (up/down microtriangles × layers; `le` matches that enumeration). TET child layout is the 8-row `tet4_refine_pattern` (4 corners + 4 octahedron tets). These are **not** the same as `sstet4_to_standard_tet4_mesh` (Kuhn/HyTeG L³ tets), even though L=2 also yields 8 children. TRI child layout is `tri3_refine_pattern`. EDGE children are `{n0, mid}` then `{mid, n1}`.

### SS lattice (already implemented, not wired to `refine()`)

`to_semistructured` / `derefine` already cover HEX, TET, QUAD/QUADSHELL, WEDGE, PYRAMID, mixed HEX+TET, and hex-dominant (WEDGE/PYRAMID in the mix). Serial + MPI tests live in `multiblock_ss_packed_test` and `parallel_to_semistructured_test`.

Explode SS → linear (needed to implement HEX-style `refine()`):

| Explode | Serial | MPI attach | Used by `refine()` |
|---------|--------|------------|--------------------|
| `sshex_to_hex8` | Yes | Yes (`attach_sshex_to_hex8` via `conversion_factor`) | HEX8 |
| `sstet_to_tet4` | Yes | Yes (same attach; L³) | No |
| `ssquad_to_quad4` | Yes | Yes (same attach; L²); preserves QUADSHELL | QUAD4 / QUADSHELL4 |
| `sswedge_to_wedge6` | Yes | Yes (same attach; L³) | WEDGE6 |
| PYRAMID | **No** | No | No |

### `promote_to` (p-refine, related)

Serial and MPI. Maps: TET4→TET10, TET4→TET15, TRI3→TRI6, TRISHELL3→TRISHELL6, QUAD4→QUAD9 (4 edge mids + unique face center), QUADSHELL4→QUADSHELL9, HEX8→HEX27 (SS L=2 lattice + PATRAN SoA pointer reorder). Same-type multiblock works. EDGE2→EDGE3 is **not** implemented.

### Sets through refine

| Object | Behavior today |
|--------|----------------|
| Mesh sideset registry | Remapped after `refine()` via `map_sideset_through_refine` (fail if non-empty and unsupported) |
| `map_sideset_through_refine` | HEX8 (face → L² children), QUAD4/QUADSHELL4 (edge → L children), WEDGE6 (TRI or QUAD face → L² children, same `lfi`), TET4 (face → 4 children/level; child `lfi` from `tet4_face_child_lfi`, not always the coarse face — octahedron child 7 on face 2 uses lfi 3), TRI3/TRISHELL3 (edge → 2 children/level), EDGE2/EDGESHELL2 (NODE1 endpoint → 1 child at that end); factor 1 is identity |
| Edgesets | `map_edgeset_through_refine` for HEX/TET/TRI/TRISHELL/QUAD/QUADSHELL/WEDGE/EDGE (`lei` unchanged; HEX/QUAD/WEDGE: L children, TET/TRI/EDGE: 2 per level). Mesh registry remapped after `refine()` |
| Nodesets | Coarse ids copied, then mid-edge closure: if both endpoints of a coarse edge are in the set, insert every new node on that edge (not face/body). TET/TRI/EDGE: per-level CRS/child mids. HEX/QUAD/WEDGE: SS lattice interiors. MPI remaps local ids by GID. Mesh registry remapped after `refine()` |
| SS `(parent, lfi)` | Intentionally not remapped; mapper prints and returns nullptr |

### Drivers / tests

- `refine.exe` / `derefine.exe`: `refine.exe` uses `initialize` (serial or MPI); `derefine.exe` is serial (`initialize_serial`). `SMESH_REFINEMENT_LEVELS` / `SMESH_DEREFINEMENT_LEVELS`.
- Tests: serial HEX/TET/QUAD/WEDGE/TRISHELL/EDGE multiblock + mixed HEX+WEDGE / QUAD4+QUADSHELL4 + HEX+TET/HEX+QUAD/PYRAMID/TET10/TRI6/SS reject (`parallel_multiblock_topology_test`); serial sideset map HEX/TET/TRI/TRISHELL/QUAD/QUADSHELL/WEDGE/EDGE + EDGE edgeset + Mesh registry HEX/TET/TRI/QUAD + HEX27/TET10 promote sets (`sideset_test`); MPI HEX/TET/TRI/TRISHELL/QUAD/WEDGE/EDGE + HEX registry + mixed HEX+WEDGE + PYRAMID/hex-dominant/TET10/SS reject + `promote_to` TET10/TET15/TRI6/TRISHELL6/QUAD9/QUADSHELL9/HEX27 + HEX promote/refine nodeset GIDs + multiblock TET10 (`parallel_transforms_test`). Converters (`db_to_raw` / `raw_to_db`) stay translation-only.

### Core files

- Unstructured kernel: `src/mesh/refinement/smesh_refine.*`
- Frontend: `src/frontend/smesh_mesh.cpp` (`refine`), `smesh_mesh_transforms_distributed.cpp`
- SS: `src/frontend/smesh_semistructured.*`, `src/mesh/semistructured/smesh_ss{hex8,tet4,quad4,wedge,pyramid}*`
- Sideset map: `src/frontend/smesh_sideset.cpp` (`map_sideset_through_refine`)
- Mixed explode: `ss_to_linear` in `smesh_semistructured.cpp`
- p-refine: `src/mesh/conversions/smesh_promotions.*`, `promote_to` in `smesh_mesh.cpp`
- Drivers: `src/drivers/refinement/refine.exe.cpp`, `derefine.exe.cpp`

---

## Gap matrix

### Unstructured `refine()` by family

| Family | `refine()` | Natural split | SS explode available | Gap |
|--------|------------|---------------|----------------------|-----|
| HEX8 | Yes (SS explode) | 8^ℓ hexes | Yes + MPI | — |
| TET4 | Yes (Bey 8-split) | 8^ℓ tets | Yes, **different** connectivity | Keep Bey as contract; do not silently switch to Kuhn |
| TRI3 / TRISHELL3 | Yes (4-split) | 4^ℓ tris | N/A (no TRI SS family) | — |
| QUAD4 / QUADSHELL4 | Yes (SS explode) | 4^ℓ quads (L=2^ℓ) | Yes + MPI | — |
| WEDGE6 | Yes (SS explode) | 8^ℓ wedges (L=2^ℓ) | Yes + MPI | — |
| PYRAMID5 | **Yes** (SS explode, 2 output blocks) | sspyramid_n_pyr(L) pyramids + sspyramid_n_tet(L) tets | Yes + MPI | Emits a sibling TET4 block `{name}_tets`. L=1 → identity (no tet block). |
| EDGE2 / EDGESHELL2 | Yes (2-split) | 2^ℓ edges | N/A | — |
| HEX27 / TET10 / TRI6 / QUAD9 / TET15 / TET20 | **No** | Demote then h-refine, or refuse | — | Refused (P1-only). Demote is not a `refine()` path. |
| PROTEUS_* (SS) | **No** | Use `to_semistructured` from linear, or explode then `refine()` | — | Refused. Explode (`sshex_to_hex8` / …) then `refine()`, or regenerate SS at a higher L. |
| Mixed-type blocks | HEX/TET/WEDGE/PYRAMID (any mix), QUAD4+QUADSHELL4 | Shared SS node pool then per-block explode | HEX+QUAD rejected (no mixed SS); HEX+TET now accepted (Kuhn TET); PYRAMID creates extra TET4 block |
| NODE1 / BEAM2 / MACRO_* | **No** | — | — | Stay unsupported with diagnostic |

Hard stops that should stay explicit until a phase lands:

- `refine()`: HEX+QUAD (no shared SS lattice for volume+surface), other unsupported mixes

### Parallel / MPI

| Item | Status | Gap |
|------|--------|-----|
| HEX8 `refine` | Done | Child GIDs via SS attach |
| TET4 / TRI3 `refine` | Done | MPI parity vs serial in `parallel_transforms_test` |
| QUAD `refine` | Done | MPI attach via `conversion_factor` L² |
| WEDGE `refine` | Done | Explode + attach factor L³ |
| Auto-map Mesh sideset registry | Done | `refine()` remaps sidesets/edgesets/nodesets |
| Mixed HEX+WEDGE `refine` | Done | Hex-dominant SS + `ss_to_linear`; MPI attach factor L³ |
| Mixed HEX/TET/WEDGE/PYRAMID `refine` | Done | SS lattice + per-block explode; PYRAMID → PYRAMID5+TET4; Kuhn TET; MPI attach via `attach_ss_to_linear` |
| `refine.exe` MPI | Done | `initialize` + library `refine()` |

### Related transforms (do not block `refine()` completeness)

| Path | Status | Note |
|------|--------|------|
| SS HEX/TET/QUAD/WEDGE/PYRAMID + mixed | Done | Not a substitute for unstructured `refine()` |
| `sstet_to_tet4` | Done | Optional alternate TET split; would break current sideset child table |
| `promote_to` HEX/EDGE | Missing | Separate workstream |
| `promote_to` MPI | Done | Unique edge/face/elem nodes via `unique_inc_tuples`; HEX/EDGE still missing |

---

## Work plan

### Phase 0 — Contract, diagnostics, test holes (done)

Lock the current HEX/TET/TRI contract so later families do not change it.

#### 0.1 Documented contract (this file)

- HEX: one-shot SS L=2^levels, L³ children, SSHEX microhex order.
- TET: iterated Bey 8-split (`tet4_refine_pattern`), **not** `sstet4_to_standard_tet4_mesh`.
- TRI: iterated 4-split (`tri3_refine_pattern`).
- SS input, higher-order, mixed types: fail with type name + workaround (`to_semistructured`, `sshex_to_hex8`, `sswedge_to_wedge6`, `convert_to`).

#### 0.2 Diagnostics

- `mesh_refine` prints `type_to_string` and returns `SMESH_FAILURE` (TET4/TRI3 kernel only).
- `refine()` / distributed `refine()`: `refine_print_unsupported` / `refine_print_mixed_types` then `nullptr` (no abort). Workarounds for PYRAMID, SS, P1-only.
- Removed stale “TODO fill p2 node indices” comments in `smesh_refine.impl.hpp`.

#### 0.3 Tests for what already exists

- MPI TET4 and TRI3 `refine(levels=1)` vs serial node/element global counts (+ sideset map) in `parallel_transforms_test`.
- Serial TRISHELL3 negative was replaced in Phase 3 by `test_serial_trishell3_refine`.

No new element types.

---

### Phase 1 — QUAD4 / QUADSHELL4 unstructured `refine()` (done)

Mirror HEX: `to_semistructured(2^levels)` then explode to linear quads.

#### 1.1 MPI attach for `ssquad_to_quad4`

- `attach_sshex_to_hex8` expands by `conversion_factor` (HEX L³, TET L³, QUAD L², WEDGE L³).
- `ssquad_block_to_quad4_block` / `ssquad_linear_type` keep QUADSHELL vs QUAD.

#### 1.2 Wire `refine()`

- Serial + MPI, single- and multi-block same-type.
- Child factor 4^levels. Sideset map: each QUAD side (EDGE) → L children, same `lfi` (`quad_edge_child_on`, child id `yi*L+xi`).

#### 1.3 Tests

- Serial square vs SS explode: `test_serial_quad4_refine`.
- Multiblock vs single: `test_serial_quad4_multiblock_refine`.
- QUADSHELL type preserved: `test_serial_quadshell4_refine`.
- MPI vs serial: `test_mpi_quad4_refine`.
- `map_sideset_through_refine`: `test_sideset_map_quad_refine`, `test_sideset_map_quadshell_refine`.

---

### Phase 2 — WEDGE6 unstructured `refine()` (done)

Same-type closed split: triangle 4-split × layer split = L³ wedges, L=2^levels.

#### 2.1 Explode

- `sswedge_to_wedge6` / `sswedge_block_to_wedge6_block` (up/down microtriangles × z layers).
- MPI attach via `conversion_factor` L³ (`sswedge_txe`).

#### 2.2 Wire `refine()` and sideset map

- Serial + MPI, single- and multi-block same-type.
- TRI faces (lfi 3,4) → L² children; QUAD faces (lfi 0–2) → L² children; same `lfi` (`wedge_child_on`).

#### 2.3 Tests

- Kernel explode: `test_sswedge_to_standard_wedge6`.
- Serial vs SS explode: `test_serial_wedge6_refine`.
- Multiblock vs single: `test_serial_wedge6_multiblock_refine`.
- MPI vs serial: `test_mpi_wedge6_refine`.
- Sideset map QUAD face: `test_sideset_map_wedge_refine_quad_face`; TRI face: `test_sideset_map_wedge_refine_tri_face`.

---

### Phase 3 — Thin aliases and EDGE (done)

#### 3.1 Shell aliases in `refine()`

- TRISHELL3 → same kernel as TRI3 (`tri3_refine_pattern`), keep TRISHELL3.
- QUADSHELL already in Phase 1.

#### 3.2 EDGE2 / EDGESHELL2

- One midpoint per edge, 2 children `{n0, mid}` / `{mid, n1}`. Serial: n2n upper-triangular + `mesh_refine`. MPI: `refine_edges_once` edge unique.
- Sideset: NODE1 endpoints stay size 1 (start → first child, end → last child of the 2^ℓ block).
- Edgeset: `map_edgeset_through_refine`, `lei` 0 → 2 children per level.

Skip BEAM2 unless a consumer appears (same topology as EDGE2 but different type id).

#### 3.3 Tests

- Serial TRISHELL3 type preserved: `test_serial_trishell3_refine`.
- TRISHELL3 multiblock: `test_serial_trishell3_multiblock_refine`.
- Serial EDGE2 / EDGESHELL2: `test_serial_edge2_refine`, `test_serial_edgeshell2_refine`.
- EDGE2 multiblock: `test_serial_edge2_multiblock_refine`.
- MPI vs serial: `test_mpi_trishell3_refine`, `test_mpi_edge2_refine`.
- Sideset map TRISHELL: `test_sideset_map_trishell_refine`; EDGE NODE1: `test_sideset_map_edge_refine`.
- Edgeset map EDGE: `test_edgeset_map_edge_refine`.

---

### Phase 4 — Mesh-owned sets through `refine()` (done)

Callers should not have to remember `map_sideset_through_refine`.

#### 4.1 Sidesets

- After a successful `refine()`, remap the Mesh registry with the existing mapper.
- Empty registry stays empty. If the registry is non-empty and the mapper does not support the pair, `refine()` fails (`nullptr`).

#### 4.2 Edgesets / nodesets

- Edgesets: parent expansion like sidesets (`lei` child table) for HEX/TET/TRI/QUAD/WEDGE/EDGE. `lei` is unchanged.
- Nodesets: coarse node ids remain valid on serial meshes (refine only **adds** nodes). MPI remaps local ids by GID. Mid-edge nodes are inserted when both coarse endpoints are already in the set.

#### 4.3 Tests

- HEX/TET/TRI/QUAD registry round-trip: `test_refine_registry_*` in `sideset_test`.
- TET Bey child lfi (octahedron face 2 uses lfi 3): `test_tet4_face_child_lfi_matches_pattern`, `test_sideset_map_tet_refine_skin`.
- Empty registry stays empty: `test_refine_registry_empty`.
- MPI HEX registry: `test_mpi_hex8_refine` registers `left` / `left_nodes` and checks the fine mesh.

---

### Phase 5 — Mixed-type blocks (done)

`refine()` no longer requires one type for every combination.

- HEX8+WEDGE6: one `to_semistructured(2^levels)` (hex-dominant SS, shared node pool) then `ss_to_linear` (per-block explode). Serial + MPI. Child factor L³ on both families.
- QUAD4+QUADSHELL4: same SS family; `to_semistructured` then `ssquad_to_quad4` (type preserved per block).
- HEX+TET: now accepted via `refine_mixed_volume_ss`. TET children are Kuhn (L³ per macro-tet), not Bey. Same-type TET still uses Bey; this Kuhn path is only for mixed-volume meshes.
- HEX/WEDGE+QUAD: **rejected** (no mixed volume+surface SS lattice). Split by family.
- PYRAMID in the mix: now accepted (Phase 9). Pyramid SS block → PYRAMID5 block + new TET4 `{name}_tets` block.

#### Tests

- Serial HEX+WEDGE vs SS explode + sideset registry: `test_serial_hex_wedge_refine`.
- Serial QUAD4+QUADSHELL4: `test_serial_quad_quadshell_refine`.
- Serial HEX+TET now succeeds (Kuhn TET): `test_serial_hex_tet_refine_rejected` updated.
- Serial HEX+QUAD reject: `test_serial_hex_quad_refine_rejected`.
- MPI HEX+WEDGE vs serial counts: `test_mpi_hex_wedge_refine`.

---

### Phase 6 / Phase 9 — PYRAMID5 (done: SS explode, 2 output blocks)

A pyramid does not tile into pyramids only. The SS lattice is used and the children are split into two blocks:

- `refine()` on PYRAMID5 (or any mix containing PYRAMID5) goes through `to_semistructured(L)` then `ss_to_linear`.
- The PYRAMID SS block explodes into a **PYRAMID5 block** (upward + inverted pyramids) and a **new TET4 block** (`{name}_tets`, seam tets).
- L=1 is identity (1 pyramid, 0 tets → no extra block).
- Layer stencil (for each layer k=0..L-1, s=L-k): s² upward pyramids, (s-1)² inverted pyramids, 2s(s-1) seam tets.
- Count formulas: `sspyramid_n_pyr(L) = L² + L(L-1)(2L-1)/3`; `sspyramid_n_tet(L) = 2L(L²-1)/3`.
- L=2: 6 pyramids + 4 tets on 14 lattice nodes. L=4: 44 + 40.
- `convert_to(TET4)` still splits each macro pyramid into 2 tets (`mesh_pyramid5_to_2x_tet4`), unchanged.

Sideset remap: quad base face (lfi=4) on PYRAMID5 → L² child upward pyramids. Tri faces not yet remapped (returns nullptr gracefully with a diagnostic).

#### Tests

- Serial PYRAMID: `test_serial_pyramid_refine` — 2 output blocks, correct counts, `to_semistructured` still works.
- MPI vs serial pyramid strip: `test_mpi_pyramid_refine` — `n_nodes_global` / `n_elements_global` match serial; per-block owned ×6 / ×4; owned node GIDs unique.
- MPI vs serial hex-dominant: `test_mpi_hex_dominant_refine` — 4→5 blocks, global counts match serial, per-block factors 8/6/4/8/8.
- Hex-dominant SS + refine in `parallel_to_semistructured_test` also compares serial vs MPI globals.
- Count/kernel unit tests: `test_sspyramid_n_pyr_n_tet_counts`, `test_sspyramid_explode_l2`, `test_sspyramid_explode_l1_identity` in `smesh_sswedge_pyramid_mesh_test`.

---

### Phase 7 — Higher-order and SS input (done)

Diagnostics for PROTEUS_* and P1-only types landed in Phase 0 (`refine_print_unsupported`). Locked:

- `refine()` stays P1-only. TET10/TET15/TRI6/HEX27 and other higher-order types print a diagnostic and return `nullptr`. Dropping mid-nodes then h-refine is **lossy** and is **not** a `refine()` path. Use the linear parent, or `promote_to` after `refine()`.
- SS / PROTEUS_* input is refused. `refine()` does not return PROTEUS_* types. Workaround: explode to linear (`sshex_to_hex8` / `sstet_to_tet4` / `ssquad_to_quad4` / `sswedge_to_wedge6`) then `refine()`, or regenerate SS from the linear parent at a higher L.
- Do not h-refine higher-order connectivity in-place.

#### Tests

- Serial TET10/TET15/TRI6 and SS HEX `refine()` → `nullptr`; explode-then-refine still works: `test_serial_higher_order_ss_refine_rejected`.
- MPI SS HEX and TET10 `refine()` → `nullptr`: `test_mpi_higher_order_ss_refine_rejected`.

---

### Phase 8 — Kernel hygiene (done)

- One child-pattern table in `smesh_refine.hpp`: `tet4_refine_pattern` / `tri3_refine_pattern` / `edge2_refine_pattern`, plus Exodus edge lists and the TET/TRI face/edge-child maps. Used by serial `mesh_refine`, MPI `refine_edges_once`, sideset/edgeset mappers, and MACRO_TET4 n2e.
- HEX/QUAD/WEDGE stay SS-mediated; not reimplemented in `mesh_refine`.
- `refine.exe` uses `initialize` (MPI-capable). Null mesh/`refine()` is a hard failure, not a crash.

---

## Execution order (recommended)

```
Phase 0  contract + MPI TET/TRI tests + diagnostics (done)
    → Phase 1  QUAD/QUADSHELL refine (SS explode + MPI attach + sideset map) (done)
    → Phase 2  WEDGE explode + refine + sideset map (done)
    → Phase 3  TRISHELL3 / EDGE2 (done)
    → Phase 4  Mesh registry / edgeset / nodeset through refine (done)
    → Phase 5  mixed blocks HEX+WEDGE and QUAD4+QUADSHELL4 (done)
    → Phase 6  PYRAMID: keep unsupported in refine() (SS-only; done → superseded by Phase 9)
    → Phase 7  higher-order / SS input: refuse (no demote; done)
    → Phase 8  dedupe kernels + MPI driver (done)
    → Phase 9  mixed-volume SS: PYRAMID+TET+HEX+WEDGE in one lattice; pyramid → 2 output blocks (done)
```

Ship each phase with tests green under `build` (serial) and `build_mpi_debug` (MPI).

---

## Non-goals (for now)

- Changing HEX child ordering or replacing TET Bey split with Kuhn (`sstet_to_tet4`).
- Returning PROTEUS_* from `refine()`.
- Adaptive / local (hanging-node) refinement.
- EDGE2→EDGE3 p-refine (separate from h-refine). `promote_to` MPI covers TET10/TET15/TRI6/TRISHELL6/QUAD9/QUADSHELL9/HEX27.
- PYRAMID tri-face sideset remap through refine (quad base is supported; tri sides are not yet remapped).
- Hanging-node adaptivity.
- Higher-order demote (drop mid-nodes) inside `refine()`.

---

## Definition of done

Unstructured refine is complete when:

1. **Types:** `refine()` works for HEX8, TET4, TRI3, QUAD4/QUADSHELL4, WEDGE6, TRISHELL3, EDGE2/EDGESHELL2, PYRAMID5 (2-block output), and any mix of HEX/TET/WEDGE/PYRAMID, serial and MPI.
2. **Rejected clearly:** mixed HEX+QUAD (no shared SS lattice for volume+surface), higher-order, SS input, mixed unsupported families.
3. **Sets:** Mesh sideset registry remaps for supported pairs; edgesets remap; nodesets keep coarse ids and add mid-edge nodes when both endpoints were in the set.
4. **Parity:** every MPI refine path compares serial vs `mpiexec -np 2` global node/element counts (and owned uniqueness). Child layout matches the sideset mapper.
5. **SS remains separate:** `to_semistructured` / `derefine` still own lattice levels; HEX/QUAD/WEDGE `refine()` may use them internally but output is linear.

---

## Cross-references

- Sideset remap through unstructured refine: [SIDESET.md](SIDESET.md) Phase 1.2, “Refinement / hierarchies”
- Multiblock `Mesh::refine`: [MULTIBLOCK.md](MULTIBLOCK.md) §1
- SS mixed conversion: [MULTIBLOCK.md](MULTIBLOCK.md) §4
