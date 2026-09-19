# FEATURES.md

Element-type coverage checklist. Sibling-family gaps only: a feature exists for some types and is missing for the natural counterparts. `NODE1` / `BEAM2` / `MACRO_*` are omitted unless they are the counterpart.

Marks: `[x]` done · `[ ]` missing · `[~]` partial · **I** intentional (not a parity bug).

Related: [2D.md](2D.md), [REFINE.md](REFINE.md), [PARAM.md](PARAM.md), [SIDESET.md](SIDESET.md), [MULTIBLOCK.md](MULTIBLOCK.md).

---

## `promote_to` (p-refine)

Serial and MPI, same-type multiblock. `elem_higher_order` already maps HEX8→HEX27 and EDGE2→EDGE3.

- [x] TET4 → TET10
- [x] TET4 → TET15
- [x] TRI3 → TRI6
- [x] TRISHELL3 → TRISHELL6
- [x] QUAD4 → QUAD9 (4 edge mids + unique face center)
- [x] QUADSHELL4 → QUADSHELL9
- [x] HEX8 → HEX27 (SS L=2 + PATRAN pointer reorder)
- [ ] EDGE2 → EDGE3 / EDGESHELL2 → EDGESHELL3
- [x] MPI (`comm_size > 1`; EDGE still missing)

**I:** mixed-type blocks rejected.

---

## `refine()` (h-refine)

Linear P1, serial + MPI. Higher-order input is refused (**I**).

- [x] HEX8
- [x] TET4 (Bey; not Kuhn `sstet_to_tet4`)
- [x] TRI3 / TRISHELL3
- [x] QUAD4 / QUADSHELL4
- [x] WEDGE6
- [x] PYRAMID5 (SS explode → PYRAMID5 + TET4 blocks)
- [x] EDGE2 / EDGESHELL2
- [x] Mixed HEX / TET / WEDGE / PYRAMID
- [x] Mixed QUAD4 + QUADSHELL4
- [ ] Mixed HEX + QUAD (**I:** no shared SS lattice)
- [ ] TET10 / TRI6 / QUAD9 / HEX27 / PROTEUS_* input (**I:** P1-only)

---

## Semistructured lattice

`to_semistructured` / `derefine`. Explode: `sshex_to_hex8`, `sstet_to_tet4`, `ssquad_to_quad4`, `sswedge_to_wedge6`, `ss_to_linear` (pyramid).

- [x] HEX
- [x] TET
- [x] QUAD / QUADSHELL
- [x] WEDGE
- [x] PYRAMID
- [x] Mixed HEX + TET
- [x] Hex-dominant (HEX + WEDGE + PYRAMID ± TET)
- [ ] TRI SS family (**I:** no `PROTEUS_TRI*`)
- [ ] Mixed QUAD + TRI
- [ ] GLL nodes except HEX SS (QUAD / TET / WEDGE / PYRAMID / mixed refuse)
- [ ] Dedicated `create_semistructured_{tet,wedge,pyramid}_*` (use `to_semistructured` on linear)

Factories:

- [x] `create_semistructured_hex_cube`
- [x] `create_semistructured_quad_square`
- [x] `create_square(PROTEUS_QUAD*)`
- [x] `create_cube(PROTEUS_HEX*)`

---

## `extrude`

Serial + MPI. Planar points padded to `z = 0`.

- [x] QUAD4 / QUADSHELL4 → HEX8
- [x] TRI3 / TRISHELL3 → WEDGE6
- [ ] EDGE* → quad strip
- [ ] Higher-order faces
- [ ] Mixed-type blocks (**I:** all blocks must share type)

---

## `convert_to`

- [x] HEX8 → TET4 (6×)
- [x] WEDGE6 → TET4 (3×)
- [x] PYRAMID5 → TET4 (2×)
- [x] QUAD4 → TRI3 (2×)
- [x] TET15 → HEX8 (4×)
- [ ] QUADSHELL4 → TRISHELL3
- [x] PROTEUS_HEX* → HEX8 (not `PROTEUS_HEX4913`)
- [x] PROTEUS_TET* → TET4 (including L=16)
- [x] PROTEUS_QUAD* / QUADSHELL* → linear (not L=16 `*289`)
- [x] PROTEUS_WEDGE* → WEDGE6
- [ ] PROTEUS_PYRAMID* (use `ss_to_linear`)
- [ ] PROTEUS_HEX4913 → HEX8
- [ ] PROTEUS_QUAD289 / PROTEUS_QUADSHELL289 → linear

---

## Jacobians / FFF

Compact 2D: adj SoA 4, FFF SoA 3. 3D: 9 / 6. Buffers sized by family.

- [x] TRI3 (planar)
- [x] QUAD4 (planar, qp `(1/2, 1/2)`)
- [x] TET4 / TET10
- [x] HEX8
- [x] HEX SS macros
- [ ] WEDGE6
- [ ] PYRAMID5
- [ ] HEX27
- [ ] TRISHELL* / QUADSHELL* surface J (3×2) (**I** this pass, [2D.md](2D.md))
- [ ] Non-HEX SS (`PROTEUS_QUAD*` / TET / WEDGE / PYRAMID): `adjugate_fill` / `fill_fff` always call `sshex8_macro_*`

---

## Restrict / prolongate

- [x] Host SS restrict: HEX
- [x] Host SS restrict: QUAD
- [x] Host SS restrict: TET (+ mixed HEX+TET)
- [ ] Host SS restrict: WEDGE / PYRAMID
- [ ] Mixed SS restrict with QUAD
- [x] Device SS restrict: HEX
- [ ] Device SS restrict: QUAD / TET / WEDGE / PYRAMID
- [x] Unstructured restrict: TET10 / MACRO_TET4 → TET4
- [x] Unstructured restrict: TRI6 / MACRO_TRI3 → TRI3
- [ ] Unstructured restrict: QUAD9 → QUAD4
- [ ] Unstructured restrict: HEX27 → HEX8
- [x] SS prolong kernels: HEX, TET, QUAD
- [ ] SS prolong kernels: WEDGE, PYRAMID
- [ ] Mesh-level `Prolongate` (stub)

Unstructured multi-block restrict is not implemented.

---

## Create dispatchers

`create_square`:

- [x] QUAD4
- [x] TRI3
- [x] TRI6 (via `promote_to`)
- [ ] QUAD9 (promote exists)
- [x] PROTEUS_QUAD*
- [ ] QUADSHELL* / PROTEUS_QUADSHELL* / TRISHELL*

`create_cube`:

- [x] HEX8
- [x] HEX27 (L=2 SS + node reorder, not `promote_to`)
- [x] TET4
- [x] TET10 (via `promote_to`)
- [x] PROTEUS_HEX*
- [ ] WEDGE6 / PYRAMID5
- [ ] PROTEUS_TET* / PROTEUS_WEDGE* / PROTEUS_PYRAMID*

Specialty (HEX-only is **I** unless noted):

- [x] `create_quad4_ring` (serial + MPI)
- [ ] TRI / HEX ring
- [x] `create_half_sphere` HEX8 / TET4
- [x] HEX nozzle / L-shape / checkerboard / bidomain / hex-dominant
- [x] MPI create for nozzle / L-shape
- [x] Mixed HEX+TET cube / hex-dominant cylinder
- [ ] Mixed QUAD+TRI factory

---

## Sidesets / edgesets / nodesets

Skin, selector create, and SS extract are in place for HEX/TET/TRI/QUAD/WEDGE/PYRAMID and their SS families ([SIDESET.md](SIDESET.md)).

`map_sideset_through_refine`:

- [x] HEX8
- [x] QUAD4 / QUADSHELL4
- [x] WEDGE6
- [x] TET4
- [x] TRI3 / TRISHELL3
- [x] EDGE2 / EDGESHELL2
- [~] PYRAMID5 (quad base `lfi==4` only; tri faces `lfi` 0–3 missing)

`map_edgeset_through_refine`:

- [x] HEX / TET / TRI / QUAD / WEDGE / EDGE
- [ ] PYRAMID5

Nodesets through refine/promote insert mid-edge nodes when both coarse endpoints are in the set (face/body nodes are not added). SS `(parent, lfi)` not remapped (**I**).

---

## SFC / parametrization / 2D contract

SFC:

- [x] Single-block, `sdim==2` (zero-z into 3D encoders)
- [x] Multiblock (`n_blocks() > 1`: union encoding, per-block permute, first-touch nodes; SS hierarchical numbering remains `n_blocks()==1` only)
- [x] MPI (`is_distributed()`: global-bbox encode, permute within owned/shared/aura; nodes and ghost import maps unchanged)

Parametrization is nodeset-based, not type-gated ([PARAM.md](PARAM.md)):

- [x] Circle, including planar `sdim==2`
- [x] Sphere / polynomial surface (`sdim>=3` required)
- [x] Identity

Planar `sdim==2` ([2D.md](2D.md)):

- [x] 2-component points, barycenters, SFC, MPI extrude pad
- [x] TRI3 / QUAD4 J/FFF
- [x] In-plane circle
- [x] Shells stay 3-component points

---

## Geometric map (per block)

One `GeomMap` on each `Mesh::Block`: `Affine`, `AxisAligned` (HEX/QUAD families only, including shells and Proteus), `IsoParametric` (default).

- [x] Stored on the block; `AXIS_ALIGNED` rejected on other families
- [x] Cartesian HEX/QUAD factories set `AxisAligned`; tet/tri cubes `Affine`; nozzle / ring / sphere-hex / hump-hex / cylinder-hex `IsoParametric`; linear tet sphere/hump stay `Affine`
- [x] clone / split / convert / promote / extrude / SS inherit (`AxisAligned` → `Affine` when the destination is not HEX/QUAD)
- [x] `meta.yaml` `geom_map:` round-trip (missing key = `IsoParametric`)
- [x] `detect_geom_map` / `detect_and_set_geom_map` from coordinates (MPI AND; empty ranks do not constrain)
