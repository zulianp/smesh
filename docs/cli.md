# Command-line tools

Binaries land in the CMake build directory. Paths below assume `./build/`.

## cube

```text
cube <element_type> <nx> <ny> <nz> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> <output_folder>
```

Example: `./build/cube HEX8 8 8 8 0 0 0 1 1 1 hex_cube`

## square

```text
square <element_type> <nx> <ny> <xmin> <ymin> <xmax> <ymax> <output_folder>
```

Planar mesh (`TRI3`, `QUAD4`, …).

## nozzle

FDA HEX8 nozzle (mm). `n_core` even and ≥ 2. `expansion < 0` is a pipe (no sudden expansion).

```text
nozzle <n_core> <n_bore> <n_outer> <expansion> <n_axial> <output_folder>
```

`n_axial` is one integer (uniform) or four comma-separated counts (graded), e.g. `6,5,10,30`. Optional `SMESH_NOZZLE_CORE_FRACTION` (default `0.5`).

```bash
./build/nozzle 2 1 1 3 6,5,10,30 fda_nozzle
./build/nozzle 2 1 1 -1 8 pipe
```

## lshape

```text
lshape <nx> <ny> <nz> <xmax> <ymax> <zmax> <step_x> <step_y> <output_folder>
```

HEX8 backward-facing step. The step must fall on a grid line.

## hump

```text
hump <element_type> <nx> <ny> <nz> <length> <height> <width> <hump_start> <hump_length> <hump_height> <output_folder>
```

`HEX8` or `TET4`. Defaults of the factory are `32 12 4 9 3 1 0.65 1 0.128`.

## hex_tet_cube

```text
hex_tet_cube <nx> <ny> <nz> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> <output_folder>
```

HEX8 cube; the second half of hexes becomes 6 TET4 each.

## hex_dominant

```text
hex_dominant <output_folder>
```

Fixed serial unit: HEX8 + PYRAMID5 + WEDGE6 + TET4. For a cylinder, use `cylinder`.

## sshex_cube / ssquad_square

```text
sshex_cube <micro_per_dim> <nx> <ny> <nz> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> <output_folder>
ssquad_square <micro_per_dim> <nx> <ny> <xmin> <ymin> <xmax> <ymax> <output_folder>
```

## Other create binaries

`cylinder`, `half_sphere`, `ring2`, `checkerboard_cube`, `bidomain_cube`.

## refine

```text
refine <mesh_folder> <output_folder>
```

Levels: environment `SMESH_REFINEMENT_LEVELS` (default `1`). Linear types only.

## adapt_refine

```text
adapt_refine <mesh_folder> <output_folder>
```

Curvature-driven conforming adapt for single-block serial `TET4` and `TRI`/`QUAD` (including shells). Not uniform `refine()`.

Environment:

- `SMESH_CELLS_PER_RADIUS` (default `8`)
- `SMESH_SHARP_EDGES_THRESHOLD` (default `0.15`)
- `SMESH_ADAPT_LEVELS` (default `8`)
- `SMESH_SMOOTH_ITERS` (default `10`)

## improve

```text
improve <mesh_folder> <output_folder>
```

Quality remesh for single-block serial `TET4` and `TRI`/`QUAD` (including shells). Not uniform `refine()` and not `adapt_refine()`.

Environment:

- `SMESH_Q_MIN` (default `0.3`)
- `SMESH_MAX_ABS_DEV` (default `0` → `0.02 * bbox_diag`)
- `SMESH_MAX_NORMAL_DEV` (default `0` → same as abs)
- `SMESH_SHARP_EDGES_THRESHOLD` (default `0.15`)
- `SMESH_IMPROVE_PASSES` (default `8`)
- `SMESH_SMOOTH_ITERS` (default `8`)

## mesh_promote

```text
mesh_promote <to_element_type> <input_mesh> <output_mesh>
```

Example: `./build/mesh_promote TET10 tet_cube tet10`

## mesh_convert

```text
mesh_convert <to_element_type> <input_mesh> <output_mesh>
```

Topology conversion (for example HEX8 → TET4), not p-refinement.

## skin

```text
skin <mesh_folder> <output_folder>
```

Surface mesh plus `parent_sideset` in the output folder.

## create_sideset

```text
create_sideset <mesh_folder> <x> <y> <z> <angle_threshold> <output_folder>
```

Starts from the skin at `(x,y,z)`. `angle_threshold` is a cosine bound in `[0, 1]`.

The tree has more drivers (graphs, ordering, extrude, …). They are not required to use the library.
