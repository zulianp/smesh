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

## refine

```text
refine <mesh_folder> <output_folder>
```

Levels: environment `SMESH_REFINEMENT_LEVELS` (default `1`). Linear types only.

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
