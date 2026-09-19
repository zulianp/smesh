# Getting started

Assume a serial build in `build/` ([install.md](install.md)).

## Create a cube

```bash
./build/cube HEX8 8 8 8 0 0 0 1 1 1 hex_cube
```

`hex_cube/` is a folder mesh: `meta.yaml`, connectivity `i0`–`i7`, coordinates `x`/`y`/`z`. See [format.md](format.md).

`TET4` works the same way (six tets per hex cell).

## Refine and promote

```bash
./build/refine hex_cube hex_refined
SMESH_REFINEMENT_LEVELS=2 ./build/refine hex_cube hex_refined_2

./build/mesh_promote HEX27 hex_cube hex27
```

`refine` is h-refinement (more linear elements). `mesh_promote` is p-refinement (same elements, extra nodes). They are not interchangeable.

## Convert type

```bash
./build/mesh_convert TET4 hex_cube hex_as_tets
```

## Export to VTK or Exodus

Converters are translation-only. They do not refine or promote.

```bash
python3 python/smesh/raw_to_db.py hex_cube hex_cube.vtk
python3 python/smesh/db_to_raw.py hex_cube.vtk hex_cube_roundtrip
```

Exodus (`.exo` / `.e` / `.ex2`) uses the same two scripts. Details: [python.md](python.md).

## Skin

```bash
./build/skin hex_cube hex_skin
```

Writes a surface mesh and `hex_skin/parent_sideset`.
