# Folder format

On-disk meshes are directories, not a single binary. C++ `Mesh::write` / `Mesh::read` and the Python converters share this layout.

## Single block

```
mesh/
  meta.yaml
  i0.int32 …          # one file per local node (SoA connectivity)
  x.float32           # or float64, matching geom_t
  y.float32
  z.float32           # omitted when spatial_dimension is 2
```

`meta.yaml` names the files, dtypes, element type, and counts. Integer and float suffixes follow the array dtype (`int32`, `int64`, `float32`, `float64`, …).

## Multiple blocks

Shared coordinates at the mesh root. Connectivity per block:

```
mesh/
  meta.yaml
  x.float32 y.float32 z.float32
  blocks/
    fluid/i0.int32 …
    solid/i0.int32 …
```

## Sets

```
sidesets/<name>/{parent,lfi,meta.yaml}
sidesets/<name>/<block_id>/…    # when one name spans blocks
nodesets/<name>/{nodes,meta.yaml}
```

`parent` is the 0-based local element index in that block. `lfi` is 0-based.

## Converters

`python/smesh/db_to_raw.py` and `raw_to_db.py` map this folder to VTK/Exodus and back. They do not change topology (no refine, no promote).
