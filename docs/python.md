# Python converters

Scripts in `python/smesh/` translate a folder mesh to VTK or Exodus and back. They do not refine, promote, or convert element families.

There is no installable `smesh` package and no C++ Python bindings in this tree. Run the scripts from `python/smesh/` so their `common/` imports resolve.

## Dependencies

```bash
pip install -r python/requirements.txt
```

That is `numpy`, `meshio`, `netCDF4`, and `pyyaml`. Extra scientific packages used elsewhere in the Python tree are listed in `python/requirements-extra.txt` and are not required here.

## Folder → VTK / Exodus

```bash
python3 python/smesh/raw_to_db.py hex_cube hex_cube.vtk
python3 python/smesh/raw_to_db.py hex_cube hex_cube.exo
```

`raw_to_exodusII.py` is the Exodus writer; `raw_to_db.py` dispatches to it for `.exo` / `.e` / `.ex2`.

## VTK / Exodus → folder

```bash
python3 python/smesh/db_to_raw.py hex_cube.vtk hex_cube
python3 python/smesh/db_to_raw.py hex_cube.exo hex_cube
```

`db_to_raw.py` uses `exodusII_to_raw.py` for Exodus paths.

## Tests

`python/smesh/tests/` are script-level pytest cases (synthetic meshes). Run them from `python/smesh/`:

```bash
cd python/smesh && python3 -m pytest tests
```
