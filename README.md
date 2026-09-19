# SMESH

C++ mesh library for serial and MPI workflows. Meshes are stored as folders of SoA arrays plus `meta.yaml`. Python converters translate that layout to and from VTK and Exodus; they do not refine or promote.

## Build

CMake 3.20+, C++17. A first build without MPI:

```bash
cmake -S . -B build \
  -DSMESH_ENABLE_MPI=OFF \
  -DSMESH_ENABLE_DEV_MODE=OFF
cmake --build build -j
```

## First mesh

```bash
./build/cube HEX8 8 8 8 0 0 0 1 1 1 hex_cube
```

That writes a folder mesh. Next steps: [Getting started](docs/getting-started.md).

## Documentation

- [Install](docs/install.md)
- [Getting started](docs/getting-started.md)
- [Concepts](docs/concepts.md)
- [Features and element types](docs/features.md)
- [Folder format](docs/format.md)
- [Command-line tools](docs/cli.md)
- [C++](docs/cpp.md)
- [Python converters](docs/python.md)
- [Performance](docs/PERFORMANCE.md) (optional, HPC)

License: [BSD-3-Clause](LICENSE).
