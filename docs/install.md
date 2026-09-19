# Install

## Requirements

- CMake 3.20 or newer
- A C++17 compiler
- Git (CMake fetches [matrixio](https://github.com/zulianp/matrix.io) at configure time)

Optional:

- MPI, plus the `external/mpi-sort` submodule
- CUDA Toolkit (`-DSMESH_ENABLE_CUDA=ON`)
- OpenMP (`-DSMESH_ENABLE_OPENMP=ON`)
- Doxygen (`-DSMESH_ENABLE_INSTALL_DOCS=ON`, CMake target `docs` → `api/html`)
- Python 3 with `numpy`, `meshio`, `netCDF4`, `pyyaml` for VTK/Exodus converters ([python.md](python.md))

## Configure and build

Serial (recommended first build):

```bash
cmake -S . -B build \
  -DSMESH_ENABLE_MPI=OFF \
  -DSMESH_ENABLE_DEV_MODE=OFF
cmake --build build -j
```

MPI:

```bash
git submodule update --init --recursive external/mpi-sort
cmake -S . -B build \
  -DSMESH_ENABLE_MPI=ON \
  -DSMESH_ENABLE_DEV_MODE=OFF
cmake --build build -j
```

CUDA (device buffers and SS restrict/prolong kernels):

```bash
cmake -S . -B build \
  -DSMESH_ENABLE_MPI=OFF \
  -DSMESH_ENABLE_DEV_MODE=OFF \
  -DSMESH_ENABLE_CUDA=ON
cmake --build build -j
```

Default `CMAKE_CUDA_ARCHITECTURES` is `90`. Override it for the GPU you have.

## Tests

Tests are on by default at the top level (`SMESH_ENABLE_TESTING=ON`):

```bash
ctest --test-dir build --output-on-failure
```

## CMake options

Stock defaults are aimed at in-tree development, not a first clone:

| Option | Default | Notes |
|--------|---------|--------|
| `SMESH_ENABLE_MPI` | `ON` | Needs MPI and `external/mpi-sort` |
| `SMESH_ENABLE_DEV_MODE` | `ON` | `-Wall -Wextra -Werror` |
| `SMESH_ENABLE_TRACE` | `ON` | Writes `smesh.trace.csv` |
| `SMESH_ENABLE_TESTING` | `ON` (top-level) | `ctest` |
| `SMESH_ENABLE_OPENMP` | `OFF` | |
| `SMESH_ENABLE_CUDA` | `OFF` | CUDA Toolkit; default arch `90` |
| `SMESH_ENABLE_CUDA_LINEINFO` | `OFF` | Device-line profiling |
| `CMAKE_BUILD_TYPE` | `Release` if unset | |

For a first build, pass `-DSMESH_ENABLE_MPI=OFF` and `-DSMESH_ENABLE_DEV_MODE=OFF`.

C++ Python bindings (`SMESH_ENABLE_PYTHON`) are not enabled in this tree. Use the scripts in `python/smesh/` instead.

Markdown user docs live in this `docs/` folder. The CMake target `docs` is Doxygen API HTML, not these pages.
