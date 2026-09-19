# C++

Headers live under `src/` (the `smesh` target adds those include paths). C++17.

```cpp
#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"

int main(int argc, char **argv) {
    auto ctx = smesh::initialize_serial(argc, argv);
    auto mesh = smesh::Mesh::create_hex8_cube(ctx->communicator(), 8, 8, 8);
    auto refined = smesh::refine(mesh, 1);
    return refined->write(smesh::Path("hex_refined"));
}
```

Link `smesh` (exported as `smesh::smesh`). After install:

```cmake
find_package(smesh REQUIRED)
target_link_libraries(app PRIVATE smesh::smesh)
```

Useful entry points in `smesh_mesh.hpp`:

- Create: `Mesh::create_hex8_cube`, `create_tet4_cube`, `create_from_file`
- Transform: `refine`, `promote_to`, `convert_to`
- IO: `write`, `write_with_xdmf` (adds `mesh.xdmf` in the folder)

For MPI, call `smesh::initialize` instead of `initialize_serial` and pass the communicator into factories. `create_from_file` reads a serial folder on every rank and distributes it.
