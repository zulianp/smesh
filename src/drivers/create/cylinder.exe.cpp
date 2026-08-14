#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_sideset.hpp"
#include "smesh_tracer.hpp"

#include <cmath>
#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("cylinder.exe");

    auto ctx = initialize_serial(argc, argv);

    if (argc != 8) {
        fprintf(stderr,
                "Usage: %s <radius> <zmin> <zmax> <nr> <ntheta> <nz> "
                "<output_folder>\n"
                "  HEX-dominant solid cylinder (axis +z). HEX8 core + polar HEX8 "
                "annulus; WEDGE6 in square-to-circle fillets.\n"
                "  ntheta >= 8 and ntheta %% 4 == 0.\n",
                argv[0]);
        return SMESH_FAILURE;
    }

    const geom_t    radius = static_cast<geom_t>(std::atof(argv[1]));
    const geom_t    zmin   = static_cast<geom_t>(std::atof(argv[2]));
    const geom_t    zmax   = static_cast<geom_t>(std::atof(argv[3]));
    const ptrdiff_t nr     = std::atoi(argv[4]);
    const ptrdiff_t ntheta = std::atoi(argv[5]);
    const ptrdiff_t nz     = std::atoi(argv[6]);
    const Path      output_folder = Path(argv[7]);
    const geom_t    height        = zmax - zmin;

    auto mesh = Mesh::create_hex_dominant_cylinder(
            ctx->communicator(), radius, height, nr, ntheta, nz, zmin);
    if (!mesh) {
        return SMESH_FAILURE;
    }

    mesh->write(output_folder);

    if (Env::read<bool>("SMESH_CREATE_SIDESETS", true)) {
        create_directory(output_folder / "sidesets");

        for (const auto &sideset : Sideset::create_from_plane(mesh, 0, 0, 1, zmin)) {
            sideset->write(output_folder / "sidesets/bottom");
        }
        for (const auto &sideset : Sideset::create_from_plane(mesh, 0, 0, 1, zmax)) {
            sideset->write(output_folder / "sidesets/top");
        }

        const geom_t r2_min = radius * radius * static_cast<geom_t>(0.99 * 0.99);
        int          wall_i = 0;
        for (const auto &sideset : Sideset::create_from_selector(
                     mesh, [r2_min](const geom_t x, const geom_t y, const geom_t) {
                         return x * x + y * y >= r2_min;
                     })) {
            if (wall_i == 0) {
                sideset->write(output_folder / "sidesets/wall");
            } else {
                sideset->write(output_folder / "sidesets/wall_extra");
            }
            ++wall_i;
        }
    }

    return SMESH_SUCCESS;
}
