#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace smesh;

// FDA benchmark nozzle (mm): 12 mm inlet, 20-degree cone, 4 mm throat, optional
// sudden expansion back to 12 mm at x = 0. expansion < 0: no expansion (pipe).
// n_axial: one count for every segment, or comma-separated per-segment grading.
static int parse_n_axial(const char *arg, const ptrdiff_t n_segments, std::vector<ptrdiff_t> *out) {
    out->clear();
    std::stringstream ss(arg);
    std::string       tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) {
            return SMESH_FAILURE;
        }
        const ptrdiff_t n = static_cast<ptrdiff_t>(std::atol(tok.c_str()));
        if (n <= 0) {
            return SMESH_FAILURE;
        }
        out->push_back(n);
    }
    if (out->size() == 1) {
        const ptrdiff_t n = (*out)[0];
        out->assign(static_cast<size_t>(n_segments), n);
    }
    return out->size() == static_cast<size_t>(n_segments) ? SMESH_SUCCESS : SMESH_FAILURE;
}

int main(int argc, char **argv) {
    SMESH_TRACE_SCOPE("nozzle.exe");
    auto ctx  = initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc != 7) {
        if (!comm->rank()) {
            fprintf(stderr,
                    "Usage: %s <n_core> <n_bore> <n_outer> <expansion> <n_axial> "
                    "<output_folder>\n"
                    "  FDA HEX8 nozzle (mm). n_core even and >= 2.\n"
                    "  expansion: station index of the sudden expansion, or < 0 for a pipe.\n"
                    "  n_axial: cells per segment (one integer = uniform; 6,5,10,30 = graded).\n",
                    argv[0]);
        }
        return SMESH_FAILURE;
    }

    const std::vector<geom_t>    x_breaks    = {-100, -62.685f, -40, 0, 120};
    const std::vector<geom_t>    bore_radius = {6, 6, 2, 2, 2};
    const geom_t                 r_exp       = 6;
    const ptrdiff_t              n_segments  = 4;
    const ptrdiff_t              n_core      = std::atoi(argv[1]);
    const ptrdiff_t              n_bore      = std::atoi(argv[2]);
    const ptrdiff_t              n_outer     = std::atoi(argv[3]);
    const ptrdiff_t              expansion   = std::atoi(argv[4]);
    std::vector<ptrdiff_t>       n_axial;
    if (parse_n_axial(argv[5], n_segments, &n_axial) != SMESH_SUCCESS) {
        if (!comm->rank()) {
            fprintf(stderr, "nozzle: n_axial must be one integer or %td comma-separated counts\n",
                    n_segments);
        }
        return SMESH_FAILURE;
    }
    const geom_t core_fraction = Env::read<geom_t>("SMESH_NOZZLE_CORE_FRACTION", static_cast<geom_t>(0.5));

    auto mesh = Mesh::create_hex8_nozzle(comm,
                                         x_breaks,
                                         bore_radius,
                                         n_axial,
                                         expansion,
                                         r_exp,
                                         n_core,
                                         n_bore,
                                         n_outer,
                                         core_fraction);
    if (!mesh) {
        return SMESH_FAILURE;
    }
    return mesh->write(Path(argv[6]));
}
