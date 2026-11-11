// Lossy coordinates bounded by the element height (compression mode h): heights against analytic
// values, then strongly graded HEX8 and TET4 cubes written with h:E through Mesh::write and read
// back. Every node must move by at most E * h_n per component, and every simplex (for hexes:
// every corner tet) must keep its orientation with a relative volume change of at most 4 sqrt(3) E
// (|dV| / V <= sum over the 4 nodes of |dx| / height <= 4 sqrt(3) E).

#include "smesh_communicator.hpp"
#include "smesh_compression_options.hpp"
#include "smesh_mesh.hpp"
#include "smesh_node_height.hpp"
#include "smesh_smz_format.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_test.hpp"

#include <cmath>
#include <filesystem>
#include <vector>

using namespace smesh;
using namespace smesh::compression;
namespace fs = std::filesystem;

static const std::string kDir = "smz_node_height_test_tmp";

static int test_analytic_heights() {
    const double hex[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    SMESH_TEST_APPROXEQ(element_min_height(HEX8, hex), 1.0 / std::sqrt(3.0), 1e-14);
    const double a      = 2.0;  // regular tet of edge a: height a sqrt(2/3)
    const double tet[4][3] = {{0, 0, 0}, {a, 0, 0}, {a / 2, a * std::sqrt(3.0) / 2, 0},
                              {a / 2, a * std::sqrt(3.0) / 6, a * std::sqrt(2.0 / 3.0)}};
    SMESH_TEST_APPROXEQ(element_min_height(TET4, tet), a * std::sqrt(2.0 / 3.0), 1e-14);
    const double flat[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.3, 0.3, 1e-6}};  // sliver: tiny height
    SMESH_TEST_ASSERT(element_min_height(TET4, flat) < 2e-6);
    const double tri[3][3] = {{0, 0, 0}, {4, 0, 0}, {0, 3, 0}};  // 2A / longest edge = 12 / 5
    SMESH_TEST_APPROXEQ(element_min_height(TRI3, tri), 12.0 / 5.0, 1e-14);
    SMESH_TEST_ASSERT(!node_height_supported(MACRO_TET4));
    return SMESH_TEST_SUCCESS;
}

static double tet_vol6(const double *a, const double *b, const double *c, const double *d) {
    double u[3], v[3], w[3];
    for (int k = 0; k < 3; k++) u[k] = b[k] - a[k], v[k] = c[k] - a[k], w[k] = d[k] - a[k];
    return u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}

static const int kHexCorner[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5}, {2, 3, 1, 6}, {3, 0, 2, 7},
                                     {4, 5, 7, 0}, {5, 6, 4, 1}, {6, 7, 5, 2}, {7, 4, 6, 3}};

// Graded and sheared: cells near the origin are (1/n)^4 of the domain.
static void grade(const std::shared_ptr<Mesh> &mesh) {
    geom_t **p = mesh->points()->data();
    for (ptrdiff_t i = 0; i < mesh->n_nodes(); i++) {
        const double x = p[0][i], y = p[1][i], z = p[2][i];
        p[0][i]        = geom_t(std::pow(x, 4.0));
        p[1][i]        = geom_t(y * y + 0.3 * x);
        p[2][i]        = geom_t(z + 0.1 * std::sin(3 * y));
    }
}

static int check_mesh(const std::string &name, std::shared_ptr<Mesh> mesh, double eb) {
    auto       comm = Communicator::self();
    const Path out(kDir + "/" + name);
    CompressionOptions c;
    char               spec[64];
    snprintf(spec, sizeof(spec), "h:%g,chunk:999,verify", eb);
    SMESH_TEST_EQ(CompressionOptions::parse(spec, &c), SMESH_SUCCESS);
    MeshWriteOptions o = MeshWriteOptions::all(CompressionOptions::lossless());
    o.coordinates      = c;
    SMESH_TEST_EQ(mesh->write(out, o), SMESH_SUCCESS);
    auto back = Mesh::create_from_file(comm, out);
    SMESH_TEST_ASSERT(back != nullptr);
    SMESH_TEST_EQ(back->n_nodes(), mesh->n_nodes());

    const auto           block = mesh->blocks()[0];
    const enum ElemType  et    = block->element_type();
    const ptrdiff_t      ne    = block->elements()->extent(1);
    idx_t              **el    = block->elements()->data();
    geom_t             **p     = mesh->points()->data();
    geom_t             **q     = back->points()->data();
    std::vector<double>  h(size_t(mesh->n_nodes()), INFINITY);
    SMESH_TEST_EQ((accumulate_node_heights<idx_t, geom_t>(et, ne, el, 3, p, h.data())), SMESH_SUCCESS);

    double worst = 0, hmin = INFINITY, hmax = 0;
    for (ptrdiff_t i = 0; i < mesh->n_nodes(); i++) {
        hmin = std::min(hmin, h[size_t(i)]), hmax = std::max(hmax, h[size_t(i)]);
        for (int d = 0; d < 3; d++) {
            const double err = std::fabs(double(q[d][i]) - double(p[d][i]));
            SMESH_TEST_ASSERT(err <= eb * h[size_t(i)]);
            worst = std::max(worst, err / (eb * h[size_t(i)]));
        }
    }
    // Orientation and volume of every simplex / corner tet.
    double     worst_dv = 0;
    const bool hex      = et == HEX8;
    const int  ntet     = hex ? 8 : 1;
    for (ptrdiff_t e = 0; e < ne; e++) {
        for (int t = 0; t < ntet; t++) {
            double a[4][3], b[4][3];
            for (int k = 0; k < 4; k++) {
                const idx_t n = el[hex ? kHexCorner[t][k] : k][e];
                for (int d = 0; d < 3; d++) a[k][d] = p[d][n], b[k][d] = q[d][n];
            }
            const double v0 = tet_vol6(a[0], a[1], a[2], a[3]), v1 = tet_vol6(b[0], b[1], b[2], b[3]);
            SMESH_TEST_ASSERT(v0 * v1 > 0);
            const double dv = std::fabs(v1 - v0) / std::fabs(v0);
            SMESH_TEST_ASSERT(dv <= 4 * std::sqrt(3.0) * eb * 1.01);
            worst_dv = std::max(worst_dv, dv / eb);
        }
    }
    SmzInfo info;
    SMESH_TEST_EQ(smz_read_info(Path(out.to_string() + "/x.float32.smz"), &info), SMESH_SUCCESS);
    const double raw = 3.0 * double(mesh->n_nodes()) * sizeof(geom_t);
    double       sz  = 0;
    for (const char *ax : {"x", "y", "z"}) sz += double(fs::file_size(out.to_string() + "/" + ax + ".float32.smz"));
    printf("%s eb=%g: h in [%.2e, %.2e], max err %.3f of the bound, max |dV|/V %.2f eb, coordinates %.2fx\n",
           name.c_str(), eb, hmin, hmax, worst, worst_dv, raw / sz);
    return SMESH_TEST_SUCCESS;
}

static int test_graded_meshes() {
    auto comm = Communicator::self();
    for (double eb : {1e-2, 1e-3, 1e-4}) {
        auto hex = Mesh::create_hex8_cube(comm, 24, 20, 16);
        grade(hex);
        if (check_mesh("hex8", hex, eb) != SMESH_TEST_SUCCESS) return SMESH_TEST_FAILURE;
        auto tet = Mesh::create_tet4_cube(comm, 12, 10, 8);
        grade(tet);
        if (check_mesh("tet4", tet, eb) != SMESH_TEST_SUCCESS) return SMESH_TEST_FAILURE;
    }
    return SMESH_TEST_SUCCESS;
}

// Element types without a height: the coordinates are written losslessly.
static int test_unsupported_falls_back() {
    SMESH_TEST_ASSERT(!node_height_supported(PROTEUS_HEX8));
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    fs::remove_all(kDir);
    fs::create_directories(kDir);
    SMESH_RUN_TEST(test_analytic_heights);
    SMESH_RUN_TEST(test_graded_meshes);
    SMESH_RUN_TEST(test_unsupported_falls_back);
    if (!err) fs::remove_all(kDir);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
