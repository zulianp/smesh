// create_hex8_nozzle: orientation, volume, orphans and the topology of the expansion.
//
// The geometry is the FDA benchmark nozzle in millimetres: a 12 mm inlet, a 20-degree cone
// of 22.685 mm, a 40 mm long 4 mm throat and a sudden expansion back to 12 mm at x = 0.

#include <cmath>
#include <vector>

#include "smesh_mesh.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

using namespace smesh;

namespace {

    const std::vector<geom_t>    X  = {-100, -62.685f, -40, 0, 120};
    const std::vector<geom_t>    RB = {6, 6, 2, 2, 2};
    const std::vector<ptrdiff_t> NA = {6, 5, 10, 30};
    const ptrdiff_t              EXPANSION = 3;
    const geom_t                 R_EXP     = 6;

    // Exact volume of the smooth nozzle: frusta upstream, a full pipe downstream.
    double exact_volume() {
        double v = 0;
        for (size_t s = 0; s + 1 < X.size(); ++s) {
            const double L = (double)X[s + 1] - (double)X[s];
            if ((ptrdiff_t)s >= EXPANSION) {
                v += M_PI * (double)R_EXP * (double)R_EXP * L;
            } else {
                const double a = RB[s], b = RB[s + 1];
                v += M_PI * L * (a * a + a * b + b * b) / 3.0;
            }
        }
        return v;
    }

    // Trilinear Jacobian determinant at a corner, from its three edges. Positive at all
    // eight corners is the orientation test a consumer's quadrature relies on.
    double corner_det(const geom_t *const *p, const idx_t *const *el, const ptrdiff_t e, const int a) {
        // Neighbours of each corner along the local xi, eta, zeta directions, ordered so that
        // a right-handed element gives a positive triple product at every corner.
        static const int nb[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                     {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};
        const idx_t o = el[a][e];
        double      d[3][3];
        for (int k = 0; k < 3; ++k) {
            const idx_t q = el[nb[a][k]][e];
            for (int c = 0; c < 3; ++c) d[k][c] = (double)p[c][q] - (double)p[c][o];
        }
        return d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
               d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
               d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
    }

    // Volume by the 2x2x2 Gauss rule on the trilinear map, which is what an isoparametric
    // consumer integrates.
    double hex_volume(const geom_t *const *p, const idx_t *const *el, const ptrdiff_t e) {
        const double g = 1.0 / std::sqrt(3.0);
        const int    s[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                                {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
        double       v = 0;
        for (int qa = 0; qa < 8; ++qa) {
            const double xi = s[qa][0] * g, et = s[qa][1] * g, ze = s[qa][2] * g;
            double       J[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            for (int a = 0; a < 8; ++a) {
                const double dN[3] = {0.125 * s[a][0] * (1 + s[a][1] * et) * (1 + s[a][2] * ze),
                                      0.125 * s[a][1] * (1 + s[a][0] * xi) * (1 + s[a][2] * ze),
                                      0.125 * s[a][2] * (1 + s[a][0] * xi) * (1 + s[a][1] * et)};
                for (int c = 0; c < 3; ++c)
                    for (int k = 0; k < 3; ++k) J[c][k] += (double)p[c][el[a][e]] * dN[k];
            }
            v += J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                 J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                 J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        }
        return v;
    }

    struct Measured {
        double    volume_rel_err;
        ptrdiff_t n_bad_corners;
        ptrdiff_t n_orphans;
        ptrdiff_t n_axis_nodes;
        ptrdiff_t n_elements;
        ptrdiff_t n_skin_faces;
    };

    bool measure(const ptrdiff_t n_core, const ptrdiff_t n_bore, const ptrdiff_t n_outer, Measured &out) {
        auto mesh = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n_core,
                                             n_bore, n_outer);
        if (!mesh) return false;
        const auto *const *p  = mesh->points()->data();
        const auto *const *el = mesh->elements(0)->data();
        const ptrdiff_t    ne = mesh->n_elements(0);
        const ptrdiff_t    nn = mesh->n_nodes();

        double            vol = 0;
        ptrdiff_t         bad = 0;
        std::vector<char> seen((size_t)nn, 0);
        for (ptrdiff_t e = 0; e < ne; ++e) {
            vol += hex_volume(p, el, e);
            for (int a = 0; a < 8; ++a) {
                if (!(corner_det(p, el, e, a) > 0)) ++bad;
                seen[(size_t)el[a][e]] = 1;
            }
        }
        ptrdiff_t orphans = 0, axis = 0;
        for (ptrdiff_t i = 0; i < nn; ++i) {
            if (!seen[(size_t)i]) ++orphans;
            if (p[1][i] == 0 && p[2][i] == 0) ++axis;
        }
        auto skin = skin_sideset(mesh);

        out.volume_rel_err = vol / exact_volume() - 1.0;
        out.n_bad_corners  = bad;
        out.n_orphans      = orphans;
        out.n_axis_nodes   = axis;
        out.n_elements     = ne;
        out.n_skin_faces   = skin ? skin->parent()->size() : -1;
        return true;
    }

}  // namespace

int test_nozzle_topology() {
    const ptrdiff_t n = 4, nb = 3, no = 2;
    Measured        m;
    SMESH_TEST_ASSERT(measure(n, nb, no, m));

    ptrdiff_t cells_up = 0, cells_down = 0;
    for (size_t s = 0; s < NA.size(); ++s) ((ptrdiff_t)s >= EXPANSION ? cells_down : cells_up) += NA[s];
    const ptrdiff_t perim = 4 * n;

    const ptrdiff_t zero = 0;
    SMESH_TEST_EQ(m.n_bad_corners, zero);
    SMESH_TEST_EQ(m.n_orphans, zero);
    // The core's node lines on both symmetry planes meet on the axis in every plane.
    SMESH_TEST_EQ(m.n_axis_nodes, cells_up + cells_down + 1);
    SMESH_TEST_EQ(m.n_elements, (cells_up + cells_down) * (n * n + nb * perim) + cells_down * no * perim);

    // The skin, face by face: the inlet disc, the outlet disc, the bore wall upstream, the
    // expanded wall downstream -- and the annulus at the expansion, which exists only if the
    // notch was carved rather than the outer ring being squeezed onto the bore.
    const ptrdiff_t inlet   = n * n + nb * perim;
    const ptrdiff_t outlet  = n * n + (nb + no) * perim;
    const ptrdiff_t walls   = perim * (cells_up + cells_down);
    const ptrdiff_t annulus = no * perim;
    SMESH_TEST_EQ(m.n_skin_faces, inlet + outlet + walls + annulus);
    return SMESH_TEST_SUCCESS;
}

int test_nozzle_volume_converges() {
    // A polygon of 4 n_core sides inscribed in the circle loses area like (2 pi / 4n)^2 / 6,
    // so doubling the cells across should cut the volume error by about four.
    Measured coarse, fine;
    SMESH_TEST_ASSERT(measure(8, 3, 3, coarse));
    SMESH_TEST_ASSERT(measure(16, 6, 6, fine));
    SMESH_TEST_ASSERT(coarse.volume_rel_err < 0);
    SMESH_TEST_ASSERT(std::fabs(coarse.volume_rel_err) < 1e-2);
    SMESH_TEST_ASSERT(std::fabs(fine.volume_rel_err) < 0.3 * std::fabs(coarse.volume_rel_err));
    SMESH_TEST_EQ(fine.n_bad_corners, (ptrdiff_t)0);
    return SMESH_TEST_SUCCESS;
}

// No test of the input validation: SMESH_ERROR aborts the process, so a rejected input
// cannot be observed from inside a test.

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_nozzle_topology);
    SMESH_RUN_TEST(test_nozzle_volume_converges);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
