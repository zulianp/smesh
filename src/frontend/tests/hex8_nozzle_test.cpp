// create_hex8_nozzle: orientation, volume, orphans and the topology of the expansion.
//
// The geometry is the FDA benchmark nozzle in millimetres: a 12 mm inlet, a 20-degree cone
// of 22.685 mm, a 40 mm long 4 mm throat and a sudden expansion back to 12 mm at x = 0.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "smesh_mesh.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
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
    double hex_volume_nodes(const geom_t *const *p, const idx_t node[8]) {
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
                    for (int k = 0; k < 3; ++k) J[c][k] += (double)p[c][node[a]] * dN[k];
            }
            v += J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                 J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                 J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        }
        return v;
    }

    double hex_volume(const geom_t *const *p, const idx_t *const *el, const ptrdiff_t e) {
        idx_t node[8];
        for (int a = 0; a < 8; ++a) node[a] = el[a][e];
        return hex_volume_nodes(p, node);
    }

    // Volume and non-positive micro-corner count of a semi-structured hex mesh, micro cell by
    // micro cell, in the standard corner order.
    double ss_volume(const std::shared_ptr<Mesh> &ss, const int L, ptrdiff_t &bad_corners) {
        const auto *const *p  = ss->points()->data();
        const auto *const *el = ss->elements(0)->data();
        double              v = 0;
        bad_corners           = 0;
        static const int nb[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                     {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};
        for (ptrdiff_t e = 0; e < ss->n_elements(0); ++e)
            for (int zi = 0; zi < L; ++zi)
                for (int yi = 0; yi < L; ++yi)
                    for (int xi = 0; xi < L; ++xi) {
                        const idx_t node[8] = {el[sshex8_lidx(L, xi, yi, zi)][e],
                                               el[sshex8_lidx(L, xi + 1, yi, zi)][e],
                                               el[sshex8_lidx(L, xi + 1, yi + 1, zi)][e],
                                               el[sshex8_lidx(L, xi, yi + 1, zi)][e],
                                               el[sshex8_lidx(L, xi, yi, zi + 1)][e],
                                               el[sshex8_lidx(L, xi + 1, yi, zi + 1)][e],
                                               el[sshex8_lidx(L, xi + 1, yi + 1, zi + 1)][e],
                                               el[sshex8_lidx(L, xi, yi + 1, zi + 1)][e]};
                        v += hex_volume_nodes(p, node);
                        for (int a = 0; a < 8; ++a) {
                            double d[3][3];
                            for (int k = 0; k < 3; ++k)
                                for (int c = 0; c < 3; ++c)
                                    d[k][c] = (double)p[c][node[nb[a][k]]] - (double)p[c][node[a]];
                            const double det = d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
                                               d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
                                               d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
                            if (!(det > 0)) ++bad_corners;
                        }
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

    // The grading arguments default to 0, so every call that predates grading measures exactly
    // the mesh it always did.
    bool measure(const ptrdiff_t n_core, const ptrdiff_t n_bore, const ptrdiff_t n_outer, Measured &out,
                 const geom_t radial_grading = 0, const geom_t axial_grading = 0) {
        auto mesh = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n_core,
                                             n_bore, n_outer, 0.5, radial_grading, axial_grading);
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

int test_nozzle_semistructured_warp() {
    // A coarse nozzle refined to level L must, once warped, BE the nozzle at L times the
    // resolution -- same node set, same volume -- and not the polygon its macro corners span.
    const int       L = 4;
    const ptrdiff_t n = 2, nb = 1, no = 1;
    auto coarse = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n, nb, no);
    SMESH_TEST_ASSERT(coarse != nullptr);
    auto ss = to_semistructured(L, coarse, true, false);
    SMESH_TEST_ASSERT(ss != nullptr);

    std::vector<ptrdiff_t> na_fine;
    for (const auto a : NA) na_fine.push_back(a * L);
    auto fine = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, na_fine, EXPANSION, R_EXP, n * L, nb * L,
                                         no * L);
    SMESH_TEST_ASSERT(fine != nullptr);
    double v_fine = 0;
    {
        const auto *const *p  = fine->points()->data();
        const auto *const *el = fine->elements(0)->data();
        for (ptrdiff_t e = 0; e < fine->n_elements(0); ++e) v_fine += hex_volume(p, el, e);
    }

    ptrdiff_t    bad     = 0;
    const double v_chord = ss_volume(ss, L, bad);
    // On the chords the lattice misses the curved volume by the macro polygon's deficit.
    SMESH_TEST_ASSERT(std::fabs(v_chord / v_fine - 1.0) > 1e-3);

    SMESH_TEST_EQ(Mesh::warp_semistructured_hex8_nozzle(ss, X, RB, NA, EXPANSION, R_EXP, n, nb, no),
                  (int)SMESH_SUCCESS);
    const double v_warp = ss_volume(ss, L, bad);
    SMESH_TEST_EQ(bad, (ptrdiff_t)0);
    SMESH_TEST_ASSERT(std::fabs(v_warp / v_fine - 1.0) < 1e-5);
    SMESH_TEST_EQ(ss->n_nodes(), fine->n_nodes());

    // Node for node: every warped node is a node of the fine nozzle.
    const auto *const     *pf = fine->points()->data();
    const auto *const     *ps = ss->points()->data();
    std::vector<ptrdiff_t> order((size_t)fine->n_nodes());
    for (ptrdiff_t i = 0; i < fine->n_nodes(); ++i) order[(size_t)i] = i;
    std::sort(order.begin(), order.end(), [&](ptrdiff_t a, ptrdiff_t b) { return pf[0][a] < pf[0][b]; });
    double worst = 0;
    for (ptrdiff_t i = 0; i < ss->n_nodes(); ++i) {
        const double x  = ps[0][i];
        auto         lo = std::lower_bound(order.begin(), order.end(), x - 1e-3,
                                           [&](ptrdiff_t a, double v) { return pf[0][a] < v; });
        double best = 1e300;
        for (auto it = lo; it != order.end() && pf[0][*it] <= x + 1e-3; ++it) {
            const double dx = (double)pf[0][*it] - x;
            const double dy = (double)pf[1][*it] - (double)ps[1][i];
            const double dz = (double)pf[2][*it] - (double)ps[2][i];
            best            = std::min(best, dx * dx + dy * dy + dz * dz);
        }
        worst = std::max(worst, std::sqrt(best));
    }
    SMESH_TEST_ASSERT(worst < 1e-3);
    return SMESH_TEST_SUCCESS;
}

int test_nozzle_grading() {
    // Grading moves nodes; it must move nothing else.
    //
    // beta = 0 is the identity, and not approximately: nozzle_stretch returns its argument
    // before any arithmetic, so a mesh asked for with zero grading is the ungraded mesh bit for
    // bit. Asserting that here is what lets every number ever recorded on this geometry stand.
    const ptrdiff_t n = 4, nb = 3, no = 2;
    Measured        plain, zero, graded;
    SMESH_TEST_ASSERT(measure(n, nb, no, plain));
    SMESH_TEST_ASSERT(measure(n, nb, no, zero, 0, 0));
    SMESH_TEST_EQ(zero.n_elements, plain.n_elements);
    SMESH_TEST_EQ(zero.n_skin_faces, plain.n_skin_faces);
    SMESH_TEST_EQ(zero.n_axis_nodes, plain.n_axis_nodes);
    SMESH_TEST_ASSERT(zero.volume_rel_err == plain.volume_rel_err);

    // With grading on the topology is untouched -- same elements, same skin, same axis -- and
    // the mesh stays valid: the stretch is monotone, so no cell folds and no node is orphaned.
    SMESH_TEST_ASSERT(measure(n, nb, no, graded, 2.0, 1.5));
    SMESH_TEST_EQ(graded.n_elements, plain.n_elements);
    SMESH_TEST_EQ(graded.n_skin_faces, plain.n_skin_faces);
    SMESH_TEST_EQ(graded.n_axis_nodes, plain.n_axis_nodes);
    SMESH_TEST_EQ(graded.n_bad_corners, (ptrdiff_t)0);
    SMESH_TEST_EQ(graded.n_orphans, (ptrdiff_t)0);

    // The volume is the SAME volume -- measured, not assumed.
    //
    // Grading redistributes nodes and changes nothing about how much space the mesh occupies.
    // The cross-section's polygonal area telescopes: splitting the annulus at different radii
    // leaves the area between the innermost and outermost rings unchanged, so the 4*n_core-gon's
    // deficit is fixed by the ANGULAR resolution alone. Axially the segments are frusta,
    // integrated exactly at any plane spacing. Measured here: -2.550465e-02 relative, graded and
    // ungraded alike, ratio 1.0000. The tolerance is relative rather than exact only because the
    // volumes are summed in a different order once the nodes move.
    std::printf("grading: volume_rel_err plain %.6e  graded %.6e  ratio %.4f\n", plain.volume_rel_err,
                graded.volume_rel_err, graded.volume_rel_err / plain.volume_rel_err);
    SMESH_TEST_ASSERT(graded.volume_rel_err < 0);
    SMESH_TEST_ASSERT(std::fabs(graded.volume_rel_err / plain.volume_rel_err - 1.0) < 1e-6);

    // Grading must nevertheless MOVE something, or every assertion above would hold on a no-op.
    // Node positions are what it moves, so that is what this checks.
    auto m_plain  = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n, nb, no);
    auto m_graded = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n, nb, no,
                                             0.5, 2.0, 1.5);
    SMESH_TEST_ASSERT(m_plain != nullptr && m_graded != nullptr);
    SMESH_TEST_EQ(m_graded->n_nodes(), m_plain->n_nodes());
    {
        const auto *const *a     = m_plain->points()->data();
        const auto *const *b     = m_graded->points()->data();
        double             worst = 0;
        for (ptrdiff_t i = 0; i < m_plain->n_nodes(); ++i)
            for (int c = 0; c < 3; ++c)
                worst = std::max(worst, std::fabs((double)a[c][i] - (double)b[c][i]));
        SMESH_TEST_ASSERT(worst > 1e-3);
    }

    // And beta = 0 moves NOTHING: bit for bit, not to within a tolerance. This is the assertion
    // that lets every number recorded on the ungraded nozzle stand unchanged.
    {
        auto m_zero = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n, nb, no,
                                               0.5, 0, 0);
        SMESH_TEST_ASSERT(m_zero != nullptr);
        SMESH_TEST_EQ(m_zero->n_nodes(), m_plain->n_nodes());
        const auto *const *a         = m_plain->points()->data();
        const auto *const *b         = m_zero->points()->data();
        ptrdiff_t          differing = 0;
        for (ptrdiff_t i = 0; i < m_plain->n_nodes(); ++i)
            for (int c = 0; c < 3; ++c)
                if (a[c][i] != b[c][i]) ++differing;
        SMESH_TEST_EQ(differing, (ptrdiff_t)0);
    }
    return SMESH_TEST_SUCCESS;
}

int test_nozzle_graded_warp() {
    // The warp identity, with grading on.
    //
    // This is why NozzlePlanes::at() evaluates the axial map at a CONTINUOUS plane coordinate.
    // Grading only the planes would leave the micro nodes between them uniformly spaced --
    // piecewise-linear grading, with a kink at every macro cell boundary -- and a warped lattice
    // would then NOT be the fine graded nozzle. This asserts it is.
    const int       L = 4;
    const ptrdiff_t n = 2, nb = 1, no = 1;
    const geom_t    br = 2.0, ba = 1.5;

    auto coarse = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, NA, EXPANSION, R_EXP, n, nb, no,
                                           0.5, br, ba);
    SMESH_TEST_ASSERT(coarse != nullptr);
    auto ss = to_semistructured(L, coarse, true, false);
    SMESH_TEST_ASSERT(ss != nullptr);

    std::vector<ptrdiff_t> na_fine;
    for (const auto a : NA) na_fine.push_back(a * L);
    auto fine = Mesh::create_hex8_nozzle(Communicator::self(), X, RB, na_fine, EXPANSION, R_EXP, n * L,
                                         nb * L, no * L, 0.5, br, ba);
    SMESH_TEST_ASSERT(fine != nullptr);
    double v_fine = 0;
    {
        const auto *const *p  = fine->points()->data();
        const auto *const *el = fine->elements(0)->data();
        for (ptrdiff_t e = 0; e < fine->n_elements(0); ++e) v_fine += hex_volume(p, el, e);
    }

    ptrdiff_t bad = 0;
    SMESH_TEST_EQ(Mesh::warp_semistructured_hex8_nozzle(ss, X, RB, NA, EXPANSION, R_EXP, n, nb, no, 0.5, br, ba),
                  (int)SMESH_SUCCESS);
    const double v_warp = ss_volume(ss, L, bad);
    SMESH_TEST_EQ(bad, (ptrdiff_t)0);
    SMESH_TEST_ASSERT(std::fabs(v_warp / v_fine - 1.0) < 1e-5);
    SMESH_TEST_EQ(ss->n_nodes(), fine->n_nodes());
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_nozzle_topology);
    SMESH_RUN_TEST(test_nozzle_volume_converges);
    SMESH_RUN_TEST(test_nozzle_semistructured_warp);
    SMESH_RUN_TEST(test_nozzle_grading);
    SMESH_RUN_TEST(test_nozzle_graded_warp);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
