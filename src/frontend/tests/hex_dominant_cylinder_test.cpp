#include <cmath>
#include <vector>

#include "smesh_mesh.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

using namespace smesh;

static geom_t tet_signed_volume(const geom_t *const *const p, const idx_t a, const idx_t b,
                                const idx_t c, const idx_t d) {
    const geom_t abx = p[0][b] - p[0][a];
    const geom_t aby = p[1][b] - p[1][a];
    const geom_t abz = p[2][b] - p[2][a];
    const geom_t acx = p[0][c] - p[0][a];
    const geom_t acy = p[1][c] - p[1][a];
    const geom_t acz = p[2][c] - p[2][a];
    const geom_t adx = p[0][d] - p[0][a];
    const geom_t ady = p[1][d] - p[1][a];
    const geom_t adz = p[2][d] - p[2][a];
    return (abx * (acy * adz - acz * ady) + aby * (acz * adx - acx * adz) +
            abz * (acx * ady - acy * adx)) /
           static_cast<geom_t>(6);
}

static geom_t hex_signed_volume(const geom_t *const *const p, idx_t *const *const el,
                                const ptrdiff_t e) {
    const idx_t i[8] = {el[0][e], el[1][e], el[2][e], el[3][e],
                        el[4][e], el[5][e], el[6][e], el[7][e]};
    const idx_t tets[6][4] = {{i[0], i[1], i[3], i[7]}, {i[0], i[1], i[7], i[5]},
                              {i[0], i[4], i[5], i[7]}, {i[1], i[2], i[3], i[6]},
                              {i[1], i[3], i[7], i[6]}, {i[1], i[5], i[6], i[7]}};
    geom_t v = 0;
    for (int t = 0; t < 6; ++t) {
        v += tet_signed_volume(p, tets[t][0], tets[t][1], tets[t][2], tets[t][3]);
    }
    return v;
}

static geom_t wedge_signed_volume(const geom_t *const *const p, idx_t *const *const el,
                                  const ptrdiff_t e) {
    const idx_t i0 = el[0][e];
    const idx_t i1 = el[1][e];
    const idx_t i2 = el[2][e];
    const idx_t i3 = el[3][e];
    const idx_t i4 = el[4][e];
    const idx_t i5 = el[5][e];
    return tet_signed_volume(p, i0, i1, i2, i3) + tet_signed_volume(p, i1, i4, i5, i3) +
           tet_signed_volume(p, i2, i1, i5, i3);
}

static int check_half_face_tables(Mesh &mesh, const int require_inter_block) {
    int n_inter_block = 0;
    for (size_t b = 0; b < mesh.n_blocks(); ++b) {
        const block_idx_t   bid  = static_cast<block_idx_t>(b);
        const enum ElemType type = mesh.element_type(bid);
        const int           ns   = elem_num_sides(type);
        const ptrdiff_t     n_e  = mesh.n_elements(bid);
        auto                hft  = mesh.half_face_table(bid);
        auto                hnbb = mesh.half_face_neighbor_block(bid);
        SMESH_TEST_ASSERT(hft != nullptr);
        SMESH_TEST_ASSERT(hnbb != nullptr);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(hft->size()), n_e * ns);
        SMESH_TEST_EQ(static_cast<ptrdiff_t>(hnbb->size()), n_e * ns);

        for (ptrdiff_t e = 0; e < n_e; ++e) {
            for (int s = 0; s < ns; ++s) {
                const element_idx_t neighbor = hft->data()[e * ns + s];
                const block_idx_t   nb       = hnbb->data()[e * ns + s];
                if (neighbor != invalid_idx<element_idx_t>()) {
                    SMESH_TEST_ASSERT(nb >= 0);
                    SMESH_TEST_ASSERT(static_cast<size_t>(nb) < mesh.n_blocks());
                    SMESH_TEST_ASSERT(neighbor < mesh.n_elements(nb));
                    if (nb != bid) {
                        n_inter_block++;
                    }
                }
            }
        }
    }
    if (require_inter_block) {
        SMESH_TEST_ASSERT(n_inter_block > 0);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_hex_dominant_cylinder_mesh(const geom_t radius, const geom_t height,
                                           const ptrdiff_t nr, const ptrdiff_t ntheta,
                                           const ptrdiff_t nz, const geom_t zmin) {
    auto mesh = Mesh::create_hex_dominant_cylinder(
            Communicator::self(), radius, height, nr, ntheta, nz, zmin);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->n_blocks() >= 1);
    SMESH_TEST_EQ(mesh->element_type(0), HEX8);
    SMESH_TEST_ASSERT(mesh->n_elements(0) > 0);

    const ptrdiff_t n_hex = mesh->n_elements(0);
    ptrdiff_t       n_wedge = 0;
    if (mesh->n_blocks() > 1) {
        SMESH_TEST_EQ(mesh->element_type(1), WEDGE6);
        n_wedge = mesh->n_elements(1);
        SMESH_TEST_ASSERT(n_wedge > 0);
    }
    SMESH_TEST_ASSERT(n_hex > n_wedge);

    auto              pts = mesh->points()->data();
    const geom_t      rmax = radius * static_cast<geom_t>(1.0001);
    const geom_t      zmax = zmin + height;
    for (ptrdiff_t n = 0; n < mesh->n_nodes(); ++n) {
        const geom_t r2 = pts[0][n] * pts[0][n] + pts[1][n] * pts[1][n];
        SMESH_TEST_ASSERT(r2 <= rmax * rmax);
        SMESH_TEST_ASSERT(pts[2][n] >= zmin - static_cast<geom_t>(1e-6));
        SMESH_TEST_ASSERT(pts[2][n] <= zmax + static_cast<geom_t>(1e-6));
    }

    geom_t vol = 0;
    {
        auto hex = mesh->elements(0)->data();
        for (ptrdiff_t e = 0; e < n_hex; ++e) {
            const geom_t v = hex_signed_volume(pts, hex, e);
            SMESH_TEST_ASSERT(v > static_cast<geom_t>(0));
            vol += v;
        }
    }
    if (n_wedge > 0) {
        auto wedge = mesh->elements(1)->data();
        for (ptrdiff_t e = 0; e < n_wedge; ++e) {
            const geom_t v = wedge_signed_volume(pts, wedge, e);
            SMESH_TEST_ASSERT(v > static_cast<geom_t>(0));
            vol += v;
        }
    }

    const geom_t prism = static_cast<geom_t>(
            0.5 * static_cast<double>(ntheta) * static_cast<double>(radius) *
            static_cast<double>(radius) * std::sin(2.0 * M_PI / static_cast<double>(ntheta)) *
            static_cast<double>(height));
    const geom_t rel = std::fabs(static_cast<double>(vol - prism)) /
                       static_cast<geom_t>(std::fabs(static_cast<double>(prism)));
    SMESH_TEST_ASSERT(rel < static_cast<geom_t>(1e-4));

    {
        const geom_t ztol = static_cast<geom_t>(1e-6);
        int n_corners = 0;
        int n_ray_hits = 0;
        for (ptrdiff_t n = 0; n < mesh->n_nodes(); ++n) {
            if (std::fabs(static_cast<double>(pts[2][n] - zmin)) > ztol) {
                continue;
            }
            const double x = pts[0][n];
            const double y = pts[1][n];
            const double r = std::sqrt(x * x + y * y);
            if (r < 1e-12) {
                continue;
            }
            if (std::fabs(std::fabs(x) - std::fabs(y)) > 1e-4 * static_cast<double>(radius)) {
                continue;
            }
            if (r > 0.8 * static_cast<double>(radius)) {
                continue;
            }
            n_corners++;
            const double th = std::atan2(y, x);
            for (ptrdiff_t m = 0; m < mesh->n_nodes(); ++m) {
                if (m == n || std::fabs(static_cast<double>(pts[2][m] - zmin)) > ztol) {
                    continue;
                }
                const double xm = pts[0][m];
                const double ym = pts[1][m];
                const double rm = std::sqrt(xm * xm + ym * ym);
                if (rm <= r + 1e-6) {
                    continue;
                }
                const double thm = std::atan2(ym, xm);
                double dth = thm - th;
                while (dth > M_PI) {
                    dth -= 2 * M_PI;
                }
                while (dth < -M_PI) {
                    dth += 2 * M_PI;
                }
                if (std::fabs(dth) < 1e-3) {
                    n_ray_hits++;
                    break;
                }
            }
        }
        SMESH_TEST_ASSERT(n_corners >= 4);
        SMESH_TEST_EQ(n_ray_hits, n_corners);
    }

    SMESH_TEST_EQ(check_half_face_tables(*mesh, n_wedge > 0 ? 1 : 0), SMESH_TEST_SUCCESS);

    auto skins = skin_sidesets(mesh);
    SMESH_TEST_ASSERT(!skins.empty());

    return SMESH_TEST_SUCCESS;
}

static int test_hex_dominant_cylinder() {
    SMESH_TEST_EQ(test_hex_dominant_cylinder_mesh(1, 2, 2, 16, 3, 0), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(test_hex_dominant_cylinder_mesh(2, 1, 1, 8, 2, -0.5), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(test_hex_dominant_cylinder_mesh(1, 1, 3, 24, 2, 0), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_hex_dominant_cylinder);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}

