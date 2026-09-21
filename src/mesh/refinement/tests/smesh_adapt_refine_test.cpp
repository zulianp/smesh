#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#include "smesh_adapt_refine.hpp"
#include "smesh_curvature.hpp"
#include "smesh_extractions.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_parametrization.hpp"
#include "smesh_smooth.hpp"
#include "smesh_test.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace smesh;

namespace {

constexpr geom_t kTol = sizeof(geom_t) == 8 ? static_cast<geom_t>(2e-1)
                                            : static_cast<geom_t>(3.5e-1);

int on_unit_cube(const geom_t x, const geom_t y, const geom_t z, const geom_t tol) {
    auto near = [&](geom_t v, geom_t t) -> int {
        return v > t - tol && v < t + tol;
    };
    auto in = [&](geom_t v) -> int { return v > -tol && v < static_cast<geom_t>(1) + tol; };
    if (!in(x) || !in(y) || !in(z)) {
        return 0;
    }
    return near(x, 0) || near(x, 1) || near(y, 0) || near(y, 1) || near(z, 0) || near(z, 1);
}

int cube_surface_ok(const std::shared_ptr<Mesh> &m, const geom_t tol) {
    geom_t **p = m->points()->data();
    const int sdim = m->spatial_dimension();
    for (ptrdiff_t i = 0; i < m->n_nodes(); ++i) {
        const geom_t z = sdim >= 3 ? p[2][i] : static_cast<geom_t>(0);
        if (!on_unit_cube(p[0][i], p[1][i], z, tol)) {
            return SMESH_TEST_FAILURE;
        }
    }
    return SMESH_TEST_SUCCESS;
}

int manifold_surface(const std::shared_ptr<Mesh> &m) {
    auto            n2n = m->node_to_node_graph_upper_triangular();
    const count_t  *rp  = n2n->rowptr()->data();
    const idx_t    *ci  = n2n->colidx()->data();
    const ptrdiff_t ne  = m->n_elements(0);
    idx_t **const   el  = m->elements(0)->data();
    const enum ElemType et  = m->element_type(0);
    const int           nxe = (et == QUAD4 || et == QUADSHELL4) ? 4 : 3;
    std::vector<int>    cnt((size_t)n2n->nnz(), 0);
    auto                slot = [&](idx_t a, idx_t b) -> ptrdiff_t {
        if (a > b) {
            const idx_t t = a;
            a             = b;
            b             = t;
        }
        for (count_t k = rp[a]; k < rp[a + 1]; ++k) {
            if (ci[k] == b) {
                return (ptrdiff_t)k;
            }
        }
        return -1;
    };
    for (ptrdiff_t e = 0; e < ne; ++e) {
        for (int s = 0; s < nxe; ++s) {
            const idx_t     a  = el[s][e];
            const idx_t     b  = el[(s + 1) % nxe][e];
            const ptrdiff_t sl = slot(a, b);
            if (sl < 0) {
                return SMESH_TEST_FAILURE;
            }
            cnt[(size_t)sl] += 1;
        }
    }
    for (size_t k = 0; k < cnt.size(); ++k) {
        if (cnt[k] != 0 && cnt[k] != 1 && cnt[k] != 2) {
            return SMESH_TEST_FAILURE;
        }
    }
    return SMESH_TEST_SUCCESS;
}

u64 connectivity_checksum(const std::shared_ptr<Mesh> &m) {
    const ptrdiff_t     ne  = m->n_elements(0);
    const int           nxe = m->n_nodes_per_element(0);
    idx_t **const       el  = m->elements(0)->data();
    u64                 h   = 1469598103934665603ull;
    h ^= (u64)ne;
    h *= 1099511628211ull;
    h ^= (u64)m->n_nodes();
    h *= 1099511628211ull;
    for (ptrdiff_t e = 0; e < ne; ++e) {
        idx_t v[8];
        for (int d = 0; d < nxe; ++d) {
            v[d] = el[d][e];
        }
        for (int i = 1; i < nxe; ++i) {
            const idx_t x = v[i];
            int         j = i;
            while (j > 0 && v[j - 1] > x) {
                v[j] = v[j - 1];
                --j;
            }
            v[j] = x;
        }
        for (int d = 0; d < nxe; ++d) {
            h ^= (u64)(u32)v[d];
            h *= 1099511628211ull;
        }
    }
    return h;
}

std::shared_ptr<Nodeset> nodes_on_sphere(const std::shared_ptr<Mesh> &m,
                                         const geom_t                 radius,
                                         const geom_t                 frac) {
    geom_t **const  p  = m->points()->data();
    const ptrdiff_t nn = m->n_nodes();
    const geom_t    lo = frac * radius;
    ptrdiff_t       n  = 0;
    for (ptrdiff_t i = 0; i < nn; ++i) {
        const geom_t r = std::sqrt(p[0][i] * p[0][i] + p[1][i] * p[1][i] + p[2][i] * p[2][i]);
        n += r > lo;
    }
    auto ids = create_host_buffer<idx_t>((size_t)n);
    ptrdiff_t k = 0;
    for (ptrdiff_t i = 0; i < nn; ++i) {
        const geom_t r = std::sqrt(p[0][i] * p[0][i] + p[1][i] * p[1][i] + p[2][i] * p[2][i]);
        if (r > lo) {
            ids->data()[k++] = (idx_t)i;
        }
    }
    return Nodeset::create(m->comm(), ids);
}

int attach_sphere_param(const std::shared_ptr<Mesh> &m, const geom_t radius) {
    auto ns = nodes_on_sphere(m, radius, static_cast<geom_t>(0.85));
    if (!ns || ns->size() < 8) {
        return SMESH_TEST_FAILURE;
    }
    auto p = SphereParametrization::create(ns, 0, 0, 0, radius);
    m->add_parametrization("sphere", p);
    return p->apply(*m);
}

struct SphereErr {
    geom_t    mx;
    geom_t    rms;
    geom_t    lmax;
    ptrdiff_t n_faces;
    ptrdiff_t n_samp;
};

SphereErr sphere_pl_error(const std::shared_ptr<Mesh> &surf, const geom_t radius) {
    const ptrdiff_t ne  = surf->n_elements(0);
    const int       nxe = surf->n_nodes_per_element(0);
    idx_t **const   el  = surf->elements(0)->data();
    geom_t **const  p   = surf->points()->data();
    const geom_t    lo  = static_cast<geom_t>(0.97) * radius;
    const geom_t    zlo = static_cast<geom_t>(0.08) * radius;
    SphereErr       out{};
    auto rad = [&](geom_t x, geom_t y, geom_t z) -> geom_t {
        return std::sqrt(x * x + y * y + z * z);
    };
    auto acc = [&](geom_t x, geom_t y, geom_t z) {
        const geom_t e = std::abs(rad(x, y, z) - radius);
        if (e > out.mx) {
            out.mx = e;
        }
        out.rms += e * e;
        ++out.n_samp;
    };
    for (ptrdiff_t e = 0; e < ne; ++e) {
        int    on = 1;
        geom_t cx = 0, cy = 0, cz = 0;
        for (int d = 0; d < nxe; ++d) {
            const idx_t i = el[d][e];
            on &= rad(p[0][i], p[1][i], p[2][i]) > lo;
            cx += p[0][i];
            cy += p[1][i];
            cz += p[2][i];
        }
        if (!on) {
            continue;
        }
        const geom_t inv = static_cast<geom_t>(1) / static_cast<geom_t>(nxe);
        if (cz * inv < zlo) {
            continue;
        }
        ++out.n_faces;
        acc(cx * inv, cy * inv, cz * inv);
        for (int s = 0; s < nxe; ++s) {
            const idx_t a = el[s][e];
            const idx_t b = el[(s + 1) % nxe][e];
            const geom_t mx = static_cast<geom_t>(0.5) * (p[0][a] + p[0][b]);
            const geom_t my = static_cast<geom_t>(0.5) * (p[1][a] + p[1][b]);
            const geom_t mz = static_cast<geom_t>(0.5) * (p[2][a] + p[2][b]);
            const geom_t dx = p[0][a] - p[0][b];
            const geom_t dy = p[1][a] - p[1][b];
            const geom_t dz = p[2][a] - p[2][b];
            const geom_t elen = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (elen > out.lmax) {
                out.lmax = elen;
            }
            acc(mx, my, mz);
        }
    }
    if (out.n_samp > 0) {
        out.rms = std::sqrt(out.rms / static_cast<geom_t>(out.n_samp));
    }
    return out;
}

int tets_positive(const std::shared_ptr<Mesh> &m) {
    if (m->element_type(0) != TET4) {
        return SMESH_TEST_SUCCESS;
    }
    idx_t **const  el = m->elements(0)->data();
    geom_t **const p  = m->points()->data();
    for (ptrdiff_t e = 0; e < m->n_elements(0); ++e) {
        const idx_t  a   = el[0][e], b = el[1][e], c = el[2][e], d = el[3][e];
        const geom_t abx = p[0][b] - p[0][a], aby = p[1][b] - p[1][a], abz = p[2][b] - p[2][a];
        const geom_t acx = p[0][c] - p[0][a], acy = p[1][c] - p[1][a], acz = p[2][c] - p[2][a];
        const geom_t adx = p[0][d] - p[0][a], ady = p[1][d] - p[1][a], adz = p[2][d] - p[2][a];
        const geom_t v =
                abx * (acy * adz - acz * ady) + aby * (acz * adx - acx * adz) + abz * (acx * ady - acy * adx);
        const geom_t s = std::fabs(abx) + std::fabs(aby) + std::fabs(abz) + std::fabs(acx) + std::fabs(acy) +
                         std::fabs(acz) + std::fabs(adx) + std::fabs(ady) + std::fabs(adz);
        const geom_t rel = v / (s * s * s + std::numeric_limits<geom_t>::min());
        SMESH_TEST_ASSERT(rel > static_cast<geom_t>(-2e-3));
    }
    return SMESH_TEST_SUCCESS;
}

}  // namespace

static int test_curvature_cube_skin_flat() {
    auto mesh = Mesh::create_cube(Communicator::self(), HEX8, 4, 4, 4, 0, 0, 0, 1, 1, 1);
    auto surf = skin(mesh);
    SMESH_TEST_ASSERT(surf != nullptr);
    geom_t *kappa = (geom_t *)calloc((size_t)surf->n_nodes(), sizeof(geom_t));
    SMESH_TEST_ASSERT(kappa != nullptr);
    const int kerr = mesh_node_curvature<idx_t, geom_t>(surf->element_type(0),
                                                        surf->n_elements(0),
                                                        surf->elements(0)->data(),
                                                        surf->spatial_dimension(),
                                                        surf->n_nodes(),
                                                        surf->points()->data(),
                                                        nullptr,
                                                        kappa);
    SMESH_TEST_EQ(kerr, SMESH_SUCCESS);
    geom_t    mean = 0;
    ptrdiff_t n    = 0;
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        if (kappa[i] < static_cast<geom_t>(0.5)) {
            mean += kappa[i];
            ++n;
        }
    }
    SMESH_TEST_ASSERT(n > 0);
    mean /= (geom_t)n;
    SMESH_TEST_ASSERT(mean < static_cast<geom_t>(0.25));
    free(kappa);
    return SMESH_TEST_SUCCESS;
}

static int test_curvature_sphere_skin() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 6, 6, 3);
    auto surf = skin(mesh);
    SMESH_TEST_ASSERT(surf != nullptr);
    geom_t *kappa = (geom_t *)calloc((size_t)surf->n_nodes(), sizeof(geom_t));
    SMESH_TEST_ASSERT(kappa != nullptr);
    const int kerr = mesh_node_curvature<idx_t, geom_t>(surf->element_type(0),
                                                        surf->n_elements(0),
                                                        surf->elements(0)->data(),
                                                        surf->spatial_dimension(),
                                                        surf->n_nodes(),
                                                        surf->points()->data(),
                                                        nullptr,
                                                        kappa);
    SMESH_TEST_EQ(kerr, SMESH_SUCCESS);
    geom_t **p   = surf->points()->data();
    geom_t   acc = 0;
    ptrdiff_t n  = 0;
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        if (p[2][i] > static_cast<geom_t>(0.2)) {
            acc += kappa[i];
            ++n;
        }
    }
    SMESH_TEST_ASSERT(n > 4);
    const geom_t mean = acc / (geom_t)n;
    SMESH_TEST_ASSERT(mean > static_cast<geom_t>(0.3));
    SMESH_TEST_ASSERT(mean < static_cast<geom_t>(3));
    free(kappa);
    return SMESH_TEST_SUCCESS;
}

static int test_size_field_sphere_uniform() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 5, 5, 2);
    auto surf = skin(mesh);
    SMESH_TEST_ASSERT(surf != nullptr);
    const ptrdiff_t nn    = surf->n_nodes();
    geom_t         *kappa = (geom_t *)calloc((size_t)nn, sizeof(geom_t));
    geom_t         *h     = (geom_t *)calloc((size_t)nn, sizeof(geom_t));
    SMESH_TEST_ASSERT(kappa != nullptr && h != nullptr);
    const int kerr = mesh_node_curvature<idx_t, geom_t>(surf->element_type(0),
                                                        surf->n_elements(0),
                                                        surf->elements(0)->data(),
                                                        surf->spatial_dimension(),
                                                        nn,
                                                        surf->points()->data(),
                                                        nullptr,
                                                        kappa);
    SMESH_TEST_EQ(kerr, SMESH_SUCCESS);
    const int herr = mesh_size_from_curvature(nn,
                                              kappa,
                                              static_cast<geom_t>(8),
                                              static_cast<geom_t>(0),
                                              static_cast<geom_t>(0),
                                              static_cast<geom_t>(2),
                                              h);
    SMESH_TEST_EQ(herr, SMESH_SUCCESS);
    geom_t **p     = surf->points()->data();
    geom_t   hmin  = h[0], hmax = h[0];
    ptrdiff_t ncap = 0;
    for (ptrdiff_t i = 0; i < nn; ++i) {
        if (p[2][i] <= static_cast<geom_t>(0.25)) {
            continue;
        }
        if (h[i] < hmin) {
            hmin = h[i];
        }
        if (h[i] > hmax) {
            hmax = h[i];
        }
        ++ncap;
    }
    SMESH_TEST_ASSERT(ncap > 4);
    SMESH_TEST_ASSERT(hmin > static_cast<geom_t>(0));
    SMESH_TEST_ASSERT(hmax < static_cast<geom_t>(3));
    auto n2n = surf->node_to_node_graph();
    const int gerr = mesh_grade_size_2to1<idx_t, count_t, geom_t>(
            nn, n2n->rowptr()->data(), n2n->colidx()->data(), 4, h);
    SMESH_TEST_EQ(gerr, SMESH_SUCCESS);
    free(kappa);
    free(h);
    return SMESH_TEST_SUCCESS;
}

static int test_size_from_curvature_error() {
    const ptrdiff_t n = 3;
    geom_t          kappa[3] = {0, 1, 4};
    geom_t          h[3]     = {0, 0, 0};
    const geom_t    eps      = static_cast<geom_t>(0.01);
    const int herr = mesh_size_from_curvature_error(n,
                                                    kappa,
                                                    static_cast<geom_t>(8),
                                                    eps,
                                                    static_cast<geom_t>(0),
                                                    static_cast<geom_t>(0),
                                                    static_cast<geom_t>(1),
                                                    h);
    SMESH_TEST_EQ(herr, SMESH_SUCCESS);
    SMESH_TEST_ASSERT(h[1] > static_cast<geom_t>(0));
    SMESH_TEST_ASSERT(h[2] < h[1]);
    SMESH_TEST_ASSERT(h[0] > h[1] * static_cast<geom_t>(10));
    const geom_t eta1 = kappa[1] * h[1] * h[1] / static_cast<geom_t>(8);
    SMESH_TEST_APPROXEQ(eta1, eps, static_cast<geom_t>(1e-5));
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_sphere_curvature_error_only() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 3, 3, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto lmax_skin = [](const std::shared_ptr<Mesh> &s) -> geom_t {
        idx_t **el  = s->elements(0)->data();
        geom_t **p  = s->points()->data();
        const int nxe = s->n_nodes_per_element(0);
        geom_t m = 0;
        for (ptrdiff_t e = 0; e < s->n_elements(0); ++e) {
            geom_t cz = 0;
            for (int d = 0; d < nxe; ++d) {
                cz += p[2][el[d][e]];
            }
            if (cz / (geom_t)nxe < static_cast<geom_t>(0.08)) {
                continue;
            }
            for (int sidx = 0; sidx < nxe; ++sidx) {
                const idx_t a = el[sidx][e], b = el[(sidx + 1) % nxe][e];
                geom_t l2 = 0;
                for (int d = 0; d < 3; ++d) {
                    const geom_t t = p[d][b] - p[d][a];
                    l2 += t * t;
                }
                const geom_t l = std::sqrt(l2);
                if (l > m) {
                    m = l;
                }
            }
        }
        return m;
    };
    const geom_t l0 = lmax_skin(skin(mesh));
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(24);
    opt.geom_error          = static_cast<geom_t>(0.004);
    opt.h_max               = 0;
    opt.max_levels          = 6;
    opt.smooth_iters        = 0;
    opt.use_parametrization = false;
    auto fine = adapt_refine(mesh, opt);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_ASSERT(fine->n_elements() > mesh->n_elements() * 2);
    const geom_t l1 = lmax_skin(skin(fine));
    SMESH_TEST_ASSERT(l1 < l0);
    SMESH_TEST_ASSERT(l1 < static_cast<geom_t>(0.35));
    {
        auto sk0 = skin(mesh);
        auto sk1 = skin(fine);
        SMESH_TEST_ASSERT(sk1->n_elements() > sk0->n_elements());
        SMESH_TEST_ASSERT(fine->n_elements() > sk1->n_elements());
        SMESH_TEST_ASSERT(fine->n_elements() < sk1->n_elements() * 20);
    }
    return SMESH_TEST_SUCCESS;
}

static int test_smooth_cube_corners() {
    auto mesh = Mesh::create_cube(Communicator::self(), HEX8, 1, 1, 1, 0, 0, 0, 1, 1, 1);
    auto surf = skin(mesh);
    AdaptRefineOptions opt;
    opt.smooth_iters        = 20;
    opt.smooth_lambda       = static_cast<geom_t>(0.4);
    opt.use_parametrization = false;
    SMESH_TEST_EQ(smooth_enhance(*surf, opt), SMESH_SUCCESS);
    SMESH_TEST_EQ(cube_surface_ok(surf, static_cast<geom_t>(1e-3)), SMESH_TEST_SUCCESS);
    auto se = extract_sharp_edges(*surf, static_cast<geom_t>(0.15));
    SMESH_TEST_ASSERT(se != nullptr);
    SMESH_TEST_EQ(se->size(), (ptrdiff_t)12);
    auto co = extract_sharp_corners(*surf, se, false);
    SMESH_TEST_ASSERT(co != nullptr);
    SMESH_TEST_EQ(co->size(), (ptrdiff_t)8);
    geom_t **p = surf->points()->data();
    const idx_t *ids = co->nodes()->data();
    for (ptrdiff_t i = 0; i < co->size(); ++i) {
        const idx_t id = ids[i];
        const geom_t x = p[0][id], y = p[1][id], z = p[2][id];
        SMESH_TEST_ASSERT(x < static_cast<geom_t>(0.01) || x > static_cast<geom_t>(0.99));
        SMESH_TEST_ASSERT(y < static_cast<geom_t>(0.01) || y > static_cast<geom_t>(0.99));
        SMESH_TEST_ASSERT(z < static_cast<geom_t>(0.01) || z > static_cast<geom_t>(0.99));
    }
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_hex_rejected() {
    auto mesh = Mesh::create_cube(Communicator::self(), HEX8, 2, 2, 2);
    auto out  = adapt_refine(mesh);
    SMESH_TEST_ASSERT(out == nullptr);
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_cube_skin_flats() {
    auto hex  = Mesh::create_cube(Communicator::self(), HEX8, 3, 3, 3);
    auto surf = skin(hex);
    SMESH_TEST_ASSERT(surf != nullptr);
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(8);
    opt.max_levels          = 3;
    opt.smooth_iters        = 10;
    opt.use_parametrization = false;
    auto fine = adapt_refine(surf, opt);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(manifold_surface(fine), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(cube_surface_ok(fine, static_cast<geom_t>(1e-3)), SMESH_TEST_SUCCESS);
    {
        geom_t **p0 = surf->points()->data();
        geom_t **p1 = fine->points()->data();
        const ptrdiff_t n0 = surf->n_nodes();
        const geom_t    tol = static_cast<geom_t>(1e-5);
        SMESH_TEST_ASSERT(fine->n_nodes() >= n0);
        for (ptrdiff_t i = 0; i < n0; ++i) {
            const geom_t dx = p1[0][i] - p0[0][i];
            const geom_t dy = p1[1][i] - p0[1][i];
            const geom_t dz = p1[2][i] - p0[2][i];
            SMESH_TEST_ASSERT(dx * dx + dy * dy + dz * dz < tol * tol);
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_tri_hump() {
    auto tet = Mesh::create_wall_mounted_hump(Communicator::self(), TET4, 12, 6, 2);
    SMESH_TEST_ASSERT(tet != nullptr);
    auto surf = skin(tet);
    SMESH_TEST_ASSERT(surf != nullptr);
    const ptrdiff_t n0 = surf->n_elements();
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(16);
    opt.max_levels          = 4;
    opt.smooth_iters        = 4;
    opt.use_parametrization = false;
    auto fine = adapt_refine(surf, opt);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_ASSERT(fine->n_elements() >= n0);
    SMESH_TEST_EQ(manifold_surface(fine), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_quad_square() {
    auto hex  = Mesh::create_cube(Communicator::self(), HEX8, 8, 4, 1, 0, 0, 0, 2, 1, 0.2);
    auto surf = skin(hex);
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(12);
    opt.max_levels          = 2;
    opt.smooth_iters        = 2;
    opt.use_parametrization = false;
    auto fine = adapt_refine(surf, opt);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(manifold_surface(fine), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_tet_half_sphere() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 4, 4, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const ptrdiff_t n0 = mesh->n_elements();
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(20);
    opt.max_levels          = 3;
    opt.smooth_iters        = 2;
    opt.use_parametrization = false;
    auto fine = adapt_refine(mesh, opt);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->element_type(0), TET4);
    SMESH_TEST_ASSERT(fine->n_elements() >= n0);
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_sphere_param() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 4, 4, 2);
    auto surf = skin(mesh);
    auto ids  = create_host_buffer<idx_t>((size_t)surf->n_nodes());
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        ids->data()[i] = (idx_t)i;
    }
    auto ns = Nodeset::create(surf->comm(), ids);
    auto p  = SphereParametrization::create(ns, 0, 0, 0, 1);
    surf->add_parametrization("sphere", p);
    p->apply(*surf);
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(16);
    opt.max_levels          = 2;
    opt.smooth_iters        = 0;
    opt.use_parametrization = true;
    auto fine = adapt_refine(surf, opt);
    SMESH_TEST_ASSERT(fine != nullptr);
    geom_t **pts = fine->points()->data();
    ptrdiff_t n_on = 0;
    for (ptrdiff_t i = 0; i < fine->n_nodes(); ++i) {
        if (pts[2][i] <= static_cast<geom_t>(0.05)) {
            continue;
        }
        const geom_t r =
                std::sqrt(pts[0][i] * pts[0][i] + pts[1][i] * pts[1][i] + pts[2][i] * pts[2][i]);
        SMESH_TEST_APPROXEQ(r, 1, kTol);
        ++n_on;
    }
    SMESH_TEST_ASSERT(n_on > 0);
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_tet_sphere_param_error() {
    const geom_t radius = static_cast<geom_t>(1);
    auto         mesh   = Mesh::create_half_sphere(Communicator::self(), TET4, radius, 3, 3, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(attach_sphere_param(mesh, radius), SMESH_SUCCESS);
    SMESH_TEST_EQ(tets_positive(mesh), SMESH_TEST_SUCCESS);

    const SphereErr e0 = sphere_pl_error(skin(mesh), radius);
    SMESH_TEST_ASSERT(e0.n_faces > 4);
    SMESH_TEST_ASSERT(e0.rms > static_cast<geom_t>(1e-4));

    SphereErr e_prev = e0;
    auto      current = mesh;
    for (int step = 0; step < 2; ++step) {
        AdaptRefineOptions opt;
        opt.cells_per_radius    = static_cast<geom_t>(24);
        opt.geom_error          = static_cast<geom_t>(0.004);
        opt.h_max               = 0;
        opt.max_levels          = 6;
        opt.smooth_iters        = 2;
        opt.smooth_lambda       = static_cast<geom_t>(0.5);
        opt.use_parametrization = true;
        const ptrdiff_t n0      = current->n_elements();
        const ptrdiff_t ns0     = skin(current)->n_elements();
        auto            fine    = adapt_refine(current, opt);
        SMESH_TEST_ASSERT(fine != nullptr);
        SMESH_TEST_EQ(fine->element_type(0), TET4);
        SMESH_TEST_ASSERT(fine->n_elements() > n0);
        SMESH_TEST_ASSERT(skin(fine)->n_elements() > ns0);
        SMESH_TEST_EQ(tets_positive(fine), SMESH_TEST_SUCCESS);

        geom_t **const pts = fine->points()->data();
        std::shared_ptr<Nodeset> ns;
        for (const auto &kv : fine->parametrizations()) {
            if (kv.first == "sphere" && kv.second) {
                ns = kv.second->nodeset();
                break;
            }
        }
        SMESH_TEST_ASSERT(ns != nullptr && ns->size() > 0);
        const idx_t *ids = ns->nodes()->data();
        const geom_t vtol =
                sizeof(geom_t) == 8 ? static_cast<geom_t>(1e-8) : static_cast<geom_t>(2e-4);
        ptrdiff_t n_on = 0;
        geom_t    rmin = radius;
        for (ptrdiff_t i = 0; i < ns->size(); ++i) {
            const idx_t  id = ids[i];
            const geom_t r =
                    std::sqrt(pts[0][id] * pts[0][id] + pts[1][id] * pts[1][id] + pts[2][id] * pts[2][id]);
            if (r < rmin) {
                rmin = r;
            }
            n_on += std::fabs(r - radius) <= vtol;
        }
        SMESH_TEST_ASSERT(n_on > (ns->size() * 3) / 4);
        SMESH_TEST_ASSERT(rmin > static_cast<geom_t>(0.90) * radius);

        const SphereErr e = sphere_pl_error(skin(fine), radius);
        if (e.n_faces > e_prev.n_faces) {
            SMESH_TEST_ASSERT(e.rms < e_prev.rms);
            SMESH_TEST_ASSERT(e.mx < e_prev.mx);
            SMESH_TEST_ASSERT(e.lmax < e_prev.lmax);
        } else {
            SMESH_TEST_ASSERT(e.rms <= e_prev.rms);
        }
        e_prev  = e;
        current = fine;
    }
    SMESH_TEST_ASSERT(e_prev.n_faces > e0.n_faces);
    SMESH_TEST_ASSERT(e_prev.rms < static_cast<geom_t>(0.35) * e0.rms);
    SMESH_TEST_ASSERT(e_prev.lmax < static_cast<geom_t>(0.20));
    SMESH_TEST_ASSERT(skin(current)->n_elements() > skin(mesh)->n_elements() * 2);
    SMESH_TEST_ASSERT(current->n_elements() < skin(current)->n_elements() * 20);
    return SMESH_TEST_SUCCESS;
}

static int test_adapt_omp_determinism() {
    auto make = []() {
        return Mesh::create_wall_mounted_hump(Communicator::self(), TET4, 8, 4, 2);
    };
    AdaptRefineOptions opt;
    opt.cells_per_radius    = static_cast<geom_t>(12);
    opt.max_levels          = 2;
    opt.smooth_iters        = 0;
    opt.use_parametrization = false;
#ifdef _OPENMP
    const int saved = omp_get_max_threads();
    omp_set_num_threads(1);
    auto a = adapt_refine(skin(make()), opt);
    omp_set_num_threads(saved > 1 ? saved : 4);
    auto b = adapt_refine(skin(make()), opt);
    omp_set_num_threads(saved);
#else
    auto a = adapt_refine(skin(make()), opt);
    auto b = adapt_refine(skin(make()), opt);
#endif
    SMESH_TEST_ASSERT(a != nullptr && b != nullptr);
    SMESH_TEST_EQ(a->n_nodes(), b->n_nodes());
    SMESH_TEST_EQ(a->n_elements(), b->n_elements());
    SMESH_TEST_EQ(connectivity_checksum(a), connectivity_checksum(b));
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_curvature_cube_skin_flat);
    SMESH_RUN_TEST(test_curvature_sphere_skin);
    SMESH_RUN_TEST(test_size_field_sphere_uniform);
    SMESH_RUN_TEST(test_size_from_curvature_error);
    SMESH_RUN_TEST(test_smooth_cube_corners);
    SMESH_RUN_TEST(test_adapt_hex_rejected);
    SMESH_RUN_TEST(test_adapt_cube_skin_flats);
    SMESH_RUN_TEST(test_adapt_tri_hump);
    SMESH_RUN_TEST(test_adapt_quad_square);
    SMESH_RUN_TEST(test_adapt_tet_half_sphere);
    SMESH_RUN_TEST(test_adapt_sphere_param);
    SMESH_RUN_TEST(test_adapt_sphere_curvature_error_only);
    SMESH_RUN_TEST(test_adapt_tet_sphere_param_error);
    SMESH_RUN_TEST(test_adapt_omp_determinism);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
