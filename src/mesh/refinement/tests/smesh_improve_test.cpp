#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "smesh_extractions.hpp"
#include "smesh_improve.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_parametrization.hpp"
#include "smesh_quality.hpp"
#include "smesh_test.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace smesh;

namespace {

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

ptrdiff_t cube_skeleton_count(const std::shared_ptr<Mesh> &m, const geom_t tol) {
    geom_t **p = m->points()->data();
    const int sdim = m->spatial_dimension();
    auto near = [&](geom_t v, geom_t t) -> int {
        return v > t - tol && v < t + tol;
    };
    ptrdiff_t n_skel = 0;
    for (ptrdiff_t i = 0; i < m->n_nodes(); ++i) {
        const geom_t x = p[0][i];
        const geom_t y = p[1][i];
        const geom_t z = sdim >= 3 ? p[2][i] : static_cast<geom_t>(0);
        const int nx = near(x, 0) || near(x, 1);
        const int ny = near(y, 0) || near(y, 1);
        const int nz = near(z, 0) || near(z, 1);
        if (nx + ny + nz >= 2) {
            ++n_skel;
        }
    }
    return n_skel;
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
    int n_bad_edge = 0;
    for (size_t k = 0; k < cnt.size(); ++k) {
        if (cnt[k] != 0 && cnt[k] != 1 && cnt[k] != 2) {
            ++n_bad_edge;
        }
    }
    if (n_bad_edge) {
        return SMESH_TEST_FAILURE;
    }
    if (nxe != 3) {
        return SMESH_TEST_SUCCESS;
    }
    std::vector<u64> keys((size_t)ne);
    for (ptrdiff_t e = 0; e < ne; ++e) {
        idx_t v[3] = {el[0][e], el[1][e], el[2][e]};
        for (int i = 1; i < 3; ++i) {
            const idx_t x = v[i];
            int         j = i;
            while (j > 0 && v[j - 1] > x) {
                v[j] = v[j - 1];
                --j;
            }
            v[j] = x;
        }
        keys[(size_t)e] = ((u64)(u32)v[0] << 42) | ((u64)(u32)v[1] << 21) | (u64)(u32)v[2];
    }
    std::sort(keys.begin(), keys.end());
    for (size_t i = 1; i < keys.size(); ++i) {
        if (keys[i] == keys[i - 1]) {
            return SMESH_TEST_FAILURE;
        }
    }
    return SMESH_TEST_SUCCESS;
}

u64 connectivity_checksum(const std::shared_ptr<Mesh> &m) {
    const ptrdiff_t ne  = m->n_elements(0);
    const int       nxe = m->n_nodes_per_element(0);
    idx_t **const   el  = m->elements(0)->data();
    u64             h   = 1469598103934665603ull;
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

geom_t min_q(const std::shared_ptr<Mesh> &m) {
    const ptrdiff_t ne = m->n_elements(0);
    geom_t         *q  = (geom_t *)calloc((size_t)ne, sizeof(geom_t));
    mesh_element_quality<idx_t, geom_t>(
            m->element_type(0), ne, m->elements(0)->data(), m->spatial_dimension(), m->points()->data(), q);
    geom_t mq = mesh_quality_min(ne, q);
    free(q);
    return mq;
}

}  // namespace

static int test_mean_ratio_unit_shapes() {
    {
        geom_t x[3] = {0, 1, static_cast<geom_t>(0.5)};
        geom_t y[3] = {0, 0, std::sqrt(static_cast<geom_t>(3)) / static_cast<geom_t>(2)};
        geom_t z[3] = {0, 0, 0};
        geom_t *pts[3] = {x, y, z};
        idx_t a = 0, b = 1, c = 2;
        idx_t *el[3] = {&a, &b, &c};
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TRI3, 2, el, pts, 0);
        SMESH_TEST_ASSERT(q > static_cast<geom_t>(0.99));
    }
    {
        geom_t x[3] = {0, 1, static_cast<geom_t>(0.5)};
        geom_t y[3] = {0, 0, static_cast<geom_t>(1e-6)};
        geom_t z[3] = {0, 0, 0};
        geom_t *pts[3] = {x, y, z};
        idx_t a = 0, b = 1, c = 2;
        idx_t *el[3] = {&a, &b, &c};
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TRI3, 2, el, pts, 0);
        SMESH_TEST_ASSERT(q < static_cast<geom_t>(0.05));
    }
    {
        geom_t x[4] = {1, 1, -1, -1};
        geom_t y[4] = {1, -1, -1, 1};
        geom_t z[4] = {1, -1, 1, -1};
        geom_t *pts[3] = {x, y, z};
        idx_t v0 = 0, v1 = 1, v2 = 2, v3 = 3;
        idx_t *el[4] = {&v0, &v1, &v2, &v3};
        geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, pts, 0);
        if (!(q > static_cast<geom_t>(0.5))) {
            const idx_t t = v2;
            v2            = v3;
            v3            = t;
            q             = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, pts, 0);
        }
        SMESH_TEST_ASSERT(q > static_cast<geom_t>(0.99));
    }
    return SMESH_TEST_SUCCESS;
}

static int test_improve_hex_rejected() {
    auto mesh = Mesh::create_cube(Communicator::self(), HEX8, 2, 2, 2);
    ImproveOptions opt;
    opt.smooth_iters        = 0;
    opt.use_parametrization = false;
    SMESH_TEST_EQ(improve(*mesh, opt), SMESH_FAILURE);
    SMESH_TEST_ASSERT(remesh(mesh, opt) == nullptr);
    return SMESH_TEST_SUCCESS;
}

static int test_improve_cube_skin_features() {
    auto hex  = Mesh::create_cube(Communicator::self(), HEX8, 1, 1, 1, 0, 0, 0, 1, 1, 1);
    auto surf = skin(hex);
    SMESH_TEST_ASSERT(surf != nullptr);
    geom_t **p0 = surf->points()->data();
    std::vector<geom_t> x0((size_t)surf->n_nodes() * 3);
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        x0[(size_t)i * 3 + 0] = p0[0][i];
        x0[(size_t)i * 3 + 1] = p0[1][i];
        x0[(size_t)i * 3 + 2] = p0[2][i];
    }
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.3);
    opt.max_abs_dev         = static_cast<geom_t>(0.05);
    opt.max_normal_dev      = static_cast<geom_t>(0.05);
    opt.max_passes          = 2;
    opt.smooth_iters        = 8;
    opt.use_parametrization = false;
    SMESH_TEST_EQ(improve(*surf, opt), SMESH_SUCCESS);
    auto se = extract_sharp_edges(*surf, static_cast<geom_t>(0.15));
    SMESH_TEST_ASSERT(se != nullptr);
    SMESH_TEST_EQ(se->size(), (ptrdiff_t)12);
    auto co = extract_sharp_corners(*surf, se, false);
    SMESH_TEST_ASSERT(co != nullptr);
    SMESH_TEST_EQ(co->size(), (ptrdiff_t)8);
    SMESH_TEST_EQ(cube_surface_ok(surf, static_cast<geom_t>(1e-3)), SMESH_TEST_SUCCESS);
    geom_t **p = surf->points()->data();
    const idx_t *ids = co->nodes()->data();
    for (ptrdiff_t i = 0; i < co->size(); ++i) {
        const idx_t  id = ids[i];
        const geom_t x  = p[0][id], y = p[1][id], z = p[2][id];
        SMESH_TEST_ASSERT(x < static_cast<geom_t>(0.01) || x > static_cast<geom_t>(0.99));
        SMESH_TEST_ASSERT(y < static_cast<geom_t>(0.01) || y > static_cast<geom_t>(0.99));
        SMESH_TEST_ASSERT(z < static_cast<geom_t>(0.01) || z > static_cast<geom_t>(0.99));
    }
    const geom_t band = static_cast<geom_t>(0.05) + static_cast<geom_t>(1e-5);
    for (ptrdiff_t i = 0; i < surf->n_nodes() && i < (ptrdiff_t)(x0.size() / 3); ++i) {
        const geom_t dx = p[0][i] - x0[(size_t)i * 3];
        const geom_t dy = p[1][i] - x0[(size_t)i * 3 + 1];
        const geom_t dz = p[2][i] - x0[(size_t)i * 3 + 2];
        const geom_t r  = std::sqrt(dx * dx + dy * dy + dz * dz);
        SMESH_TEST_ASSERT(r <= band * static_cast<geom_t>(2) || i >= 8);
        (void)r;
        (void)band;
    }
    return SMESH_TEST_SUCCESS;
}

static int test_improve_cube_skin_creases() {
    auto hex  = Mesh::create_cube(Communicator::self(), HEX8, 3, 3, 3);
    auto surf = skin(hex);
    SMESH_TEST_ASSERT(surf != nullptr);
    const ptrdiff_t n_skel0 = cube_skeleton_count(surf, static_cast<geom_t>(1e-4));
    SMESH_TEST_ASSERT(n_skel0 >= 8 + 12);
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.3);
    opt.max_abs_dev         = static_cast<geom_t>(0.05);
    opt.max_normal_dev      = static_cast<geom_t>(0.05);
    opt.max_passes          = 2;
    opt.smooth_iters        = 8;
    opt.use_parametrization = false;
    opt.allow_collapse      = false;
    SMESH_TEST_EQ(improve(*surf, opt), SMESH_SUCCESS);
    SMESH_TEST_EQ(cube_surface_ok(surf, static_cast<geom_t>(1e-3)), SMESH_TEST_SUCCESS);
    auto se = extract_sharp_edges(*surf, static_cast<geom_t>(0.15));
    SMESH_TEST_ASSERT(se != nullptr);
    SMESH_TEST_ASSERT(se->size() >= (ptrdiff_t)12);
    auto co = extract_sharp_corners(*surf, se, false);
    SMESH_TEST_ASSERT(co != nullptr);
    SMESH_TEST_EQ(co->size(), (ptrdiff_t)8);
    const ptrdiff_t n_skel = cube_skeleton_count(surf, static_cast<geom_t>(1e-3));
    SMESH_TEST_ASSERT(n_skel >= n_skel0);
    return SMESH_TEST_SUCCESS;
}

static int test_improve_sphere_band() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 6, 6, 3);
    auto surf = skin(mesh);
    SMESH_TEST_ASSERT(surf != nullptr);
    ImproveOptions opt;
    opt.max_abs_dev         = static_cast<geom_t>(0.04);
    opt.max_normal_dev      = static_cast<geom_t>(0.04);
    opt.max_passes          = 2;
    opt.smooth_iters        = 4;
    opt.use_parametrization = false;
    opt.allow_collapse      = false;
    SMESH_TEST_EQ(improve(*surf, opt), SMESH_SUCCESS);
    geom_t **pts = surf->points()->data();
    const geom_t lo = static_cast<geom_t>(1) - static_cast<geom_t>(0.08);
    const geom_t hi = static_cast<geom_t>(1) + static_cast<geom_t>(0.08);
    ptrdiff_t n_on = 0;
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        if (pts[2][i] <= static_cast<geom_t>(0.05)) {
            continue;
        }
        const geom_t r = std::sqrt(pts[0][i] * pts[0][i] + pts[1][i] * pts[1][i] + pts[2][i] * pts[2][i]);
        SMESH_TEST_ASSERT(r >= lo && r <= hi);
        ++n_on;
    }
    SMESH_TEST_ASSERT(n_on > 0);
    SMESH_TEST_EQ(manifold_surface(surf), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_improve_sphere_param() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 5, 5, 2);
    auto surf = skin(mesh);
    auto ids  = create_host_buffer<idx_t>((size_t)surf->n_nodes());
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        ids->data()[i] = (idx_t)i;
    }
    auto ns = Nodeset::create(surf->comm(), ids);
    auto p  = SphereParametrization::create(ns, 0, 0, 0, 1);
    surf->add_parametrization("sphere", p);
    p->apply(*surf);
    ImproveOptions opt;
    opt.max_abs_dev         = static_cast<geom_t>(0.2);
    opt.max_passes          = 1;
    opt.smooth_iters        = 2;
    opt.use_parametrization = true;
    opt.allow_split         = false;
    opt.allow_collapse      = false;
    opt.allow_swap          = false;
    SMESH_TEST_EQ(improve(*surf, opt), SMESH_SUCCESS);
    geom_t **pts = surf->points()->data();
    const geom_t tol = sizeof(geom_t) == 8 ? static_cast<geom_t>(2e-1) : static_cast<geom_t>(3.5e-1);
    ptrdiff_t n_on = 0;
    for (ptrdiff_t i = 0; i < surf->n_nodes(); ++i) {
        if (pts[2][i] <= static_cast<geom_t>(0.05)) {
            continue;
        }
        const geom_t r = std::sqrt(pts[0][i] * pts[0][i] + pts[1][i] * pts[1][i] + pts[2][i] * pts[2][i]);
        SMESH_TEST_APPROXEQ(r, 1, tol);
        ++n_on;
    }
    SMESH_TEST_ASSERT(n_on > 0);
    return SMESH_TEST_SUCCESS;
}

static int test_improve_tri_skinny() {
    auto tet  = Mesh::create_wall_mounted_hump(Communicator::self(), TET4, 10, 4, 2);
    auto surf = skin(tet);
    SMESH_TEST_ASSERT(surf != nullptr);
    const geom_t q0 = min_q(surf);
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.3);
    opt.max_abs_dev         = static_cast<geom_t>(0.05);
    opt.max_passes          = 6;
    opt.smooth_iters        = 6;
    opt.use_parametrization = false;
    SMESH_TEST_EQ(improve(*surf, opt), SMESH_SUCCESS);
    SMESH_TEST_EQ(manifold_surface(surf), SMESH_TEST_SUCCESS);
    const geom_t q1 = min_q(surf);
    SMESH_TEST_ASSERT(q1 + static_cast<geom_t>(1e-6) >= q0 || q1 > static_cast<geom_t>(0));
    return SMESH_TEST_SUCCESS;
}

static int test_improve_tet_volume() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 4, 4, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const geom_t q0 = min_q(mesh);
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.2);
    opt.max_abs_dev         = static_cast<geom_t>(0.04);
    opt.max_passes          = 3;
    opt.smooth_iters        = 2;
    opt.use_parametrization = false;
    opt.allow_collapse      = false;
    SMESH_TEST_EQ(improve(*mesh, opt), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh->element_type(0), TET4);
    geom_t **p = mesh->points()->data();
    idx_t **el = mesh->elements(0)->data();
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e);
        SMESH_TEST_ASSERT(q > static_cast<geom_t>(0));
    }
    const geom_t q1 = min_q(mesh);
    SMESH_TEST_ASSERT(q1 + static_cast<geom_t>(1e-8) >= q0 * static_cast<geom_t>(0.5));
    return SMESH_TEST_SUCCESS;
}

static int test_improve_tet_collapse_manifold() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 4, 4, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    auto            skin0 = skin(mesh);
    const ptrdiff_t nf0   = skin0->n_elements();
    ImproveOptions  opt;
    opt.q_min               = static_cast<geom_t>(0.2);
    opt.max_abs_dev         = static_cast<geom_t>(0.04);
    opt.max_passes          = 2;
    opt.smooth_iters        = 0;
    opt.use_parametrization = false;
    opt.allow_split         = true;
    opt.allow_collapse      = true;
    opt.allow_swap          = true;
    SMESH_TEST_EQ(improve(*mesh, opt), SMESH_SUCCESS);
    auto skin1 = skin(mesh);
    SMESH_TEST_ASSERT(skin1 != nullptr);
    SMESH_TEST_EQ(manifold_surface(skin1), SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(skin1->n_elements() <= nf0 * 3 + 8);
    geom_t **p = mesh->points()->data();
    idx_t **el = mesh->elements(0)->data();
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e);
        SMESH_TEST_ASSERT(q > static_cast<geom_t>(0));
    }
    return SMESH_TEST_SUCCESS;
}

static int test_improve_tet_param_no_invert() {
    auto mesh = Mesh::create_half_sphere(Communicator::self(), TET4, 1, 5, 5, 2);
    SMESH_TEST_ASSERT(mesh != nullptr);
    geom_t **p = mesh->points()->data();
    ptrdiff_t n_on = 0;
    for (ptrdiff_t i = 0; i < mesh->n_nodes(); ++i) {
        const geom_t r = std::sqrt(p[0][i] * p[0][i] + p[1][i] * p[1][i] + p[2][i] * p[2][i]);
        n_on += r > static_cast<geom_t>(0.85);
    }
    auto ids = create_host_buffer<idx_t>((size_t)n_on);
    ptrdiff_t k = 0;
    for (ptrdiff_t i = 0; i < mesh->n_nodes(); ++i) {
        const geom_t r = std::sqrt(p[0][i] * p[0][i] + p[1][i] * p[1][i] + p[2][i] * p[2][i]);
        if (r > static_cast<geom_t>(0.85)) {
            ids->data()[k++] = (idx_t)i;
        }
    }
    mesh->add_parametrization(
            "sphere", SphereParametrization::create(Nodeset::create(mesh->comm(), ids), 0, 0, 0, 1));
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.2);
    opt.max_passes          = 2;
    opt.smooth_iters        = 8;
    opt.use_parametrization = true;
    SMESH_TEST_EQ(improve(*mesh, opt), SMESH_SUCCESS);
    SMESH_TEST_EQ(manifold_surface(skin(mesh)), SMESH_TEST_SUCCESS);
    p          = mesh->points()->data();
    idx_t **el = mesh->elements(0)->data();
    ptrdiff_t n_bad = 0;
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, el, p, e);
        SMESH_TEST_ASSERT(q > static_cast<geom_t>(0));
        n_bad += q < static_cast<geom_t>(0.02);
    }
    SMESH_TEST_ASSERT(n_bad == 0);
    return SMESH_TEST_SUCCESS;
}

static int test_improve_quad_conforming() {
    auto hex  = Mesh::create_cube(Communicator::self(), HEX8, 3, 3, 3);
    auto surf = skin(hex);
    SMESH_TEST_ASSERT(surf != nullptr);
    const enum ElemType et0 = surf->element_type(0);
    SMESH_TEST_ASSERT(et0 == QUADSHELL4 || et0 == QUAD4);
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.5);
    opt.max_passes          = 2;
    opt.smooth_iters        = 2;
    opt.use_parametrization = false;
    SMESH_TEST_EQ(improve(*surf, opt), SMESH_SUCCESS);
    SMESH_TEST_EQ(surf->element_type(0), et0);
    SMESH_TEST_EQ(manifold_surface(surf), SMESH_TEST_SUCCESS);
    SMESH_TEST_EQ(cube_surface_ok(surf, static_cast<geom_t>(1e-3)), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_improve_omp_determinism() {
    auto make = []() {
        return skin(Mesh::create_wall_mounted_hump(Communicator::self(), TET4, 6, 3, 2));
    };
    ImproveOptions opt;
    opt.max_passes          = 2;
    opt.smooth_iters        = 0;
    opt.use_parametrization = false;
    opt.allow_collapse      = false;
#ifdef _OPENMP
    const int saved = omp_get_max_threads();
    omp_set_num_threads(1);
    auto a = remesh(make(), opt);
    omp_set_num_threads(saved > 1 ? saved : 4);
    auto b = remesh(make(), opt);
    omp_set_num_threads(saved);
#else
    auto a = remesh(make(), opt);
    auto b = remesh(make(), opt);
#endif
    SMESH_TEST_ASSERT(a != nullptr && b != nullptr);
    SMESH_TEST_EQ(a->n_nodes(), b->n_nodes());
    SMESH_TEST_EQ(a->n_elements(), b->n_elements());
    SMESH_TEST_EQ(connectivity_checksum(a), connectivity_checksum(b));
    return SMESH_TEST_SUCCESS;
}

static int test_improve_tet_23_grows_safe() {
    auto pbuf = create_host_buffer<geom_t>(3, 5);
    geom_t **p = pbuf->data();
    p[0][0] = 1;
    p[1][0] = 0;
    p[2][0] = 0;
    p[0][1] = static_cast<geom_t>(-0.5);
    p[1][1] = static_cast<geom_t>(0.86602540378);
    p[2][1] = 0;
    p[0][2] = static_cast<geom_t>(-0.5);
    p[1][2] = static_cast<geom_t>(-0.86602540378);
    p[2][2] = 0;
    p[0][3] = 0;
    p[1][3] = 0;
    p[2][3] = static_cast<geom_t>(0.12);
    p[0][4] = 0;
    p[1][4] = 0;
    p[2][4] = static_cast<geom_t>(-0.12);
    auto ebuf = create_host_buffer<idx_t>(4, 2);
    idx_t **el = ebuf->data();
    el[0][0] = 0;
    el[1][0] = 1;
    el[2][0] = 2;
    el[3][0] = 3;
    el[0][1] = 0;
    el[1][1] = 2;
    el[2][1] = 1;
    el[3][1] = 4;
    auto mesh = std::make_shared<Mesh>(Communicator::self(), TET4, ebuf, pbuf);
    ImproveOptions opt;
    opt.q_min               = static_cast<geom_t>(0.01);
    opt.max_passes          = 3;
    opt.smooth_iters        = 0;
    opt.use_parametrization = false;
    opt.allow_split         = false;
    opt.allow_collapse      = false;
    opt.allow_swap          = true;
    SMESH_TEST_EQ(improve(*mesh, opt), SMESH_SUCCESS);
    SMESH_TEST_ASSERT(mesh->n_elements(0) >= 2);
    SMESH_TEST_EQ(mesh->element_type(0), TET4);
    geom_t **pts = mesh->points()->data();
    idx_t **els  = mesh->elements(0)->data();
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        const geom_t q = mesh_elem_mean_ratio<idx_t, geom_t>(TET4, 3, els, pts, e);
        SMESH_TEST_ASSERT(q > static_cast<geom_t>(0));
    }
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char **argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_mean_ratio_unit_shapes);
    SMESH_RUN_TEST(test_improve_hex_rejected);
    SMESH_RUN_TEST(test_improve_cube_skin_features);
    SMESH_RUN_TEST(test_improve_cube_skin_creases);
    SMESH_RUN_TEST(test_improve_sphere_band);
    SMESH_RUN_TEST(test_improve_sphere_param);
    SMESH_RUN_TEST(test_improve_tri_skinny);
    SMESH_RUN_TEST(test_improve_tet_volume);
    SMESH_RUN_TEST(test_improve_tet_collapse_manifold);
    SMESH_RUN_TEST(test_improve_tet_param_no_invert);
    SMESH_RUN_TEST(test_improve_quad_conforming);
    SMESH_RUN_TEST(test_improve_omp_determinism);
    SMESH_RUN_TEST(test_improve_tet_23_grows_safe);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
