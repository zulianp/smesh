#include "smesh_edgeset.hpp"
#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_sideset.hpp"
#include "smesh_test.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace smesh;

static const int k_tet_face[4][3] = {{0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1}};
static const int k_tet_edge[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};

static std::shared_ptr<Mesh> make_tets(const idx_t *conn, ptrdiff_t n_elem, const geom_t *xyz, ptrdiff_t n_nodes) {
    auto elems = create_host_buffer<idx_t>(4, (size_t)n_elem);
    auto pts   = create_host_buffer<geom_t>(3, (size_t)n_nodes);
    for (int k = 0; k < 4; ++k) {
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
            elems->data()[k][e] = conn[(ptrdiff_t)k * n_elem + e];
        }
    }
    for (int d = 0; d < 3; ++d) {
        for (ptrdiff_t i = 0; i < n_nodes; ++i) {
            pts->data()[d][i] = xyz[(ptrdiff_t)d * n_nodes + i];
        }
    }
    auto mesh = std::make_shared<Mesh>(Communicator::self(), TET4, elems, pts);
    mesh->block(0)->set_name("tet");
    return mesh;
}

static std::shared_ptr<Mesh> unit_tet() {
    const idx_t  conn[4]  = {0, 1, 2, 3};
    const geom_t xyz[12]  = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    return make_tets(conn, 1, xyz, 4);
}

static std::shared_ptr<Mesh> face_pair() {
    const idx_t  conn[8] = {0, 0, 1, 2, 2, 1, 3, 4};
    const geom_t xyz[15] = {0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, -1};
    return make_tets(conn, 2, xyz, 5);
}

static std::shared_ptr<Mesh> edge_pair() {
    const idx_t  conn[8] = {0, 0, 1, 1, 2, 4, 3, 5};
    const geom_t xyz[18] = {0, 1, 0, 0, 2, 3, 0, 0, 1, 0, 3, 1, 0, 0, 0, 1, 1, 2};
    return make_tets(conn, 2, xyz, 6);
}

static double tet_volume(const geom_t *x, const geom_t *y, const geom_t *z, idx_t a, idx_t b, idx_t c, idx_t d) {
    const double ax = (double)x[a];
    const double ay = (double)y[a];
    const double az = (double)z[a];
    const double bx = (double)x[b] - ax;
    const double by = (double)y[b] - ay;
    const double bz = (double)z[b] - az;
    const double cx = (double)x[c] - ax;
    const double cy = (double)y[c] - ay;
    const double cz = (double)z[c] - az;
    const double dx = (double)x[d] - ax;
    const double dy = (double)y[d] - ay;
    const double dz = (double)z[d] - az;
    return bx * (cy * dz - cz * dy) - by * (cx * dz - cz * dx) + bz * (cx * dy - cy * dx);
}

static int perm_sign3(idx_t a, idx_t b, idx_t c) {
    idx_t v[3] = {a, b, c};
    int   sign = 1;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (v[i] > v[j]) {
                const idx_t t = v[i];
                v[i]          = v[j];
                v[j]          = t;
                sign          = -sign;
            }
        }
    }
    return sign;
}

static int check_conforming(const std::shared_ptr<Mesh> &mesh, const double volume, const double tol = 1e-8) {
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_ASSERT(mesh->n_blocks() >= 1);
    const ptrdiff_t n_nodes = mesh->n_nodes();
    const geom_t   *x       = mesh->points()->data()[0];
    const geom_t   *y       = mesh->points()->data()[1];
    const geom_t   *z       = mesh->points()->data()[2];

    double sum = 0;
    struct FaceRec {
        idx_t n[3];
        int   sign;
    };
    std::vector<FaceRec> faces;
    for (size_t bi = 0; bi < mesh->n_blocks(); ++bi) {
        SMESH_TEST_EQ(mesh->element_type(static_cast<block_idx_t>(bi)), TET4);
        const ptrdiff_t n_elem = mesh->n_elements(static_cast<block_idx_t>(bi));
        const idx_t    *e0     = mesh->elements(static_cast<block_idx_t>(bi))->data()[0];
        const idx_t    *e1     = mesh->elements(static_cast<block_idx_t>(bi))->data()[1];
        const idx_t    *e2     = mesh->elements(static_cast<block_idx_t>(bi))->data()[2];
        const idx_t    *e3     = mesh->elements(static_cast<block_idx_t>(bi))->data()[3];
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
        const idx_t v[4] = {e0[e], e1[e], e2[e], e3[e]};
        for (int k = 0; k < 4; ++k) {
            SMESH_TEST_ASSERT(v[k] >= 0 && (ptrdiff_t)v[k] < n_nodes);
        }
        const double vol = tet_volume(x, y, z, v[0], v[1], v[2], v[3]);
        SMESH_TEST_ASSERT(vol > 1e-12);
        sum += vol;
        for (int f = 0; f < 4; ++f) {
            idx_t tri[3] = {v[k_tet_face[f][0]], v[k_tet_face[f][1]], v[k_tet_face[f][2]]};
            const int sign = perm_sign3(tri[0], tri[1], tri[2]);
            std::sort(tri, tri + 3);
            faces.push_back(FaceRec{{tri[0], tri[1], tri[2]}, sign});
        }
        for (int ed = 0; ed < 6; ++ed) {
            const idx_t a = v[k_tet_edge[ed][0]];
            const idx_t b = v[k_tet_edge[ed][1]];
            const double mx = 0.5 * ((double)x[a] + (double)x[b]);
            const double my = 0.5 * ((double)y[a] + (double)y[b]);
            const double mz = 0.5 * ((double)z[a] + (double)z[b]);
            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                if (i == (ptrdiff_t)a || i == (ptrdiff_t)b) {
                    continue;
                }
                const double dx = (double)x[i] - mx;
                const double dy = (double)y[i] - my;
                const double dz = (double)z[i] - mz;
                SMESH_TEST_ASSERT(dx * dx + dy * dy + dz * dz > 1e-16);
            }
        }
    }
    }
    SMESH_TEST_APPROXEQ(sum, volume, tol * (std::fabs(volume) + 1.0));

    std::sort(faces.begin(), faces.end(), [](const FaceRec &a, const FaceRec &b) {
        if (a.n[0] != b.n[0]) {
            return a.n[0] < b.n[0];
        }
        if (a.n[1] != b.n[1]) {
            return a.n[1] < b.n[1];
        }
        return a.n[2] < b.n[2];
    });
    for (size_t i = 0; i < faces.size();) {
        size_t j = i + 1;
        while (j < faces.size() && faces[j].n[0] == faces[i].n[0] && faces[j].n[1] == faces[i].n[1] &&
               faces[j].n[2] == faces[i].n[2]) {
            ++j;
        }
        const size_t count = j - i;
        SMESH_TEST_ASSERT(count == 1 || count == 2);
        if (count == 2) {
            SMESH_TEST_ASSERT(faces[i].sign == -faces[i + 1].sign);
        }
        i = j;
    }
    return SMESH_TEST_SUCCESS;
}

static int count_shared_faces(const std::shared_ptr<Mesh> &mesh) {
    const ptrdiff_t n_elem = mesh->n_elements(0);
    const idx_t    *e0     = mesh->elements(0)->data()[0];
    const idx_t    *e1     = mesh->elements(0)->data()[1];
    const idx_t    *e2     = mesh->elements(0)->data()[2];
    const idx_t    *e3     = mesh->elements(0)->data()[3];
    struct Key {
        idx_t n[3];
    };
    std::vector<Key> faces;
    faces.reserve((size_t)n_elem * 4);
    for (ptrdiff_t e = 0; e < n_elem; ++e) {
        const idx_t v[4] = {e0[e], e1[e], e2[e], e3[e]};
        for (int f = 0; f < 4; ++f) {
            idx_t tri[3] = {v[k_tet_face[f][0]], v[k_tet_face[f][1]], v[k_tet_face[f][2]]};
            std::sort(tri, tri + 3);
            faces.push_back(Key{{tri[0], tri[1], tri[2]}});
        }
    }
    std::sort(faces.begin(), faces.end(), [](const Key &a, const Key &b) {
        if (a.n[0] != b.n[0]) {
            return a.n[0] < b.n[0];
        }
        if (a.n[1] != b.n[1]) {
            return a.n[1] < b.n[1];
        }
        return a.n[2] < b.n[2];
    });
    int shared = 0;
    for (size_t i = 0; i + 1 < faces.size(); ++i) {
        if (faces[i].n[0] == faces[i + 1].n[0] && faces[i].n[1] == faces[i + 1].n[1] &&
            faces[i].n[2] == faces[i + 1].n[2]) {
            ++shared;
            ++i;
        }
    }
    return shared;
}

static int same_mesh(const std::shared_ptr<Mesh> &a, const std::shared_ptr<Mesh> &b) {
    SMESH_TEST_ASSERT(a != nullptr && b != nullptr);
    SMESH_TEST_EQ(a->n_elements(0), b->n_elements(0));
    SMESH_TEST_EQ(a->n_nodes(), b->n_nodes());
    for (int k = 0; k < 4; ++k) {
        SMESH_ASSERT_ARRAY_EQ(a->n_elements(0), a->elements(0)->data()[k], b->elements(0)->data()[k]);
    }
    for (int d = 0; d < 3; ++d) {
        SMESH_ASSERT_ARRAY_APPROX_EQ(a->n_nodes(), a->points()->data()[d], b->points()->data()[d], 1e-12);
    }
    return SMESH_TEST_SUCCESS;
}

static std::vector<std::pair<element_idx_t, i16>> sorted_sides(const std::shared_ptr<Mesh> &mesh) {
    std::vector<std::pair<element_idx_t, i16>> keys;
    for (const auto &named : mesh->sidesets()) {
        const auto ss = named.second;
        for (ptrdiff_t i = 0; i < ss->size(); ++i) {
            keys.push_back({ss->parent()->data()[i], ss->lfi()->data()[i]});
        }
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

static int test_bey_parity() {
    auto                    mesh = unit_tet();
    const element_idx_t     mark = 0;
    auto                    adapt = refine_adaptive(mesh, &mark, 1, 1);
    auto                    uniform = refine(mesh, 1);
    SMESH_TEST_EQ(same_mesh(adapt, uniform), SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(adapt->block(0)->name() == "tet");
    SMESH_TEST_EQ(adapt->geom_map(0), mesh->geom_map(0));
    SMESH_TEST_EQ(check_conforming(adapt, 1.0), SMESH_TEST_SUCCESS);

    auto cube = Mesh::create_tet4_cube(Communicator::self(), 1, 1, 1);
    std::vector<element_idx_t> all(static_cast<size_t>(cube->n_elements(0)));
    for (ptrdiff_t i = 0; i < cube->n_elements(0); ++i) {
        all[(size_t)i] = (element_idx_t)i;
    }
    auto cube_adapt   = refine_adaptive(cube, all.data(), cube->n_elements(0), 1);
    auto cube_uniform = refine(cube, 1);
    SMESH_TEST_EQ(cube_adapt->n_elements(0), cube_uniform->n_elements(0));
    SMESH_TEST_EQ(cube_adapt->n_nodes(), cube_uniform->n_nodes());
    SMESH_TEST_EQ(same_mesh(cube_adapt, cube_uniform), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_face_neighbor() {
    auto                mesh = face_pair();
    const element_idx_t mark = 0;
    auto                fine = refine_adaptive(mesh, &mark, 1, 1);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->n_elements(0), static_cast<ptrdiff_t>(12));
    SMESH_TEST_EQ(check_conforming(fine, 2.0), SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(count_shared_faces(fine) >= 4);
    return SMESH_TEST_SUCCESS;
}

static int test_edge_neighbor() {
    auto                mesh = edge_pair();
    const element_idx_t mark = 0;
    auto                fine = refine_adaptive(mesh, &mark, 1, 1);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->n_elements(0), static_cast<ptrdiff_t>(10));
    SMESH_TEST_EQ(check_conforming(fine, 6.0), SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(count_shared_faces(fine) >= 8);
    return SMESH_TEST_SUCCESS;
}

static int test_level_two() {
    auto                mesh = face_pair();
    const element_idx_t mark = 0;
    auto                fine = refine_adaptive(mesh, &mark, 1, 2);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_ASSERT(fine->n_elements(0) > 12);
    SMESH_TEST_EQ(check_conforming(fine, 2.0), SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(count_shared_faces(fine) > 0);
    const ptrdiff_t n  = fine->n_nodes();
    const geom_t   *px = fine->points()->data()[0];
    const geom_t   *py = fine->points()->data()[1];
    const geom_t   *pz = fine->points()->data()[2];
    int             saw_quarter = 0;
    for (ptrdiff_t i = 0; i < n; ++i) {
        const double xs[3] = {(double)px[i], (double)py[i], (double)pz[i]};
        for (int d = 0; d < 3; ++d) {
            const double q = xs[d] * 4.0;
            SMESH_TEST_APPROXEQ(q, std::round(q), 1e-6);
        }
        if (std::fabs((double)px[i] - 0.25) < 1e-8 && std::fabs((double)py[i]) < 1e-8 &&
            std::fabs((double)pz[i]) < 1e-8) {
            saw_quarter = 1;
        }
    }
    SMESH_TEST_ASSERT(saw_quarter == 1);
    return SMESH_TEST_SUCCESS;
}

static int test_sets() {
    auto                mesh = face_pair();
    const element_idx_t mark = 0;
    auto                parent = create_host_buffer<element_idx_t>(1);
    auto                lfi    = create_host_buffer<i16>(1);
    parent->data()[0]          = 1;
    lfi->data()[0]             = 3;
    mesh->add_sideset("shared", Sideset::create(mesh->comm(), parent, lfi, 0));
    auto nodes = create_host_buffer<idx_t>(2);
    nodes->data()[0] = 0;
    nodes->data()[1] = 1;
    mesh->add_nodeset("edge", Nodeset::create(mesh->comm(), nodes));

    auto fine = refine_adaptive(mesh, &mark, 1, 1);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->sidesets().size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(fine->sidesets()[0].second->size(), static_cast<ptrdiff_t>(4));
    const auto ss = fine->sidesets()[0].second;
    const idx_t *e0 = fine->elements(0)->data()[0];
    const idx_t *e1 = fine->elements(0)->data()[1];
    const idx_t *e2 = fine->elements(0)->data()[2];
    const idx_t *e3 = fine->elements(0)->data()[3];
    const geom_t *z = fine->points()->data()[2];
    for (ptrdiff_t i = 0; i < ss->size(); ++i) {
        const element_idx_t pe = ss->parent()->data()[i];
        const int           lf = ss->lfi()->data()[i];
        SMESH_TEST_ASSERT(lf >= 0 && lf < 4);
        const idx_t v[4] = {e0[pe], e1[pe], e2[pe], e3[pe]};
        for (int k = 0; k < 3; ++k) {
            SMESH_TEST_APPROXEQ((double)z[v[k_tet_face[lf][k]]], 0.0, 1e-8);
        }
    }

    const auto ns = fine->nodesets()[0].second;
    int        saw_mid = 0;
    int        saw0    = 0;
    int        saw1    = 0;
    const geom_t *x = fine->points()->data()[0];
    const geom_t *y = fine->points()->data()[1];
    for (ptrdiff_t i = 0; i < ns->size(); ++i) {
        const idx_t id = ns->nodes()->data()[i];
        if (id == 0) {
            saw0 = 1;
        }
        if (id == 1) {
            saw1 = 1;
        }
        if (std::fabs((double)x[id] - 0.5) < 1e-8 && std::fabs((double)y[id]) < 1e-8 &&
            std::fabs((double)z[id]) < 1e-8) {
            saw_mid = 1;
        }
    }
    SMESH_TEST_ASSERT(saw0 == 1 && saw1 == 1 && saw_mid == 1);

    auto                one = unit_tet();
    auto                ep  = create_host_buffer<element_idx_t>(1);
    auto                el  = create_host_buffer<i16>(1);
    ep->data()[0]           = 0;
    el->data()[0]           = 0;
    one->add_edgeset("e0", Edgeset::create(one->comm(), ep, el, 0));
    auto                sp = create_host_buffer<element_idx_t>(1);
    auto                sl = create_host_buffer<i16>(1);
    sp->data()[0]          = 0;
    sl->data()[0]          = 0;
    one->add_sideset("f0", Sideset::create(one->comm(), sp, sl, 0));
    auto nodes1 = create_host_buffer<idx_t>(2);
    nodes1->data()[0] = 0;
    nodes1->data()[1] = 1;
    one->add_nodeset("n", Nodeset::create(one->comm(), nodes1));

    const element_idx_t one_mark = 0;
    auto                adapt_one = refine_adaptive(one, &one_mark, 1, 1);
    auto                uniform   = refine(one, 1);
    SMESH_TEST_ASSERT(adapt_one != nullptr && uniform != nullptr);
    SMESH_TEST_EQ(adapt_one->edgesets()[0].second->size(), static_cast<ptrdiff_t>(2));
    const auto ak = sorted_sides(adapt_one);
    const auto uk = sorted_sides(uniform);
    SMESH_TEST_EQ(ak.size(), uk.size());
    for (size_t i = 0; i < ak.size(); ++i) {
        SMESH_TEST_EQ(ak[i].first, uk[i].first);
        SMESH_TEST_EQ(ak[i].second, uk[i].second);
    }
    std::vector<idx_t> an(adapt_one->nodesets()[0].second->nodes()->data(),
                          adapt_one->nodesets()[0].second->nodes()->data() + adapt_one->nodesets()[0].second->size());
    std::vector<idx_t> un(uniform->nodesets()[0].second->nodes()->data(),
                          uniform->nodesets()[0].second->nodes()->data() + uniform->nodesets()[0].second->size());
    std::sort(an.begin(), an.end());
    std::sort(un.begin(), un.end());
    SMESH_TEST_EQ(an.size(), un.size());
    for (size_t i = 0; i < an.size(); ++i) {
        SMESH_TEST_EQ(an[i], un[i]);
    }
    return SMESH_TEST_SUCCESS;
}

static ptrdiff_t count_centroid_x_above(const std::shared_ptr<Mesh> &mesh, double thresh) {
    const idx_t  *e0 = mesh->elements(0)->data()[0];
    const idx_t  *e1 = mesh->elements(0)->data()[1];
    const idx_t  *e2 = mesh->elements(0)->data()[2];
    const idx_t  *e3 = mesh->elements(0)->data()[3];
    const geom_t *x  = mesh->points()->data()[0];
    ptrdiff_t     n  = 0;
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        const double cx = 0.25 * ((double)x[e0[e]] + (double)x[e1[e]] + (double)x[e2[e]] + (double)x[e3[e]]);
        if (cx > thresh) {
            ++n;
        }
    }
    return n;
}

/// Level 2 grades face neighbors to depth 1 and stops there. The far half of a
/// neighbor hex must stay at the coarse resolution.
static int test_level_two_stays_local() {
    auto hex = Mesh::create_hex8_cube(Communicator::self(), 2, 1, 1);
    auto tet = convert_to(TET4, hex);
    SMESH_TEST_ASSERT(tet != nullptr);
    const ptrdiff_t n_coarse = tet->n_elements(0);
    const idx_t    *e0       = tet->elements(0)->data()[0];
    const idx_t    *e1       = tet->elements(0)->data()[1];
    const idx_t    *e2       = tet->elements(0)->data()[2];
    const idx_t    *e3       = tet->elements(0)->data()[3];
    const geom_t   *x        = tet->points()->data()[0];
    std::vector<element_idx_t> marks;
    for (ptrdiff_t e = 0; e < n_coarse; ++e) {
        const double cx = 0.25 * ((double)x[e0[e]] + (double)x[e1[e]] + (double)x[e2[e]] + (double)x[e3[e]]);
        if (cx < 0.5) {
            marks.push_back((element_idx_t)e);
        }
    }
    SMESH_TEST_ASSERT(!marks.empty());
    const ptrdiff_t far_coarse = count_centroid_x_above(tet, 0.8);
    SMESH_TEST_ASSERT(far_coarse > 0);

    auto fine1 = refine_adaptive(tet, marks.data(), (ptrdiff_t)marks.size(), 1);
    auto fine2 = refine_adaptive(tet, marks.data(), (ptrdiff_t)marks.size(), 2);
    SMESH_TEST_ASSERT(fine1 != nullptr && fine2 != nullptr);
    SMESH_TEST_ASSERT(fine2->n_elements(0) > fine1->n_elements(0));
    SMESH_TEST_EQ(count_centroid_x_above(fine2, 0.8), far_coarse);
    SMESH_TEST_EQ(check_conforming(fine2, 6.0, 1e-5), SMESH_TEST_SUCCESS);

    auto one = convert_to(TET4, Mesh::create_hex8_cube(Communicator::self(), 1, 1, 1));
    SMESH_TEST_ASSERT(one != nullptr);
    const element_idx_t one_mark = 0;
    auto one_l1   = refine_adaptive(one, &one_mark, 1, 1);
    auto one_fine = refine_adaptive(one, &one_mark, 1, 2);
    SMESH_TEST_ASSERT(one_fine != nullptr && one_l1 != nullptr);
    // Level 1 does not grade. Level 2 refines the marked tet to depth 2 (64) and
    // its three face neighbors to depth 1, then closes the remaining edge
    // neighbors. That is 205 tets. Promoting every edge neighbor is 251.
    SMESH_TEST_EQ(one_l1->n_elements(0), static_cast<ptrdiff_t>(47));
    SMESH_TEST_EQ(one_fine->n_elements(0), static_cast<ptrdiff_t>(205));
    SMESH_TEST_EQ(check_conforming(one_fine, 6.0, 1e-5), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

/// Two marked tets bisect two edges of a third. That mask is not one edge or one
/// face. The third tet must close without an 8-split of its other edges.
static int test_irregular_closure() {
    const idx_t  conn[12] = {0, 1, 0, 1, 2, 1, 5, 6, 2, 4, 7, 3};
    const geom_t xyz[24]  = {0, 1, 0, 0, 0.5f, 0.5f, 1.2f, 0.8f, 0, 0, 1, 0, -1, -0.2f, 1.2f, 0.8f,
                             0, 0, 0, 1, 0,    0.4f,  0.2f, -0.8f};
    auto         mesh = make_tets(conn, 3, xyz, 8);
    const element_idx_t marks[2] = {0, 1};
    auto                parent   = create_host_buffer<element_idx_t>(1);
    auto                lfi      = create_host_buffer<i16>(1);
    parent->data()[0]            = 2;
    lfi->data()[0]               = 3;
    mesh->add_sideset("face", Sideset::create(mesh->comm(), parent, lfi, 0));

    const idx_t  *e0 = mesh->elements(0)->data()[0];
    const idx_t  *e1 = mesh->elements(0)->data()[1];
    const idx_t  *e2 = mesh->elements(0)->data()[2];
    const idx_t  *e3 = mesh->elements(0)->data()[3];
    const geom_t *x  = mesh->points()->data()[0];
    const geom_t *y  = mesh->points()->data()[1];
    const geom_t *z  = mesh->points()->data()[2];
    double        volume = 0;
    for (ptrdiff_t e = 0; e < 3; ++e) {
        volume += tet_volume(x, y, z, e0[e], e1[e], e2[e], e3[e]);
    }

    auto fine = refine_adaptive(mesh, marks, 2, 1);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->sidesets()[0].second->size(), static_cast<ptrdiff_t>(3));
    SMESH_TEST_EQ(check_conforming(fine, volume, 1e-5), SMESH_TEST_SUCCESS);

    auto deeper = refine_adaptive(mesh, marks, 2, 2);
    SMESH_TEST_ASSERT(deeper != nullptr);
    SMESH_TEST_EQ(check_conforming(deeper, volume, 1e-5), SMESH_TEST_SUCCESS);
    SMESH_TEST_ASSERT(deeper->n_elements(0) > fine->n_elements(0));
    return SMESH_TEST_SUCCESS;
}

static int test_reject_and_clone() {
    auto                mesh = unit_tet();
    const element_idx_t mark = 0;
    auto                parent = create_host_buffer<element_idx_t>(1);
    auto                lfi    = create_host_buffer<i16>(1);
    parent->data()[0]          = 0;
    lfi->data()[0]             = 1;
    mesh->add_sideset("f", Sideset::create(mesh->comm(), parent, lfi, 0));

    auto copy = refine_adaptive(mesh, &mark, 1, 0);
    SMESH_TEST_ASSERT(copy != nullptr);
    SMESH_TEST_EQ(copy->n_elements(0), mesh->n_elements(0));
    SMESH_TEST_EQ(copy->sidesets().size(), static_cast<size_t>(1));
    SMESH_TEST_EQ(copy->sidesets()[0].second->size(), static_cast<ptrdiff_t>(1));

    auto empty = refine_adaptive(mesh, nullptr, 0, 3);
    SMESH_TEST_ASSERT(empty != nullptr);
    SMESH_TEST_EQ(empty->n_elements(0), static_cast<ptrdiff_t>(1));
    SMESH_TEST_EQ(empty->sidesets()[0].second->lfi()->data()[0], static_cast<i16>(1));

    SMESH_TEST_ASSERT(refine_adaptive(mesh, &mark, 1, -1) == nullptr);
    const element_idx_t bad = 4;
    SMESH_TEST_ASSERT(refine_adaptive(mesh, &bad, 1, 1) == nullptr);

    auto hex = Mesh::create_hex8_cube(Communicator::self(), 1, 1, 1);
    SMESH_TEST_ASSERT(refine_adaptive(hex, &mark, 1, 1) == nullptr);

    mesh->set_element_type(0, TET10);
    SMESH_TEST_ASSERT(refine_adaptive(mesh, &mark, 1, 1) == nullptr);
    mesh->set_element_type(0, PROTEUS_TET4);
    SMESH_TEST_ASSERT(refine_adaptive(mesh, &mark, 1, 1) == nullptr);
    mesh->set_element_type(0, TET4);

    auto extra = create_host_buffer<idx_t>(4, static_cast<size_t>(1));
    extra->data()[0][0] = 0;
    extra->data()[1][0] = 1;
    extra->data()[2][0] = 2;
    extra->data()[3][0] = 3;
    mesh->add_block("other", TET4, extra);
    mesh->set_element_type(1, HEX8);
    SMESH_TEST_ASSERT(refine_adaptive(mesh, &mark, 1, 1) == nullptr);
    return SMESH_TEST_SUCCESS;
}

static int test_multiblock_sets() {
    auto hex = Mesh::create_hex8_bidomain_cube(Communicator::self(), 2, 1, 1);
    SMESH_TEST_ASSERT(hex != nullptr);
    auto face = Sideset::create_from_plane(hex, 1, 0, 0, 1, static_cast<geom_t>(1e-5), {"right"});
    SMESH_TEST_ASSERT(face.size() == 1);
    hex->add_sideset("xmax", face[0]);

    auto tet = convert_to(TET4, hex);
    SMESH_TEST_ASSERT(tet != nullptr);
    SMESH_TEST_EQ(tet->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(tet->block(0)->name(), std::string("left"));
    SMESH_TEST_EQ(tet->block(1)->name(), std::string("right"));
    SMESH_TEST_EQ(tet->n_elements(0), hex->n_elements(0) * 24);
    SMESH_TEST_EQ(tet->n_elements(1), hex->n_elements(1) * 24);
    SMESH_TEST_EQ(tet->sidesets("xmax").size(), static_cast<size_t>(1));
    const auto coarse_ss = tet->sidesets("xmax")[0];
    SMESH_TEST_EQ(coarse_ss->block_id(), static_cast<block_idx_t>(1));
    SMESH_TEST_EQ(coarse_ss->size(), face[0]->size() * 4);

    const ptrdiff_t n0 = tet->n_elements(0);
    const ptrdiff_t n1 = tet->n_elements(1);
    std::vector<element_idx_t> marks(static_cast<size_t>(n1));
    for (ptrdiff_t e = 0; e < n1; ++e) {
        marks[static_cast<size_t>(e)] = static_cast<element_idx_t>(n0 + e);
    }

    auto copy = refine_adaptive(tet, marks.data(), n1, 0);
    SMESH_TEST_ASSERT(copy != nullptr);
    SMESH_TEST_EQ(copy->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(copy->block(1)->name(), std::string("right"));
    SMESH_TEST_EQ(copy->sidesets("xmax")[0]->size(), coarse_ss->size());
    SMESH_TEST_EQ(copy->sidesets("xmax")[0]->block_id(), static_cast<block_idx_t>(1));

    auto fine = refine_adaptive(tet, marks.data(), n1, 1);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->n_blocks(), static_cast<size_t>(2));
    SMESH_TEST_EQ(fine->block(0)->name(), std::string("left"));
    SMESH_TEST_EQ(fine->block(1)->name(), std::string("right"));
    SMESH_TEST_ASSERT(fine->n_elements(1) > n1);
    SMESH_TEST_EQ(fine->sidesets().size(), static_cast<size_t>(1));
    const auto ss = fine->sidesets("xmax")[0];
    SMESH_TEST_EQ(ss->block_id(), static_cast<block_idx_t>(1));
    SMESH_TEST_EQ(ss->size(), coarse_ss->size() * 4);

    const element_idx_t *parent = ss->parent()->data();
    const i16           *lfi    = ss->lfi()->data();
    const idx_t         *e0     = fine->elements(1)->data()[0];
    const idx_t         *e1     = fine->elements(1)->data()[1];
    const idx_t         *e2     = fine->elements(1)->data()[2];
    const idx_t         *e3     = fine->elements(1)->data()[3];
    const geom_t        *x      = fine->points()->data()[0];
    for (ptrdiff_t i = 0; i < ss->size(); ++i) {
        SMESH_TEST_ASSERT(parent[i] >= 0 && static_cast<ptrdiff_t>(parent[i]) < fine->n_elements(1));
        SMESH_TEST_ASSERT(lfi[i] >= 0 && lfi[i] < 4);
        const idx_t v[4] = {e0[parent[i]], e1[parent[i]], e2[parent[i]], e3[parent[i]]};
        for (int k = 0; k < 3; ++k) {
            SMESH_TEST_APPROXEQ(static_cast<double>(x[v[k_tet_face[lfi[i]][k]]]), 1.0, 1e-5);
        }
    }
    SMESH_TEST_EQ(check_conforming(fine, 6.0, 1e-5), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static int test_graded_depths() {
    const idx_t  conn[8] = {0, 4, 1, 5, 2, 6, 3, 7};
    const geom_t xyz[24] = {0, 1, 0, 0, 3, 4, 3, 3, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1};
    auto mesh = make_tets(conn, 2, xyz, 8);
    const element_idx_t marked[2] = {0, 1};
    const int           depth[2]  = {2, 1};
    auto fine = refine_adaptive(mesh, marked, depth, 2);
    SMESH_TEST_ASSERT(fine != nullptr);
    SMESH_TEST_EQ(fine->n_elements(0), static_cast<ptrdiff_t>(72));

    const idx_t *e0 = fine->elements(0)->data()[0];
    const idx_t *e1 = fine->elements(0)->data()[1];
    const idx_t *e2 = fine->elements(0)->data()[2];
    const idx_t *e3 = fine->elements(0)->data()[3];
    const geom_t *x = fine->points()->data()[0];
    std::vector<char> owner(static_cast<size_t>(fine->n_nodes()), 0);
    ptrdiff_t n_deep = 0;
    ptrdiff_t n_shallow = 0;
    for (ptrdiff_t e = 0; e < fine->n_elements(0); ++e) {
        const idx_t v[4] = {e0[e], e1[e], e2[e], e3[e]};
        const geom_t cx = (x[v[0]] + x[v[1]] + x[v[2]] + x[v[3]]) * static_cast<geom_t>(0.25);
        const char tag = cx < static_cast<geom_t>(2) ? 1 : 2;
        if (tag == 1) {
            ++n_deep;
        } else {
            ++n_shallow;
        }
        for (int k = 0; k < 4; ++k) {
            char &slot = owner[static_cast<size_t>(v[k])];
            SMESH_TEST_ASSERT(slot == 0 || slot == tag);
            slot = tag;
        }
    }
    SMESH_TEST_EQ(n_deep, static_cast<ptrdiff_t>(64));
    SMESH_TEST_EQ(n_shallow, static_cast<ptrdiff_t>(8));
    SMESH_TEST_EQ(check_conforming(fine, 2.0, 1e-5), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

static double orient_tet(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy,
                         double cz, double dx, double dy, double dz) {
    const double ux = bx - ax;
    const double uy = by - ay;
    const double uz = bz - az;
    const double vx = cx - ax;
    const double vy = cy - ay;
    const double vz = cz - az;
    const double wx = dx - ax;
    const double wy = dy - ay;
    const double wz = dz - az;
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx);
}

static ptrdiff_t count_centroids_in_tet(const std::shared_ptr<Mesh> &mesh, const double p[4][3]) {
    const idx_t  *e0 = mesh->elements(0)->data()[0];
    const idx_t  *e1 = mesh->elements(0)->data()[1];
    const idx_t  *e2 = mesh->elements(0)->data()[2];
    const idx_t  *e3 = mesh->elements(0)->data()[3];
    const geom_t *x  = mesh->points()->data()[0];
    const geom_t *y  = mesh->points()->data()[1];
    const geom_t *z  = mesh->points()->data()[2];
    const double  parent =
            orient_tet(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2], p[2][0], p[2][1], p[2][2], p[3][0],
                       p[3][1], p[3][2]);
    ptrdiff_t n = 0;
    for (ptrdiff_t e = 0; e < mesh->n_elements(0); ++e) {
        const double cx = 0.25 * ((double)x[e0[e]] + (double)x[e1[e]] + (double)x[e2[e]] + (double)x[e3[e]]);
        const double cy = 0.25 * ((double)y[e0[e]] + (double)y[e1[e]] + (double)y[e2[e]] + (double)y[e3[e]]);
        const double cz = 0.25 * ((double)z[e0[e]] + (double)z[e1[e]] + (double)z[e2[e]] + (double)z[e3[e]]);
        const double s0 = orient_tet(cx, cy, cz, p[1][0], p[1][1], p[1][2], p[2][0], p[2][1], p[2][2], p[3][0],
                                     p[3][1], p[3][2]);
        const double s1 = orient_tet(p[0][0], p[0][1], p[0][2], cx, cy, cz, p[2][0], p[2][1], p[2][2], p[3][0],
                                     p[3][1], p[3][2]);
        const double s2 = orient_tet(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2], cx, cy, cz, p[3][0],
                                     p[3][1], p[3][2]);
        const double s3 = orient_tet(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2], p[2][0], p[2][1], p[2][2],
                                     cx, cy, cz);
        const double tol = std::fabs(parent) * 1e-6;
        if (s0 * parent >= -tol && s1 * parent >= -tol && s2 * parent >= -tol && s3 * parent >= -tol) {
            ++n;
        }
    }
    return n;
}

/// A at depth 2 forces face neighbor B to depth 1. C shares a face with B and
/// none with A, so depth 2 leaves C unmarked and depth 3 raises C to 1.
static int test_grade_hops() {
    const idx_t  conn[12] = {0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5};
    const geom_t xyz[18]  = {0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 1};
    auto         mesh = make_tets(conn, 3, xyz, 6);
    const double tet_a[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const double tet_b[4][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}};
    const double tet_c[4][3] = {{0, 1, 0}, {0, 0, 1}, {1, 1, 1}, {0, 2, 1}};

    const element_idx_t mark = 0;
    auto                level2 = refine_adaptive(mesh, &mark, 1, 2);
    SMESH_TEST_ASSERT(level2 != nullptr);
    const ptrdiff_t a2 = count_centroids_in_tet(level2, tet_a);
    const ptrdiff_t b2 = count_centroids_in_tet(level2, tet_b);
    const ptrdiff_t c2 = count_centroids_in_tet(level2, tet_c);
    SMESH_TEST_EQ(a2, static_cast<ptrdiff_t>(64));
    SMESH_TEST_ASSERT(b2 > 8 && b2 < 64);
    SMESH_TEST_ASSERT(c2 > 1 && c2 < 8);
    SMESH_TEST_EQ(check_conforming(level2, 5.0, 1e-5), SMESH_TEST_SUCCESS);

    auto level3 = refine_adaptive(mesh, &mark, 1, 3);
    SMESH_TEST_ASSERT(level3 != nullptr);
    const ptrdiff_t c3 = count_centroids_in_tet(level3, tet_c);
    SMESH_TEST_ASSERT(c3 >= 8 && c3 < 64);
    SMESH_TEST_EQ(check_conforming(level3, 5.0, 1e-5), SMESH_TEST_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_bey_parity);
    SMESH_RUN_TEST(test_face_neighbor);
    SMESH_RUN_TEST(test_edge_neighbor);
    SMESH_RUN_TEST(test_level_two);
    SMESH_RUN_TEST(test_level_two_stays_local);
    SMESH_RUN_TEST(test_irregular_closure);
    SMESH_RUN_TEST(test_sets);
    SMESH_RUN_TEST(test_reject_and_clone);
    SMESH_RUN_TEST(test_multiblock_sets);
    SMESH_RUN_TEST(test_graded_depths);
    SMESH_RUN_TEST(test_grade_hops);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
