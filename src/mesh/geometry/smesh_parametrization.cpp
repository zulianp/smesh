#include "smesh_parametrization.hpp"

#include "smesh_mesh.hpp"
#include "smesh_nodeset.hpp"

#include <cmath>
#include <cstddef>

namespace smesh {

namespace {

int unit3(geom_t *x, geom_t *y, geom_t *z) {
    const geom_t n2 = (*x) * (*x) + (*y) * (*y) + (*z) * (*z);
    if (!(n2 > static_cast<geom_t>(0))) {
        return SMESH_FAILURE;
    }
    const geom_t inv = static_cast<geom_t>(1) / std::sqrt(n2);
    *x *= inv;
    *y *= inv;
    *z *= inv;
    return SMESH_SUCCESS;
}

void cross3(geom_t ax,
            geom_t ay,
            geom_t az,
            geom_t bx,
            geom_t by,
            geom_t bz,
            geom_t *cx,
            geom_t *cy,
            geom_t *cz) {
    *cx = ay * bz - az * by;
    *cy = az * bx - ax * bz;
    *cz = ax * by - ay * bx;
}

int require_nodeset(const std::shared_ptr<Nodeset> &ns, const char *who) {
    if (!ns || !ns->nodes()) {
        SMESH_ERROR("%s: missing nodeset\n", who);
        return SMESH_FAILURE;
    }
    return SMESH_SUCCESS;
}

int gather_scatter_prepare(Mesh                             &mesh,
                           const std::shared_ptr<Nodeset>   &ns,
                           const char                       *who,
                           ptrdiff_t                        *n_out,
                           const idx_t                     **ids_out,
                           geom_t                          **px,
                           geom_t                          **py,
                           geom_t                          **pz) {
    if (require_nodeset(ns, who) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    if (mesh.spatial_dimension() < 3 || !mesh.points() || !mesh.points()->data()) {
        SMESH_ERROR("%s: mesh must have 3D points\n", who);
        return SMESH_FAILURE;
    }
    geom_t **const pts = mesh.points()->data();
    const ptrdiff_t n_nodes = mesh.n_nodes();
    const ptrdiff_t n       = ns->size();
    const idx_t    *ids     = n > 0 ? ns->nodes()->data() : nullptr;
    for (ptrdiff_t i = 0; i < n; ++i) {
        const idx_t id = ids[i];
        if (id < 0 || static_cast<ptrdiff_t>(id) >= n_nodes) {
            SMESH_ERROR("%s: node id out of range\n", who);
            return SMESH_FAILURE;
        }
    }
    *n_out   = n;
    *ids_out = ids;
    *px      = pts[0];
    *py      = pts[1];
    *pz      = pts[2];
    return SMESH_SUCCESS;
}

void gather_xyz(const ptrdiff_t n,
                const idx_t *const SMESH_RESTRICT ids,
                const geom_t *const SMESH_RESTRICT px,
                const geom_t *const SMESH_RESTRICT py,
                const geom_t *const SMESH_RESTRICT pz,
                geom_t *const SMESH_RESTRICT x,
                geom_t *const SMESH_RESTRICT y,
                geom_t *const SMESH_RESTRICT z) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        const idx_t id = ids[i];
        x[i]           = px[id];
        y[i]           = py[id];
        z[i]           = pz[id];
    }
}

void scatter_xyz(const ptrdiff_t n,
                 const idx_t *const SMESH_RESTRICT ids,
                 const geom_t *const SMESH_RESTRICT x,
                 const geom_t *const SMESH_RESTRICT y,
                 const geom_t *const SMESH_RESTRICT z,
                 geom_t *const SMESH_RESTRICT px,
                 geom_t *const SMESH_RESTRICT py,
                 geom_t *const SMESH_RESTRICT pz) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        const idx_t id = ids[i];
        px[id]         = x[i];
        py[id]         = y[i];
        pz[id]         = z[i];
    }
}

void project_sphere(const ptrdiff_t n,
                    geom_t *const SMESH_RESTRICT x,
                    geom_t *const SMESH_RESTRICT y,
                    geom_t *const SMESH_RESTRICT z,
                    const geom_t                 cx,
                    const geom_t                 cy,
                    const geom_t                 cz,
                    const geom_t                 radius) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        const geom_t dx  = x[i] - cx;
        const geom_t dy  = y[i] - cy;
        const geom_t dz  = z[i] - cz;
        const geom_t len2 = dx * dx + dy * dy + dz * dz;
        if (!(len2 > static_cast<geom_t>(0))) {
            x[i] = cx + radius;
            y[i] = cy;
            z[i] = cz;
            continue;
        }
        const geom_t s = radius / std::sqrt(len2);
        x[i]           = cx + s * dx;
        y[i]           = cy + s * dy;
        z[i]           = cz + s * dz;
    }
}

void project_circle(const ptrdiff_t n,
                    geom_t *const SMESH_RESTRICT x,
                    geom_t *const SMESH_RESTRICT y,
                    geom_t *const SMESH_RESTRICT z,
                    const geom_t                 cx,
                    const geom_t                 cy,
                    const geom_t                 cz,
                    const geom_t                 ax,
                    const geom_t                 ay,
                    const geom_t                 az,
                    const geom_t                 radius) {
    geom_t fx = static_cast<geom_t>(1);
    geom_t fy = static_cast<geom_t>(0);
    geom_t fz = static_cast<geom_t>(0);
    if (std::fabs(ax) > static_cast<geom_t>(0.9)) {
        fx = static_cast<geom_t>(0);
        fy = static_cast<geom_t>(1);
        fz = static_cast<geom_t>(0);
    }
    geom_t ux, uy, uz;
    cross3(ax, ay, az, fx, fy, fz, &ux, &uy, &uz);
    if (unit3(&ux, &uy, &uz) != SMESH_SUCCESS) {
        ux = static_cast<geom_t>(0);
        uy = static_cast<geom_t>(1);
        uz = static_cast<geom_t>(0);
        if (unit3(&ux, &uy, &uz) != SMESH_SUCCESS) {
            ux = static_cast<geom_t>(1);
            uy = static_cast<geom_t>(0);
            uz = static_cast<geom_t>(0);
        }
    }

    for (ptrdiff_t i = 0; i < n; ++i) {
        const geom_t dx = x[i] - cx;
        const geom_t dy = y[i] - cy;
        const geom_t dz = z[i] - cz;
        const geom_t a  = dx * ax + dy * ay + dz * az;
        geom_t       px = dx - a * ax;
        geom_t       py = dy - a * ay;
        geom_t       pz = dz - a * az;
        const geom_t len2 = px * px + py * py + pz * pz;
        if (!(len2 > static_cast<geom_t>(0))) {
            x[i] = cx + radius * ux;
            y[i] = cy + radius * uy;
            z[i] = cz + radius * uz;
            continue;
        }
        const geom_t s = radius / std::sqrt(len2);
        x[i]           = cx + s * px;
        y[i]           = cy + s * py;
        z[i]           = cz + s * pz;
    }
}

geom_t eval_poly(const int degree, const geom_t *const SMESH_RESTRICT c, const geom_t xi, const geom_t eta) {
    geom_t s = static_cast<geom_t>(0);
    int    k = 0;
    geom_t xi_i = static_cast<geom_t>(1);
    for (int i = 0; i <= degree; ++i) {
        geom_t eta_j = static_cast<geom_t>(1);
        for (int j = 0; j <= degree - i; ++j) {
            s += c[k++] * xi_i * eta_j;
            eta_j *= eta;
        }
        xi_i *= xi;
    }
    return s;
}

void project_polynomial_surface(const ptrdiff_t n,
                                geom_t *const SMESH_RESTRICT x,
                                geom_t *const SMESH_RESTRICT y,
                                geom_t *const SMESH_RESTRICT z,
                                const geom_t                 ox,
                                const geom_t                 oy,
                                const geom_t                 oz,
                                const geom_t                 e0x,
                                const geom_t                 e0y,
                                const geom_t                 e0z,
                                const geom_t                 e1x,
                                const geom_t                 e1y,
                                const geom_t                 e1z,
                                const geom_t                 e2x,
                                const geom_t                 e2y,
                                const geom_t                 e2z,
                                const int                    degree,
                                const geom_t *const SMESH_RESTRICT coeffs) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        const geom_t rx = x[i] - ox;
        const geom_t ry = y[i] - oy;
        const geom_t rz = z[i] - oz;
        const geom_t xi  = rx * e0x + ry * e0y + rz * e0z;
        const geom_t eta = rx * e1x + ry * e1y + rz * e1z;
        const geom_t zeta = eval_poly(degree, coeffs, xi, eta);
        x[i] = ox + xi * e0x + eta * e1x + zeta * e2x;
        y[i] = oy + xi * e0y + eta * e1y + zeta * e2y;
        z[i] = oz + xi * e0z + eta * e1z + zeta * e2z;
    }
}

int apply_gathered(Mesh                           &mesh,
                   const std::shared_ptr<Nodeset> &ns,
                   const char                     *who,
                   void (*kernel)(ptrdiff_t, geom_t *, geom_t *, geom_t *, const void *),
                   const void                     *ctx) {
    ptrdiff_t    n   = 0;
    const idx_t *ids = nullptr;
    geom_t      *px = nullptr, *py = nullptr, *pz = nullptr;
    if (gather_scatter_prepare(mesh, ns, who, &n, &ids, &px, &py, &pz) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    if (n <= 0) {
        return SMESH_SUCCESS;
    }
    auto buf = create_host_buffer<geom_t>(3, static_cast<size_t>(n));
    geom_t *const x = buf->data()[0];
    geom_t *const y = buf->data()[1];
    geom_t *const z = buf->data()[2];
    gather_xyz(n, ids, px, py, pz, x, y, z);
    kernel(n, x, y, z, ctx);
    scatter_xyz(n, ids, x, y, z, px, py, pz);
    return SMESH_SUCCESS;
}

}  // namespace

CircleParametrization::CircleParametrization(const std::shared_ptr<Nodeset> &nodeset,
                                             geom_t                          cx,
                                             geom_t                          cy,
                                             geom_t                          cz,
                                             geom_t                          ax,
                                             geom_t                          ay,
                                             geom_t                          az,
                                             geom_t                          radius)
    : nodeset_(nodeset), cx_(cx), cy_(cy), cz_(cz), ax_(ax), ay_(ay), az_(az), radius_(radius) {
    if (require_nodeset(nodeset_, "CircleParametrization") != SMESH_SUCCESS) {
        return;
    }
    if (!(radius_ > static_cast<geom_t>(0))) {
        SMESH_ERROR("CircleParametrization: radius must be positive\n");
    }
    if (unit3(&ax_, &ay_, &az_) != SMESH_SUCCESS) {
        SMESH_ERROR("CircleParametrization: axis must be non-zero\n");
    }
}

std::shared_ptr<CircleParametrization> CircleParametrization::create(const std::shared_ptr<Nodeset> &nodeset,
                                                                     geom_t                          cx,
                                                                     geom_t                          cy,
                                                                     geom_t                          cz,
                                                                     geom_t                          ax,
                                                                     geom_t                          ay,
                                                                     geom_t                          az,
                                                                     geom_t                          radius) {
    return std::make_shared<CircleParametrization>(nodeset, cx, cy, cz, ax, ay, az, radius);
}

int CircleParametrization::apply(Mesh &mesh) const {
    struct Ctx {
        geom_t cx, cy, cz, ax, ay, az, radius;
    } ctx{cx_, cy_, cz_, ax_, ay_, az_, radius_};
    return apply_gathered(
            mesh,
            nodeset_,
            "CircleParametrization::apply",
            [](ptrdiff_t n, geom_t *x, geom_t *y, geom_t *z, const void *p) {
                const auto *c = static_cast<const Ctx *>(p);
                project_circle(n, x, y, z, c->cx, c->cy, c->cz, c->ax, c->ay, c->az, c->radius);
            },
            &ctx);
}

SphereParametrization::SphereParametrization(const std::shared_ptr<Nodeset> &nodeset,
                                             geom_t                          cx,
                                             geom_t                          cy,
                                             geom_t                          cz,
                                             geom_t                          radius)
    : nodeset_(nodeset), cx_(cx), cy_(cy), cz_(cz), radius_(radius) {
    if (require_nodeset(nodeset_, "SphereParametrization") != SMESH_SUCCESS) {
        return;
    }
    if (!(radius_ > static_cast<geom_t>(0))) {
        SMESH_ERROR("SphereParametrization: radius must be positive\n");
    }
}

std::shared_ptr<SphereParametrization> SphereParametrization::create(const std::shared_ptr<Nodeset> &nodeset,
                                                                     geom_t                          cx,
                                                                     geom_t                          cy,
                                                                     geom_t                          cz,
                                                                     geom_t                          radius) {
    return std::make_shared<SphereParametrization>(nodeset, cx, cy, cz, radius);
}

int SphereParametrization::apply(Mesh &mesh) const {
    struct Ctx {
        geom_t cx, cy, cz, radius;
    } ctx{cx_, cy_, cz_, radius_};
    return apply_gathered(
            mesh,
            nodeset_,
            "SphereParametrization::apply",
            [](ptrdiff_t n, geom_t *x, geom_t *y, geom_t *z, const void *p) {
                const auto *c = static_cast<const Ctx *>(p);
                project_sphere(n, x, y, z, c->cx, c->cy, c->cz, c->radius);
            },
            &ctx);
}

PolynomialSurfaceParametrization::PolynomialSurfaceParametrization(const std::shared_ptr<Nodeset> &nodeset,
                                                                   geom_t                          ox,
                                                                   geom_t                          oy,
                                                                   geom_t                          oz,
                                                                   geom_t                          e0x,
                                                                   geom_t                          e0y,
                                                                   geom_t                          e0z,
                                                                   geom_t                          e1x,
                                                                   geom_t                          e1y,
                                                                   geom_t                          e1z,
                                                                   int                             degree,
                                                                   const geom_t                   *coeffs,
                                                                   ptrdiff_t                       n_coeffs)
    : nodeset_(nodeset), ox_(ox), oy_(oy), oz_(oz), e0x_(e0x), e0y_(e0y), e0z_(e0z), e1x_(e1x), e1y_(e1y), e1z_(e1z),
      e2x_(0), e2y_(0), e2z_(0), degree_(degree), coeffs_(nullptr) {
    if (require_nodeset(nodeset_, "PolynomialSurfaceParametrization") != SMESH_SUCCESS) {
        return;
    }
    if (degree_ < 0) {
        SMESH_ERROR("PolynomialSurfaceParametrization: degree must be >= 0\n");
    }
    const ptrdiff_t expect =
            (static_cast<ptrdiff_t>(degree_) + 1) * (static_cast<ptrdiff_t>(degree_) + 2) / 2;
    if (!coeffs || n_coeffs != expect) {
        SMESH_ERROR("PolynomialSurfaceParametrization: expected %td coefficients\n", (ptrdiff_t)expect);
    }
    if (unit3(&e0x_, &e0y_, &e0z_) != SMESH_SUCCESS) {
        SMESH_ERROR("PolynomialSurfaceParametrization: e0 must be non-zero\n");
    }
    const geom_t d = e1x_ * e0x_ + e1y_ * e0y_ + e1z_ * e0z_;
    e1x_ -= d * e0x_;
    e1y_ -= d * e0y_;
    e1z_ -= d * e0z_;
    if (unit3(&e1x_, &e1y_, &e1z_) != SMESH_SUCCESS) {
        SMESH_ERROR("PolynomialSurfaceParametrization: e0 and e1 must span a plane\n");
    }
    cross3(e0x_, e0y_, e0z_, e1x_, e1y_, e1z_, &e2x_, &e2y_, &e2z_);
    coeffs_ = create_host_buffer<geom_t>(static_cast<size_t>(expect));
    for (ptrdiff_t i = 0; i < expect; ++i) {
        coeffs_->data()[i] = coeffs[i];
    }
}

std::shared_ptr<PolynomialSurfaceParametrization> PolynomialSurfaceParametrization::create(
        const std::shared_ptr<Nodeset> &nodeset,
        geom_t                          ox,
        geom_t                          oy,
        geom_t                          oz,
        geom_t                          e0x,
        geom_t                          e0y,
        geom_t                          e0z,
        geom_t                          e1x,
        geom_t                          e1y,
        geom_t                          e1z,
        int                             degree,
        const geom_t                   *coeffs,
        ptrdiff_t                       n_coeffs) {
    return std::make_shared<PolynomialSurfaceParametrization>(
            nodeset, ox, oy, oz, e0x, e0y, e0z, e1x, e1y, e1z, degree, coeffs, n_coeffs);
}

int PolynomialSurfaceParametrization::apply(Mesh &mesh) const {
    struct Ctx {
        geom_t        ox, oy, oz;
        geom_t        e0x, e0y, e0z;
        geom_t        e1x, e1y, e1z;
        geom_t        e2x, e2y, e2z;
        int           degree;
        const geom_t *coeffs;
    } ctx{ox_,  oy_,  oz_,  e0x_, e0y_, e0z_, e1x_, e1y_, e1z_, e2x_, e2y_, e2z_, degree_,
          coeffs_ ? coeffs_->data() : nullptr};
    return apply_gathered(
            mesh,
            nodeset_,
            "PolynomialSurfaceParametrization::apply",
            [](ptrdiff_t n, geom_t *x, geom_t *y, geom_t *z, const void *p) {
                const auto *c = static_cast<const Ctx *>(p);
                project_polynomial_surface(n,
                                           x,
                                           y,
                                           z,
                                           c->ox,
                                           c->oy,
                                           c->oz,
                                           c->e0x,
                                           c->e0y,
                                           c->e0z,
                                           c->e1x,
                                           c->e1y,
                                           c->e1z,
                                           c->e2x,
                                           c->e2y,
                                           c->e2z,
                                           c->degree,
                                           c->coeffs);
            },
            &ctx);
}

IdentityParametrization::IdentityParametrization(const std::shared_ptr<Nodeset> &nodeset) : nodeset_(nodeset) {
    if (require_nodeset(nodeset_, "IdentityParametrization") != SMESH_SUCCESS) {
        return;
    }
}

std::shared_ptr<IdentityParametrization>
IdentityParametrization::create(const std::shared_ptr<Nodeset> &nodeset) {
    return std::make_shared<IdentityParametrization>(nodeset);
}

int IdentityParametrization::apply(Mesh &mesh) const {
    SMESH_UNUSED(mesh);
    if (require_nodeset(nodeset_, "IdentityParametrization::apply") != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }
    return SMESH_SUCCESS;
}

}  // namespace smesh
