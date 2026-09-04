#ifndef SMESH_PARAMETRIZATION_HPP
#define SMESH_PARAMETRIZATION_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"

#include <memory>

namespace smesh {

/// Analytic map that snaps selected mesh node coordinates onto a curve,
/// surface, or volume. Topology sets stay topology-only: pass a Nodeset
/// (built from a sideset, edgeset, block, or a hand-picked id list) into
/// the subclass constructor. Public API is `apply(mesh)` only.
class Parametrization {
public:
    virtual ~Parametrization() = default;
    virtual int apply(Mesh &mesh) const = 0;
};

/// Closest-point projection onto a plane circle (center, unit axis, radius).
class CircleParametrization final : public Parametrization {
public:
    CircleParametrization(const std::shared_ptr<Nodeset> &nodeset,
                          geom_t                          cx,
                          geom_t                          cy,
                          geom_t                          cz,
                          geom_t                          ax,
                          geom_t                          ay,
                          geom_t                          az,
                          geom_t                          radius);

    static std::shared_ptr<CircleParametrization> create(const std::shared_ptr<Nodeset> &nodeset,
                                                         geom_t                          cx,
                                                         geom_t                          cy,
                                                         geom_t                          cz,
                                                         geom_t                          ax,
                                                         geom_t                          ay,
                                                         geom_t                          az,
                                                         geom_t                          radius);

    int apply(Mesh &mesh) const override;

private:
    std::shared_ptr<Nodeset> nodeset_;
    geom_t                   cx_, cy_, cz_;
    geom_t                   ax_, ay_, az_;
    geom_t                   radius_;
};

/// Radial projection onto a sphere. A node at the center is sent to
/// `center + (radius, 0, 0)`.
class SphereParametrization final : public Parametrization {
public:
    SphereParametrization(const std::shared_ptr<Nodeset> &nodeset,
                          geom_t                          cx,
                          geom_t                          cy,
                          geom_t                          cz,
                          geom_t                          radius);

    static std::shared_ptr<SphereParametrization> create(const std::shared_ptr<Nodeset> &nodeset,
                                                         geom_t                          cx,
                                                         geom_t                          cy,
                                                         geom_t                          cz,
                                                         geom_t                          radius);

    int apply(Mesh &mesh) const override;

private:
    std::shared_ptr<Nodeset> nodeset_;
    geom_t                   cx_, cy_, cz_;
    geom_t                   radius_;
};

/// Graph `ζ = p(ξ, η)` in a local orthonormal frame. `p` is a total-degree
/// polynomial; coefficients are `ξ^i η^j` with `i` outer, `j` inner,
/// `i + j <= degree` (`n = (degree+1)*(degree+2)/2`). Projection keeps
/// `(ξ, η)` and replaces the normal coordinate with `p(ξ, η)`.
class PolynomialSurfaceParametrization final : public Parametrization {
public:
    PolynomialSurfaceParametrization(const std::shared_ptr<Nodeset> &nodeset,
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
                                     ptrdiff_t                       n_coeffs);

    static std::shared_ptr<PolynomialSurfaceParametrization>
    create(const std::shared_ptr<Nodeset> &nodeset,
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
           ptrdiff_t                       n_coeffs);

    int apply(Mesh &mesh) const override;

private:
    std::shared_ptr<Nodeset> nodeset_;
    geom_t                   ox_, oy_, oz_;
    geom_t                   e0x_, e0y_, e0z_;
    geom_t                   e1x_, e1y_, e1z_;
    geom_t                   e2x_, e2y_, e2z_;
    int                      degree_;
    SharedBuffer<geom_t>     coeffs_;
};

/// Placeholder volume/slot map. `apply` does not move coordinates.
class IdentityParametrization final : public Parametrization {
public:
    explicit IdentityParametrization(const std::shared_ptr<Nodeset> &nodeset);

    static std::shared_ptr<IdentityParametrization> create(const std::shared_ptr<Nodeset> &nodeset);

    int apply(Mesh &mesh) const override;

private:
    std::shared_ptr<Nodeset> nodeset_;
};

}  // namespace smesh

#endif  // SMESH_PARAMETRIZATION_HPP
