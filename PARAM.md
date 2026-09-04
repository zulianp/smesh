# PARAM.md

## Goal

Snap selected mesh node coordinates onto an analytic curve, surface, or volume map. Topology sets stay topology-only.

## Contract

- `Parametrization` is a thin interface: `int apply(Mesh &mesh) const`.
- Concrete maps inherit it (`CircleParametrization`, `SphereParametrization`, `PolynomialSurfaceParametrization`, `IdentityParametrization`).
- Each subclass holds a `Nodeset` — that is how it knows which coordinates it may change.
- The caller builds the nodeset (`create_nodeset_from_sideset`, `create_nodeset_from_edgeset`, unique nodes of a block, or a hand-picked id list) and passes it to `create(...)`.
- Do not store a Sideset, Edgeset, or Block on the parametrization.
- Do not add fields to Sideset / Edgeset / Block / Nodeset.
- `apply` gathers those node coords from `mesh.points()`, snaps them onto the map, and scatters back. Owned and ghost points are already local under MPI.
- Mesh owns an optional named registry (`add_parametrization`, `parametrizations`). Clone copies the `shared_ptr` entries. No folder IO in this phase.
- `refine()` / `promote_to` are unchanged. After subdivision, construct a new parametrization with a remapped or rebuilt nodeset and `apply`. `map_nodeset_through_refine` does not insert new mid-edge nodes.

## Maps

| Type | Snap |
|------|------|
| `CircleParametrization` | Closest point on a plane circle (center, axis, radius). A node on the axis is sent along a default in-plane direction. |
| `SphereParametrization` | Radial from center. A node at the center is sent to `center + (radius, 0, 0)`. |
| `PolynomialSurfaceParametrization` | Graph `ζ = p(ξ, η)` in a local orthonormal frame. `p` is a total-degree polynomial; coefficients are `ξ^i η^j` with `i` outer, `j` inner, `i + j <= degree`. Projection keeps `(ξ, η)` and replaces the normal coordinate. |
| `IdentityParametrization` | No-op (placeholder volume/slot). |

Eval (`u,v,w` → xyz) is deferred.

## Files

- `src/mesh/geometry/smesh_parametrization.hpp`
- `src/mesh/geometry/smesh_parametrization.cpp`
- `src/mesh/geometry/tests/smesh_parametrization_test.cpp`

## Non-goals (this phase)

- `apply_eval` / public `eval` / UV fields
- Public `bind` / `nodeset()` / `dim()` / `project`
- Kind enum / packed-coeff switch
- Changing `refine()` / `promote_to` / SS fill-points
- NURBS / CAD
- Volume maps beyond `IdentityParametrization`
- Mesh folder IO for parametrizations
