#ifndef SMESH_FFF_HPP
#define SMESH_FFF_HPP

#include <cstddef>

#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

namespace smesh {

/// Packed FFF SoA length: 3 for planar TRI3/QUAD4 (symmetric UT), 6 for 3D families.
inline int fff_components(const enum ElemType element_type) {
    switch (element_type) {
        case TRI3:
        case QUAD4:
            return 3;
        default:
            return 6;
    }
}

template <typename FFFType>
int tri3_fff_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

template <typename FFFType>
int quad4_fff_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const geom_t qx, const geom_t qy, const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

template <typename FFFType>
int tet4_fff_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

template <typename FFFType>
int hex8_fff_fill(
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const geom_t qx, const geom_t qy, const geom_t qz, const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

/// Per-macro-element FFF from the eight SSHEX8 macro corners (isoparametric ref. point 1/2,1/2,1/2).
template <typename FFFType>
int sshex8_macro_fff_fill(
    const int level,
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    const ptrdiff_t stride,
    FFFType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

} // namespace smesh

#endif // SMESH_FFF_HPP