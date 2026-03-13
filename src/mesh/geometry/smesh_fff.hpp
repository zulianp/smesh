#ifndef SMESH_FFF_HPP
#define SMESH_FFF_HPP

#include "smesh_types.hpp"

namespace smesh {

// TODO: once we need 2D
// template <typename AdjugateType, typename DeterminantType>
// int tri3_fff_fill(
//     const ptrdiff_t nelements,
//     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
//     const ptrdiff_t stride,
//     AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

// template <typename AdjugateType, typename DeterminantType>
// int quad4_fff_fill(
//     const ptrdiff_t nelements,
//     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
//     const ptrdiff_t stride,
//     AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT fff);

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

} // namespace smesh

#endif // SMESH_FFF_HPP