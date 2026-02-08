#ifndef SMESH_SFC_HPP
#define SMESH_SFC_HPP

#include "smesh_types.hpp"

namespace smesh {
template <typename geom_t>
int encode_morton3(const ptrdiff_t n_points,
                   const geom_t *const SMESH_RESTRICT x,
                   const geom_t *const SMESH_RESTRICT y,
                   const geom_t *const SMESH_RESTRICT z,
                   u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_hilbert3(const ptrdiff_t n_points,
                    const geom_t *const SMESH_RESTRICT x,
                    const geom_t *const SMESH_RESTRICT y,
                    const geom_t *const SMESH_RESTRICT z,
                    u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_cartesian3(const ptrdiff_t n_points,
                      const geom_t *const SMESH_RESTRICT x,
                      const geom_t *const SMESH_RESTRICT y,
                      const geom_t *const SMESH_RESTRICT z, int fast, int mid,
                      int slow, u32 *const SMESH_RESTRICT encoding);

} // namespace smesh

#endif // SMESH_SFC_HPP