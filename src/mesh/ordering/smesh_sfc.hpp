#ifndef SMESH_SFC_HPP
#define SMESH_SFC_HPP

#include "smesh_types.hpp"

namespace smesh {

/// Grows the box [box_min, box_max] into a cube, every side the longest one, from the same lower
/// corner. The encoders quantise each axis over the range they are given, so a box of different
/// sides gives a thin axis as many cells as a long one: a curve over it cuts across the thin
/// direction as often as along the long ones. Partitions taken as ranges of such a curve are
/// slices through a plate's thickness, with interfaces as large as the plate (FF200's valves, a
/// box of 31.8 x 1.15 x 19.6 mm: five times the halo of a cubic box at 32 ranks).
template <typename geom_t>
inline void sfc_cubic_box(const int spatial_dim, const geom_t *const box_min,
                          geom_t *const box_max) {
  geom_t side = 0;
  for (int d = 0; d < spatial_dim; ++d) {
    const geom_t range = box_max[d] - box_min[d];
    if (range > side) side = range;
  }
  for (int d = 0; d < spatial_dim; ++d) {
    box_max[d] = box_min[d] + side;
  }
}

template <typename geom_t>
int encode_morton3(const ptrdiff_t n_points,
                   const geom_t *const SMESH_RESTRICT x,
                   const geom_t *const SMESH_RESTRICT y,
                   const geom_t *const SMESH_RESTRICT z,
                   u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_morton3(const ptrdiff_t n_points,
                   const geom_t *const SMESH_RESTRICT x,
                   const geom_t *const SMESH_RESTRICT y,
                   const geom_t *const SMESH_RESTRICT z,
                   const geom_t x_min, const geom_t x_max,
                   const geom_t y_min, const geom_t y_max,
                   const geom_t z_min, const geom_t z_max,
                   u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_hilbert3(const ptrdiff_t n_points,
                    const geom_t *const SMESH_RESTRICT x,
                    const geom_t *const SMESH_RESTRICT y,
                    const geom_t *const SMESH_RESTRICT z,
                    u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_hilbert3(const ptrdiff_t n_points,
                    const geom_t *const SMESH_RESTRICT x,
                    const geom_t *const SMESH_RESTRICT y,
                    const geom_t *const SMESH_RESTRICT z,
                    const geom_t x_min, const geom_t x_max,
                    const geom_t y_min, const geom_t y_max,
                    const geom_t z_min, const geom_t z_max,
                    u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_cartesian3(const ptrdiff_t n_points,
                      const geom_t *const SMESH_RESTRICT x,
                      const geom_t *const SMESH_RESTRICT y,
                      const geom_t *const SMESH_RESTRICT z, int fast, int mid,
                      int slow, u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_cartesian3(const ptrdiff_t n_points,
                      const geom_t *const SMESH_RESTRICT x,
                      const geom_t *const SMESH_RESTRICT y,
                      const geom_t *const SMESH_RESTRICT z,
                      const geom_t x_min, const geom_t x_max,
                      const geom_t y_min, const geom_t y_max,
                      const geom_t z_min, const geom_t z_max, int fast, int mid,
                      int slow, u32 *const SMESH_RESTRICT encoding);

template <typename geom_t>
int encode_random3(const ptrdiff_t n_points,
                  const geom_t *const SMESH_RESTRICT x,
                  const geom_t *const SMESH_RESTRICT y,
                  const geom_t *const SMESH_RESTRICT z,
                  u32 *const SMESH_RESTRICT encoding);

} // namespace smesh

#endif // SMESH_SFC_HPP
