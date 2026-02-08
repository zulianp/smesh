#ifndef SMESH_COMMON_HPP
#define SMESH_COMMON_HPP

#include "smesh_base.hpp"

#include <algorithm>
#include <cmath>

namespace smesh {
template <typename T>
static SMESH_INLINE void normalize3(T *const SMESH_RESTRICT ax,
                                    T *const SMESH_RESTRICT ay,
                                    T *const SMESH_RESTRICT az) noexcept {
  const T inv_len = 1 / sqrt(*ax * *ax + *ay * *ay + *az * *az);
  *ax *= inv_len;
  *ay *= inv_len;
  *az *= inv_len;
}

template <typename T>
static SMESH_INLINE T dot3(const T ax, const T ay, const T az, const T bx,
                           const T by, const T bz) noexcept {
  return ax * bx + ay * by + az * bz;
}

template <typename T>
static SMESH_INLINE void
normal3(const T ax, const T ay, const T az, const T bx, const T by, const T bz,
        const T cx, const T cy, const T cz, T *const SMESH_RESTRICT nx,
        T *const SMESH_RESTRICT ny, T *const SMESH_RESTRICT nz) {
  T ux = bx - ax, uy = by - ay, uz = bz - az;
  T vx = cx - ax, vy = cy - ay, vz = cz - az;

  normalize3(&ux, &uy, &uz);
  normalize3(&vx, &vy, &vz);

  *nx = uy * vz - uz * vy;
  *ny = uz * vx - ux * vz;
  *nz = ux * vy - uy * vx;

  normalize3(nx, ny, nz);
}

template <class T> static inline T clamp(T v, T lo, T hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}

template <typename T>
void minmax(const ptrdiff_t n, const T *const SMESH_RESTRICT values,
            T *const SMESH_RESTRICT min, T *const SMESH_RESTRICT max) {
  if (n <= 0) {
    return;
  }

  T min_v = values[0];
  T max_v = values[0];
#pragma omp parallel for reduction(min : min_v) reduction(max : max_v)
  for (ptrdiff_t i = 1; i < n; i++) {
    min_v = std::min(min_v, values[i]);
    max_v = std::max(max_v, values[i]);
  }

  *min = min_v;
  *max = max_v;
}
} // namespace smesh

#endif // SMESH_COMMON_HPP