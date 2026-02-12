#ifndef SMESH_POINT_TRANSFORM_IMPL_HPP
#define SMESH_POINT_TRANSFORM_IMPL_HPP

#include "smesh_point_transform.hpp"

namespace smesh {

template <typename T>
void shift(const ptrdiff_t n_nodes, const T dx, T *const SMESH_RESTRICT x) {
#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    x[i] += dx;
  }
}

template <typename T>
void shift3(const ptrdiff_t n_nodes, const T dx, const T dy, const T dz,
            T *const SMESH_RESTRICT x, T *const SMESH_RESTRICT y,
            T *const SMESH_RESTRICT z) {
  shift(n_nodes, dx, x);
  shift(n_nodes, dy, y);
  shift(n_nodes, dz, z);
}

template <typename T>
void scale(const ptrdiff_t n_nodes, const T s, T *const SMESH_RESTRICT x) {
#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    x[i] *= s;
  }
}

template <typename T>
void scale3(const ptrdiff_t n_nodes, const T sx, const T sy, const T sz,
            T *const SMESH_RESTRICT x, T *const SMESH_RESTRICT y,
            T *const SMESH_RESTRICT z) {
  scale(n_nodes, sx, x);
  scale(n_nodes, sy, y);
  scale(n_nodes, sz, z);
}

template <typename T>
void rotate2(const ptrdiff_t n_nodes, const T angle, T *const SMESH_RESTRICT x,
             T *const SMESH_RESTRICT y) {
#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    const T cos_angle = cos(angle);
    const T sin_angle = sin(angle);
    const T xi = x[i] * cos_angle - y[i] * sin_angle;
    const T yi = x[i] * sin_angle + y[i] * cos_angle;

    x[i] = xi;
    y[i] = yi;
  }
}

} // namespace smesh

#endif // SMESH_POINT_TRANSFORM_IMPL_HPP
