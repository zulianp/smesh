#include "smesh_point_transform.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_POINT_TRANSFORM(T)                          \
  template void shift<T>(const ptrdiff_t n_nodes, const T dx,                  \
                         T *const SMESH_RESTRICT x);                           \
  template void shift3<T>(const ptrdiff_t n_nodes, const T dx, const T dy,     \
                          const T dz, T *const SMESH_RESTRICT x,               \
                          T *const SMESH_RESTRICT y,                           \
                          T *const SMESH_RESTRICT z);                          \
  template void scale<T>(const ptrdiff_t n_nodes, const T s,                   \
                         T *const SMESH_RESTRICT x);                           \
  template void scale3<T>(const ptrdiff_t n_nodes, const T sx, const T sy,     \
                          const T sz, T *const SMESH_RESTRICT x,               \
                          T *const SMESH_RESTRICT y,                           \
                          T *const SMESH_RESTRICT z);                          \
  template void rotate2<T>(const ptrdiff_t n_nodes, const T angle,             \
                           T *const SMESH_RESTRICT x,                          \
                           T *const SMESH_RESTRICT y);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_POINT_TRANSFORM(f32);
SMESH_EXPLICIT_INSTANTIATE_POINT_TRANSFORM(f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_POINT_TRANSFORM