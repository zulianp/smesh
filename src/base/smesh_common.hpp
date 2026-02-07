#ifndef SMESH_COMMON_HPP
#define SMESH_COMMON_HPP

#include "smesh_base.hpp"

namespace smesh {
template <typename T>
static SMESH_INLINE void normalize3(T *const SMESH_RESTRICT ax,
                                    T *const SMESH_RESTRICT ay,
                                    T *const SMESH_RESTRICT az) noexcept {
  const T inv_len = 1/sqrt(*ax * *ax + *ay * *ay + *az * *az);
  *ax *= inv_len;
  *ay *= inv_len;
  *az *= inv_len;
}

template <typename T>
static SMESH_INLINE T dot3(const T ax, const T ay, const T az, const T bx,
                           const T by, const T bz) noexcept {
  return ax * bx + ay * by + az * bz;
}
} // namespace smesh

#endif // SMESH_COMMON_HPP