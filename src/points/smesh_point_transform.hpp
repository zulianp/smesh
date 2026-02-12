#ifndef SMESH_POINT_TRANSFORM_HPP
#define SMESH_POINT_TRANSFORM_HPP

#include "smesh_base.hpp"

namespace smesh {

    template <typename T>
    void shift(const ptrdiff_t n_nodes, const T dx, T *const SMESH_RESTRICT x);
    template <typename T>
    void shift3(const ptrdiff_t n_nodes, const T dx, const T dy, const T dz,
                T *const SMESH_RESTRICT x, T *const SMESH_RESTRICT y,
                T *const SMESH_RESTRICT z);
    template <typename T>
    void scale(const ptrdiff_t n_nodes, const T s, T *const SMESH_RESTRICT x);
    template <typename T>
    void scale3(const ptrdiff_t n_nodes, const T sx, const T sy, const T sz,
                T *const SMESH_RESTRICT x, T *const SMESH_RESTRICT y,
                T *const SMESH_RESTRICT z);
    template <typename T>
    void rotate2(const ptrdiff_t n_nodes, const T angle, T *const SMESH_RESTRICT x,
                T *const SMESH_RESTRICT y);


} // namespace smesh

#endif // SMESH_POINT_TRANSFORM_HPP