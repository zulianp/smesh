#ifndef SMESH_DEVICE_BUFFER_HPP
#define SMESH_DEVICE_BUFFER_HPP

#include "smesh_buffer.hpp"

#include <vector>

namespace smesh {

template <typename T> SharedBuffer<T> create_device_buffer(const ptrdiff_t n);

template <typename T>
SharedBuffer<T *> create_device_buffer(const ptrdiff_t n0, const ptrdiff_t n1);

template <typename T> SharedBuffer<T> to_device(const SharedBuffer<T> &in);

template <typename T>
std::vector<SharedBuffer<T>> to_device(const std::vector<SharedBuffer<T>> &in);

template <typename T> SharedBuffer<T *> to_device(const SharedBuffer<T *> &in);

template <typename T> SharedBuffer<T> to_host(const SharedBuffer<T> &in);

template <typename T> SharedBuffer<T *> to_host(const SharedBuffer<T *> &in);

} // namespace smesh

#endif // SMESH_DEVICE_BUFFER_HPP
