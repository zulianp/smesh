#ifndef SMESH_DEVICE_BUFFER_HPP
#define SMESH_DEVICE_BUFFER_HPP

#include "smesh_buffer.hpp"

#include <vector>

#ifdef SMESH_ENABLE_CUDA
#include "smesh_device_arrays.hpp"
#endif

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

template <typename T>
SharedBuffer<T *> create_2d_device(const std::vector<SharedBuffer<T>> &buffers);

template <typename T>
SharedBuffer<T *> create_2d(const std::vector<SharedBuffer<T>> &buffers) {
  if (buffers.empty()) {
    SMESH_ERROR("No buffers to create 2D buffer");
    return nullptr;
  }

  int n0 = buffers.size();
  int n1 = buffers[0]->size();

#ifdef SMESH_ENABLE_CUDA
  if (buffers[0]->mem_space() == MEMORY_SPACE_DEVICE) {
    return create_2d_device(buffers);
  }
#endif

  auto mem = (T **)SMESH_ALLOC(n0 * sizeof(T *));
  for (int i0 = 0; i0 < n0; i0++) {
    mem[i0] = buffers[i0]->data();
  }

  return Buffer<T *>::own(
      n0, n1, mem,
      [keep_alive = buffers](int /*n*/, void **ptr) {
        // Internal buffers are freed by "buffers" destructor
        SMESH_FREE(ptr);
      },
      MEMORY_SPACE_HOST);
}

#ifdef SMESH_ENABLE_CUDA
template <typename T>
inline SharedBuffer<T *>
create_2d_device(const std::vector<SharedBuffer<T>> &buffers) {
  if (buffers.empty()) {
    SMESH_ERROR("No buffers to create 2D buffer");
    return nullptr;
  }

  auto n0 = buffers.size();
  auto n1 = buffers[0]->size();
  auto mem = (T **)device::alloc<T *>(n0);

  std::vector<T *> host_dev_ptrs(n0);
  for (size_t i0 = 0; i0 < n0; i0++) {
    host_dev_ptrs[i0] = buffers[i0]->data();
  }
  device::host_to_device(n0, host_dev_ptrs.data(), mem);

  return Buffer<T *>::own(
      n0, n1, mem,
      [keep_alive = buffers](int /*n*/, void **ptr) { device::destroy(ptr); },
      MEMORY_SPACE_DEVICE);
}
#endif

} // namespace smesh

#endif // SMESH_DEVICE_BUFFER_HPP
