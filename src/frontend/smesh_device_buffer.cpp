#include "smesh_device_buffer.impl.hpp"

namespace smesh {

template <> SharedBuffer<void> to_device<void>(const SharedBuffer<void> &in) {
#ifndef SMESH_ENABLE_CUDA
  return in;
#else
  if (in->mem_space() == MEMORY_SPACE_DEVICE) {
    return in;
  }

  const size_t nbytes = in->nbytes();
  char *buff = device::alloc<char>(nbytes);
  device::host_to_device<char>(
      nbytes, static_cast<const char *>(in->void_data()), buff);

  auto out = Buffer<void>::make_empty();
  out->n_ = in->size();
  out->type_ = in->element_type();
  out->ptr_ = buff;
  out->destroy_ = &device::destroy;
  out->mem_space_ = MEMORY_SPACE_DEVICE;

  return out;
#endif
}

template <>
SharedBuffer<void *> to_device<void>(const SharedBuffer<void *> &in) {
#ifndef SMESH_ENABLE_CUDA
  return in;
#else
  if (in->mem_space() == MEMORY_SPACE_DEVICE) {
    return in;
  }

  const size_t n0 = in->extent(0);
  const size_t n1 = in->extent(1);

  auto dev_buff0 =
      reinterpret_cast<void **>(device::alloc<char>(n0 * sizeof(void *)));

  std::vector<void *> host_dev_ptrs(n0);
  for (size_t i0 = 0; i0 < n0; i0++) {
    host_dev_ptrs[i0] = device::alloc<char>(n1);
    device::host_to_device<char>(n1, static_cast<const char *>(in->data()[i0]),
                                 static_cast<char *>(host_dev_ptrs[i0]));
  }

  device::host_to_device<char>(
      n0 * sizeof(void *), reinterpret_cast<const char *>(host_dev_ptrs.data()),
      reinterpret_cast<char *>(dev_buff0));

  return std::make_shared<Buffer<void *>>(
      n0, n1, dev_buff0,
      [n0, host_dev_ptrs](int /*n*/, void **ptr) {
        for (size_t i0 = 0; i0 < n0; i0++) {
          device::destroy(host_dev_ptrs[i0]);
        }

        device::destroy(ptr);
      },
      MEMORY_SPACE_DEVICE);
#endif
}

#define SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(T)                         \
  template SharedBuffer<T> create_device_buffer<T>(const ptrdiff_t);           \
  template SharedBuffer<T> to_device<T>(const SharedBuffer<T> &);              \
  template std::vector<SharedBuffer<T>> to_device<T>(                          \
      const std::vector<SharedBuffer<T>> &);                                   \
  template SharedBuffer<T> to_host<T>(const SharedBuffer<T> &)

SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(f16);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(f32);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(f64);
template SharedBuffer<const f32>
to_host<const f32>(const SharedBuffer<const f32> &);
template SharedBuffer<const f64>
to_host<const f64>(const SharedBuffer<const f64> &);

SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(i8);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(i16);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(i32);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(i64);

SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(u8);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(u16);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(u32);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(u64);

SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(char);

#if defined(__clang__)
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D(long);
#endif

#undef SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_1D

#define SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(T)                         \
  template SharedBuffer<T *> create_device_buffer<T>(const ptrdiff_t,          \
                                                     const ptrdiff_t);         \
  template SharedBuffer<T *> to_device<T>(const SharedBuffer<T *> &);          \
  template SharedBuffer<T *> to_host<T>(const SharedBuffer<T *> &);            \
  template SharedBuffer<T *> create_2d_device<T>(                              \
      const std::vector<SharedBuffer<T>> &)

SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(f16);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(f32);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(f64);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(i32);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(i64);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(i16);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(u8);
SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D(u16);

#undef SMESH_EXPLICIT_INSTANTIATE_DEVICE_BUFFER_2D

} // namespace smesh
