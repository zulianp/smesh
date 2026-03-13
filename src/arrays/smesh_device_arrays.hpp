#ifndef SMESH_DEVICE_ARRAYS_HPP
#define SMESH_DEVICE_ARRAYS_HPP

#include "smesh_types.hpp"

namespace smesh {

namespace device {

void destroy(void *a);
void memset(void *ptr, int value, const size_t n);
bool is_ptr_device(const void *ptr);

template <typename T> T *alloc(const size_t n);

template <typename T>
void copy(const size_t n, const T *const src, T *const dest);

template <typename T>
void device_to_host(const size_t n, const T *const d, T *h);

template <typename T>
void host_to_device(const size_t n, const T *const h, T *d);

} // namespace device

} // namespace smesh

#endif // SMESH_DEVICE_ARRAYS_HPP
