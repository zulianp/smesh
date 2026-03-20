#include "smesh_device_arrays.hpp"
#include "smesh_cuda_base.cuh"
#include "smesh_types.hpp"

namespace smesh {

namespace device {
void memset(void *ptr, int value, const size_t n) {
  SMESH_CUDA_CHECK(cudaMemset(ptr, value, n));
}

bool is_ptr_device(const void *ptr) {
  cudaPointerAttributes attributes;
  cudaError_t err = cudaPointerGetAttributes(&attributes, ptr);

  if (err != cudaSuccess) {
    SMESH_ERROR("cudaPointerGetAttributes failed: %s\n",
                cudaGetErrorString(err));
  }

#if CUDART_VERSION >= 10000
  // CUDA 10.0 and newer
  return attributes.type == cudaMemoryTypeDevice;
#else
  return attributes.memoryType == cudaMemoryTypeDevice;
#endif
}

template <typename T> T *alloc(const size_t n) {
  T *ptr = nullptr;
  SMESH_CUDA_CHECK(cudaMalloc((void **)&ptr, n * sizeof(T)));
  SMESH_CUDA_CHECK(cudaMemset(ptr, 0, n * sizeof(T)));
  SMESH_ASSERT(ptr);
  return ptr;
}

void destroy(void *a) { SMESH_CUDA_CHECK(cudaFree(a)); }

template <typename T>
void copy(const size_t n, const T *const src, T *const dest) {
  SMESH_CUDA_CHECK(
      cudaMemcpy(dest, src, n * sizeof(T), cudaMemcpyDeviceToDevice));
}

template <typename T>
void device_to_host(const size_t n, const T *const d, T *h) {
  SMESH_CUDA_CHECK(cudaMemcpy(h, d, n * sizeof(T), cudaMemcpyDeviceToHost));
}

template <typename T>
void host_to_device(const size_t n, const T *const h, T *d) {
  SMESH_CUDA_CHECK(cudaMemcpy(d, h, n * sizeof(T), cudaMemcpyHostToDevice));
}

} // namespace device

} // namespace smesh

#define SMESH_INSTANTIATE_DEVICE_ARRAYS(T)                                     \
  template T *smesh::device::alloc<T>(const size_t n);                         \
  template void smesh::device::copy<T>(const size_t n, const T *const src,     \
                                       T *const dest);                         \
  template void smesh::device::device_to_host<T>(const size_t n,               \
                                                 const T *const d, T *h);      \
  template void smesh::device::host_to_device<T>(const size_t n,               \
                                                 const T *const h, T *d);      \
  template T **smesh::device::alloc<T *>(const size_t n);                      \
  template void smesh::device::copy<T *>(const size_t n, T *const *const src,  \
                                         T **const dest);                      \
  template void smesh::device::device_to_host<T *>(const size_t n,             \
                                                   T *const *const d, T **h);  \
  template void smesh::device::host_to_device<T *>(const size_t n,             \
                                                   T *const *const h, T **d);

namespace smesh {
SMESH_INSTANTIATE_DEVICE_ARRAYS(f16)
SMESH_INSTANTIATE_DEVICE_ARRAYS(char)
#if defined(__clang__)
SMESH_INSTANTIATE_DEVICE_ARRAYS(long)
#endif
SMESH_INSTANTIATE_DEVICE_ARRAYS(i8)
SMESH_INSTANTIATE_DEVICE_ARRAYS(i16)
SMESH_INSTANTIATE_DEVICE_ARRAYS(i32)
SMESH_INSTANTIATE_DEVICE_ARRAYS(i64)
SMESH_INSTANTIATE_DEVICE_ARRAYS(u8)
SMESH_INSTANTIATE_DEVICE_ARRAYS(u16)
SMESH_INSTANTIATE_DEVICE_ARRAYS(u32)
SMESH_INSTANTIATE_DEVICE_ARRAYS(u64)
SMESH_INSTANTIATE_DEVICE_ARRAYS(f32)
SMESH_INSTANTIATE_DEVICE_ARRAYS(f64)
} // namespace smesh

#undef SMESH_INSTANTIATE_DEVICE_ARRAYS
