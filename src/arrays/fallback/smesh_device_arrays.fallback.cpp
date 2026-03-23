#include "smesh_device_arrays.hpp"
#include "smesh_types.hpp"

#include <stddef.h>
#include <stdint.h>
#include <cstring>

namespace smesh {

    namespace device {
        void memset(void *ptr, int value, const size_t n) { ::memset(ptr, value, n); }

        bool is_ptr_device(const void * /*ptr*/) { return false; }

        template <typename T>
        T *alloc(const size_t n) {
            auto ptr = (T *)malloc(n * sizeof(T));
            SMESH_ASSERT(ptr);
            memset(ptr, 0, n * sizeof(T));

            return ptr;
        }

        template <typename T>
        void destroy(void *a) {
            free(a);
        }

        template <typename T>
        void copy(const size_t n, const T *const src, T *const dest) {
            memcpy(dest, src, n * sizeof(T));
        }

        template <typename T>
        void device_to_host(const size_t n, const T *const d, T *h) {
            memcpy(h, d, n * sizeof(T));
        }

        template <typename T>
        void host_to_device(const size_t n, const T *const h, T *d) {
            memcpy(d, h, n * sizeof(T));
        }

#define SMESH_INSTANTIATE_DEVICE_ARRAYS(T)                                                       \
    template T   *smesh::device::alloc<T>(const size_t n);                                       \
    template void smesh::device::copy<T>(const size_t n, const T *const src, T *const dest);     \
    template void smesh::device::device_to_host<T>(const size_t n, const T *const d, T *h);      \
    template void smesh::device::host_to_device<T>(const size_t n, const T *const h, T *d);      \
    template T  **smesh::device::alloc<T *>(const size_t n);                                     \
    template void smesh::device::copy<T *>(const size_t n, T *const *const src, T **const dest); \
    template void smesh::device::device_to_host<T *>(const size_t n, T *const *const d, T **h);  \
    template void smesh::device::host_to_device<T *>(const size_t n, T *const *const h, T **d);

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

    }  // namespace device

}  // namespace smesh

#undef SMESH_INSTANTIATE_DEVICE_ARRAYS
