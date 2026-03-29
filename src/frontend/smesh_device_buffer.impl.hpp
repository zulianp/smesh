#ifndef SMESH_DEVICE_BUFFER_IMPL_HPP
#define SMESH_DEVICE_BUFFER_IMPL_HPP

#include "smesh_device_arrays.hpp"
#include "smesh_alloc.hpp"

#include "smesh_device_buffer.hpp"

#include <cstdlib>
#include <memory>
#include <type_traits>

namespace smesh {

    template <typename T>
    SharedBuffer<T> create_device_buffer(const ptrdiff_t n) {
#ifndef SMESH_ENABLE_CUDA
        return create_host_buffer<T>(n);
#else
        auto ret = std::make_shared<Buffer<T>>(n, (T *)device::alloc<T>(n), &device::destroy, MEMORY_SPACE_DEVICE);
        return ret;
#endif  // SMESH_ENABLE_CUDA
    }

    template <typename T>
    SharedBuffer<T *> create_device_buffer(const ptrdiff_t n0, const ptrdiff_t n1) {
#ifndef SMESH_ENABLE_CUDA
        return create_host_buffer<T>(n0, n1);
#else
        T **dev_buff0 = (T **)device::alloc<T *>(n0);

        std::vector<T *> host_dev_ptrs(n0);
        for (ptrdiff_t i0 = 0; i0 < n0; i0++) {
            host_dev_ptrs[i0] = (T *)device::alloc<T>(n1);
            device::memset(host_dev_ptrs[i0], 0, n1 * sizeof(T));
        }

        device::host_to_device(n0, host_dev_ptrs.data(), dev_buff0);

        return std::make_shared<Buffer<T *>>(
                n0,
                n1,
                dev_buff0,
                [n0, host_dev_ptrs](int /*n*/, void **ptr) {
                    for (ptrdiff_t i0 = 0; i0 < n0; i0++) {
                        device::destroy(host_dev_ptrs[i0]);
                    }

                    device::destroy(ptr);
                },
                MEMORY_SPACE_DEVICE);
#endif  // SMESH_ENABLE_CUDA
    }

    template <typename T>
    SharedBuffer<T> to_device(const SharedBuffer<T> &in) {
        if (!in) {
            SMESH_ERROR("Input buffer is null");
        }

#ifndef SMESH_ENABLE_CUDA
        return in;
#else
        if (in->mem_space() == MEMORY_SPACE_DEVICE) {
            return in;
        }

        T *buff = (T *)device::alloc<T>(in->size());
        device::host_to_device(in->size(), in->data(), buff);

        return std::make_shared<Buffer<T>>(in->size(), buff, &device::destroy, MEMORY_SPACE_DEVICE);
#endif  // SMESH_ENABLE_CUDA
    }

    template <typename T>
    std::vector<SharedBuffer<T>> to_device(const std::vector<SharedBuffer<T>> &in) {
#ifndef SMESH_ENABLE_CUDA
        return in;
#else
        std::vector<SharedBuffer<T>> ret;
        for (const auto &b : in) {
            ret.push_back(to_device(b));
        }
        return ret;
#endif  // SMESH_ENABLE_CUDA
    }

    template <typename T>
    SharedBuffer<T *> to_device(const SharedBuffer<T *> &in) {
        if (!in) {
            SMESH_ERROR("Input buffer is null");
        }

#ifndef SMESH_ENABLE_CUDA
        return in;
#else
        if (in->mem_space() == MEMORY_SPACE_DEVICE) {
            return in;
        }

        size_t n0 = in->extent(0);
        size_t n1 = in->extent(1);

        T **dev_buff0 = (T **)device::alloc<T *>(n0);

        std::vector<T *> host_dev_ptrs(n0);
        for (size_t i0 = 0; i0 < n0; i0++) {
            host_dev_ptrs[i0] = (T *)device::alloc<T>(n1);
            device::host_to_device(n1, in->data()[i0], host_dev_ptrs[i0]);
        }

        device::host_to_device(n0, host_dev_ptrs.data(), dev_buff0);

        return std::make_shared<Buffer<T *>>(
                n0,
                n1,
                dev_buff0,
                [n0, host_dev_ptrs](int /*n*/, void **ptr) {
                    for (size_t i0 = 0; i0 < n0; i0++) {
                        device::destroy(host_dev_ptrs[i0]);
                    }

                    device::destroy(ptr);
                },
                MEMORY_SPACE_DEVICE);
#endif  // SMESH_ENABLE_CUDA
    }

    template <typename T>
    SharedBuffer<T> to_host(const SharedBuffer<T> &in) {
        if (!in) {
            SMESH_ERROR("Input buffer is null");
        }

#ifndef SMESH_ENABLE_CUDA
        return in;
#else
        if (in->mem_space() == MEMORY_SPACE_HOST) {
            return in;
        }

        using NonConstT = typename std::remove_const<T>::type;

        NonConstT *buff = static_cast<NonConstT *>(SMESH_ALLOC(in->size() * sizeof(NonConstT)));
        SMESH_ASSERT(buff != nullptr);
        device::device_to_host(in->size(), in->data(), buff);
        return std::make_shared<Buffer<T>>(in->size(), buff, &free, MEMORY_SPACE_HOST);
#endif  // SMESH_ENABLE_CUDA
    }

    template <typename T>
    std::shared_ptr<Buffer<T *>> to_host(const std::shared_ptr<Buffer<T *>> &in) {
        if (!in) {
            SMESH_ERROR("Input buffer is null");
        }

#ifndef SMESH_ENABLE_CUDA
        return in;
#else
        if (in->mem_space() == MEMORY_SPACE_HOST) {
            return in;
        }

        const ptrdiff_t n0 = in->extent(0);
        const ptrdiff_t n1 = in->extent(1);

        auto buffer = create_host_buffer<T>(n0, n1);

        T **dev_addr = static_cast<T **>(SMESH_ALLOC(n0 * sizeof(T *)));
        SMESH_ASSERT(dev_addr != nullptr);

        device::device_to_host(n0, in->data(), dev_addr);

        for (ptrdiff_t i = 0; i < n0; i++) {
            device::device_to_host(n1, dev_addr[i], buffer->data()[i]);
        }

        SMESH_FREE(dev_addr);

        return buffer;
#endif  // SMESH_ENABLE_CUDA
    }

}  // namespace smesh

#endif  // SMESH_DEVICE_BUFFER_IMPL_HPP
