#include "smesh_buffer.impl.hpp"

namespace smesh {

BaseBuffer::BaseBuffer(const size_t n, const PrimitiveType type,
                       void *const ptr, std::function<void(void *)> destroy,
                       MemorySpace mem_space)
    : n_(n), type_(type), ptr_(ptr), destroy_(destroy), mem_space_(mem_space) {}

BaseBuffer::~BaseBuffer() {
  if (destroy_) {
    destroy_((void *)ptr_);
  }
}

void *BaseBuffer::void_data() { return ptr_; }
const void *BaseBuffer::void_data() const { return ptr_; }

template class Buffer<f16>;
template class Buffer<f32>;
template class Buffer<f64>;

template class Buffer<i8>;
template class Buffer<i16>;
template class Buffer<i32>;
template class Buffer<i64>;

template class Buffer<u8>;
template class Buffer<u16>;
template class Buffer<u32>;
template class Buffer<u64>;

template class Buffer<mask_t>;
template class Buffer<ptrdiff_t>;

#define SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(T)                                \
  template std::shared_ptr<Buffer<T>> create_host_buffer<T>(const size_t);     \
  template std::shared_ptr<Buffer<T>> manage_host_buffer<T>(const size_t,      \
                                                            T *);              \
  template std::shared_ptr<Buffer<T>> view<T>(                                 \
      const std::shared_ptr<Buffer<T>> &, const size_t, const size_t)

SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(f16);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(f32);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(f64);

SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(i8);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(i16);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(i32);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(i64);

SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(u8);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(u16);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(u32);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(u64);

SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(mask_t);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D(ptrdiff_t);

#undef SMESH_EXPLICIT_INSTANTIATE_BUFFER_1D

#define SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D(T)                                \
  template class Buffer<T *>;                                                  \
  template std::shared_ptr<Buffer<T *>> create_host_buffer<T>(const size_t,    \
                                                              const size_t);   \
  template std::shared_ptr<Buffer<T *>> create_host_buffer_fake_SoA<T>(        \
      const size_t, const size_t);                                             \
  template std::shared_ptr<Buffer<T *>> convert_host_buffer_to_fake_SoA<T>(    \
      const size_t, const std::shared_ptr<Buffer<T>> &);                       \
  template std::shared_ptr<Buffer<T *>> manage_host_buffer<T>(                 \
      const size_t, const size_t, T **);                                       \
  template std::shared_ptr<Buffer<T *>> view<T>(                               \
      const std::shared_ptr<Buffer<T *>> &, const size_t, const size_t,        \
      const size_t, const size_t);                                             \
  template std::shared_ptr<Buffer<T>> sub<T>(                                  \
      const std::shared_ptr<Buffer<T *>> &, const size_t);                     \
  template std::shared_ptr<Buffer<T *>> copy<T>(                               \
      const std::shared_ptr<Buffer<T *>> &);                                   \
  template std::shared_ptr<Buffer<T *>> zeros_like<T>(                         \
      const std::shared_ptr<Buffer<T *>> &)

SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D(f32);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D(i32);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D(i16);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D(u8);
SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D(u16);

#undef SMESH_EXPLICIT_INSTANTIATE_BUFFER_2D
} // namespace smesh
