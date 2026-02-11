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

// TODO: explicit instantiation of functions in smesh_buffer.impl.hpp
} // namespace smesh
