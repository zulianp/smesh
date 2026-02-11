#ifndef SMESH_BUFFER_IMPL_HPP
#define SMESH_BUFFER_IMPL_HPP

#include <cstdlib>
#include <cstring>
#include <type_traits>

#include "smesh_buffer.hpp"

namespace smesh {

template <typename T>
Buffer<T>::Buffer(const size_t n, T *const ptr,
                  std::function<void(void *)> destroy, MemorySpace mem_space)
    : BaseBuffer(n, TypeToEnum<T>::value(), static_cast<void *>(ptr), destroy,
                 mem_space) {}

template <typename T>
Buffer<T>::Buffer()
    : BaseBuffer(0, SMESH_TYPE_UNDEFINED, nullptr, nullptr,
                 MEMORY_SPACE_INVALID) {}

template <typename T> T *Buffer<T>::data() {
  return static_cast<T *>(this->ptr_);
}

template <typename T> const T *Buffer<T>::data() const {
  return static_cast<const T *>(this->ptr_);
}

template <typename T> void Buffer<T>::print(std::ostream &os) const {
  if (this->mem_space_ == MEMORY_SPACE_DEVICE) {
    os << "On the device!\n";
    return;
  }

  os << "Buffer size " << this->n_ << "\n";
  const auto *d = data();
  for (size_t i = 0; i < this->n_; i++) {
    os << d[i] << " ";
  }
  os << "\n";
}

template <typename T>
std::shared_ptr<Buffer<T>> Buffer<T>::wrap(const size_t n, T *const x,
                                          enum MemorySpace mem_space) {
  return std::make_shared<Buffer<T>>(n, x, nullptr, mem_space);
}

template <typename T>
std::shared_ptr<Buffer<T>> Buffer<T>::own(const size_t n, T *x,
                                         std::function<void(void *)> destroy,
                                         enum MemorySpace mem_space) {
  return std::make_shared<Buffer<T>>(n, x, destroy, mem_space);
}

template <typename T> std::shared_ptr<Buffer<T>> Buffer<T>::make_empty() {
  return std::make_shared<Buffer<T>>();
}

template <typename T> int Buffer<T>::to_file(const Path &path) const {
  return array_write_convert_from_extension<T>(path, data(), this->n_);
}

template <typename T> Buffer<T *>::~Buffer() {
  if (destroy_) {
    destroy_(static_cast<int>(extent_[0]), reinterpret_cast<void **>(ptr_));
  }
}

template <typename T> T **Buffer<T *>::data() { return ptr_; }
template <typename T> const T **Buffer<T *>::data() const { return ptr_; }

template <typename T> void Buffer<T *>::print(std::ostream &os) {
  if (mem_space_ == MEMORY_SPACE_DEVICE) {
    os << "On the device!\n";
    return;
  }

  os << "Buffer size " << extent_[0] << ", " << extent_[1] << "\n";
  for (size_t i = 0; i < extent_[0]; i++) {
    for (size_t j = 0; j < extent_[1]; j++) {
      os << ptr_[i][j] << " ";
    }
    os << "\n";
  }
  os << "\n";
}

template <typename T>
std::shared_ptr<Buffer<T *>> Buffer<T *>::wrap(const size_t n0,
                                               const size_t n1, T **x,
                                               enum MemorySpace mem_space) {
  return std::make_shared<Buffer<T *>>(n0, n1, x, nullptr, mem_space);
}

template <typename T>
std::shared_ptr<Buffer<T *>> Buffer<T *>::own(
    const size_t n0, const size_t n1, T **x,
    std::function<void(int n, void **)> destroy, enum MemorySpace mem_space) {
  return std::make_shared<Buffer<T *>>(n0, n1, x, destroy, mem_space);
}

template <typename T> int Buffer<T *>::to_files(const Path &format) {
  int ret = SMESH_SUCCESS;
  char path[2048];
  for (int i = 0; i < static_cast<int>(extent_[0]); i++) {
    int nchars = snprintf(path, sizeof(path), format.c_str(), i);
    assert(nchars < static_cast<int>(sizeof(path)));

    if (nchars >= static_cast<int>(sizeof(path))) {
      SMESH_ERROR("Path is too long!\n");
    }

    if (array_write_convert_from_extension<T>(Path(path), ptr_[i], extent_[1])) {
      ret = SMESH_FAILURE;
    }
  }

  return ret;
}

template <typename T> void Buffer<T *>::release() {
  destroy_ = nullptr;
  ptr_ = nullptr;
}

template <typename T>
std::shared_ptr<Buffer<T>> create_host_buffer(const size_t n) {
  auto ret = std::make_shared<Buffer<T>>(n, static_cast<T *>(calloc(n, sizeof(T))),
                                         &free, MEMORY_SPACE_HOST);
  return ret;
}

template <typename T>
std::shared_ptr<Buffer<T *>> create_host_buffer(const size_t n0,
                                                const size_t n1) {
  T **data = static_cast<T **>(malloc(n0 * sizeof(T *)));
  for (size_t i = 0; i < n0; ++i) {
    data[i] = static_cast<T *>(calloc(n1, sizeof(T)));
  }

  auto ret = std::make_shared<Buffer<T *>>(
      n0, n1, data,
      [=](int n, void **x) {
        for (int i = 0; i < n; ++i) {
          free(x[i]);
        }
        free(x);
      },
      MEMORY_SPACE_HOST);
  return ret;
}

template <typename T>
std::shared_ptr<Buffer<T *>> create_host_buffer_fake_SoA(const size_t n0,
                                                         const size_t n1) {
  T *allocated = static_cast<T *>(calloc(n0 * n1, sizeof(T)));

  T **data = static_cast<T **>(malloc(n0 * sizeof(T *)));
  for (size_t i = 0; i < n0; ++i) {
    data[i] = &allocated[i * n1];
  }

  auto ret = std::make_shared<Buffer<T *>>(
      n0, n1, data,
      [=](int, void **x) {
        free(x[0]);
        free(x);
      },
      MEMORY_SPACE_HOST);
  return ret;
}

template <typename T>
std::shared_ptr<Buffer<T *>>
convert_host_buffer_to_fake_SoA(const size_t n0,
                                const std::shared_ptr<Buffer<T>> &in) {
  assert(n0 > 0);
  assert(in->size() % n0 == 0);
  const size_t n1 = in->size() / n0;

  T **data = static_cast<T **>(malloc(n0 * sizeof(T *)));
  for (size_t i = 0; i < n0; ++i) {
    data[i] = const_cast<T *>(&in->data()[i * n1]);
  }

  auto ret = std::make_shared<Buffer<T *>>(
      n0, n1, data, [lifetime = in](int, void **x) { free(x); },
      MEMORY_SPACE_HOST);
  return ret;
}

template <typename T>
static std::shared_ptr<Buffer<T>>
soa_to_aos(const size_t in_stride0, const size_t in_stride1,
           const std::shared_ptr<Buffer<T *>> &in) {
  auto n0 = in->extent(0);
  auto n1 = in->extent(1);

  auto out = smesh::create_host_buffer<T>(n0 * n1);

  auto d_in = in->data();
  auto d_out = out->data();

  for (size_t i = 0; i < n0; i++) {
    for (size_t j = 0; j < n1; j++) {
      d_out[i * in_stride0 + j * in_stride1] = d_in[i][j];
    }
  }

  return out;
}

template <typename T>
std::shared_ptr<Buffer<T>> manage_host_buffer(const size_t n, T *data) {
  return Buffer<T>::own(n, data, &free, MEMORY_SPACE_HOST);
}

template <typename T>
std::shared_ptr<Buffer<T *>> manage_host_buffer(const size_t n0,
                                                const size_t n1, T **data) {
  auto ret = std::make_shared<Buffer<T *>>(
      n0, n1, data,
      [=](int n, void **x) {
        for (int i = 0; i < n; ++i) {
          free(x[i]);
        }
        free(x);
      },
      MEMORY_SPACE_HOST);
  return ret;
}

template <typename T>
std::shared_ptr<Buffer<T>> view(const std::shared_ptr<Buffer<T>> &buffer,
                                const size_t begin, const size_t end) {
  return std::make_shared<Buffer<T>>(
      end - begin, const_cast<T *>(&buffer->data()[begin]),
      [keep_alive = buffer](void *) { (void)keep_alive; }, buffer->mem_space());
}

template <typename T>
std::shared_ptr<Buffer<T *>>
view(const std::shared_ptr<Buffer<T *>> &buffer, const size_t begin0,
     const size_t end0, const size_t begin1, const size_t end1) {
  const size_t extent0 = end0 - begin0;
  const size_t extent1 = end1 - begin1;

  T **new_buffer = static_cast<T **>(malloc(extent0 * sizeof(T *)));
  for (size_t i0 = 0; i0 < extent0; i0++) {
    new_buffer[i0] = &(buffer->data()[begin0 + i0][begin1]);
  }

  return std::make_shared<Buffer<T *>>(
      extent0, extent1, new_buffer,
      [keep_alive = buffer](int, void **buff) { free(buff); },
      buffer->mem_space());
}

template <typename T>
std::shared_ptr<Buffer<T>> sub(const std::shared_ptr<Buffer<T *>> &buffer,
                               const size_t i0) {
  return std::make_shared<Buffer<T>>(
      buffer->extent(1), buffer->data()[i0], [keep_alive = buffer](void *) {},
      buffer->mem_space());
}

template <typename O>
std::shared_ptr<Buffer<O>> astype(const std::shared_ptr<Buffer<O>> &in) {
  return in;
}

template <typename O, typename I>
std::shared_ptr<Buffer<O>> astype(const std::shared_ptr<Buffer<I>> &in) {
  const size_t size = in->size();
  auto out = create_host_buffer<O>(size);

  auto din = in->data();
  auto dout = out->data();

  for (size_t i = 0; i < size; i++) {
    dout[i] = static_cast<O>(din[i]);
  }

  return out;
}

template <typename O, typename I>
std::shared_ptr<Buffer<O *>> astype(const std::shared_ptr<Buffer<I *>> &in) {
  const size_t n0 = in->extent(0);
  const size_t n1 = in->extent(1);
  auto out = create_host_buffer<O>(n0, n1);

  auto din = in->data();
  auto dout = out->data();

  for (size_t i = 0; i < n0; i++) {
    for (size_t j = 0; j < n1; j++) {
      dout[i][j] = static_cast<O>(din[i][j]);
    }
  }

  return out;
}

template <typename T>
std::shared_ptr<Buffer<T *>> copy(const std::shared_ptr<Buffer<T *>> &buffer) {
  if (buffer->mem_space() == MEMORY_SPACE_DEVICE) {
    SMESH_IMPLEMENT_ME();
  }

  auto data = static_cast<T **>(malloc(buffer->extent(0) * sizeof(T *)));
  for (size_t i = 0; i < buffer->extent(0); i++) {
    data[i] = static_cast<T *>(malloc(buffer->extent(1) * sizeof(T)));
    std::memcpy(data[i], buffer->data()[i], buffer->extent(1) * sizeof(T));
  }

  return manage_host_buffer<T>(buffer->extent(0), buffer->extent(1), data);
}

template <typename T>
std::shared_ptr<Buffer<T *>>
zeros_like(const std::shared_ptr<Buffer<T *>> &buffer) {
  if (buffer->mem_space() == MEMORY_SPACE_DEVICE) {
    SMESH_IMPLEMENT_ME();
  }

  auto data = static_cast<T **>(malloc(buffer->extent(0) * sizeof(T *)));
  for (size_t i = 0; i < buffer->extent(0); i++) {
    data[i] = static_cast<T *>(calloc(buffer->extent(1), sizeof(T)));
  }

  return manage_host_buffer<T>(buffer->extent(0), buffer->extent(1), data);
}

} // namespace smesh

#endif // SMESH_BUFFER_IMPL_HPP