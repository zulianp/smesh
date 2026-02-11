#ifndef SMESH_BUFFER_HPP
#define SMESH_BUFFER_HPP

#include <cassert>
#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>

#include "smesh_base.hpp"
#include "smesh_read.hpp"
#include "smesh_spaces.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"

namespace smesh {

class BaseBuffer {
public:
  size_t n_;
  PrimitiveType type_;
  void *ptr_;
  std::function<void(void *)> destroy_;
  MemorySpace mem_space_;
  BaseBuffer(const size_t n, const PrimitiveType type, void *const ptr,
             std::function<void(void *)> destroy, MemorySpace mem_space);

  virtual ~BaseBuffer();
  void *void_data();
  const void *void_data() const;

  inline PrimitiveType element_type() const { return type_; }
  inline size_t size() const { return n_; }
  inline size_t nbytes() const { return n_ * num_bytes(type_); }
  inline MemorySpace mem_space() const { return mem_space_; }
};

template <typename T> class Buffer final : public BaseBuffer {
public:
  Buffer(const size_t n, T *const ptr, std::function<void(void *)> destroy,
         MemorySpace mem_space);

  Buffer();

  T *data();
  const T *data() const;

  void print(std::ostream &os = std::cout) const;

  static std::shared_ptr<Buffer<T>>
  wrap(const size_t n, T *const x,
       enum MemorySpace mem_space = MEMORY_SPACE_INVALID);

  static std::shared_ptr<Buffer<T>>
  own(const size_t n, T *x, std::function<void(void *)> destroy,
      enum MemorySpace mem_space = MEMORY_SPACE_INVALID);

  static std::shared_ptr<Buffer<T>> make_empty();

  int to_file(const Path &path) const;
};

template <typename T> class Buffer<T *> {
public:
  Buffer(const size_t n0, const size_t n1, T **const ptr,
         std::function<void(int n, void **)> destroy, MemorySpace mem_space)
      : extent_{n0, n1}, ptr_(ptr), destroy_(destroy), mem_space_(mem_space) {}

  ~Buffer();

  T **data();
  const T **data() const;
  inline size_t extent(int i) const { return extent_[i]; }
  inline MemorySpace mem_space() const { return mem_space_; }

  inline size_t nbytes() const { return extent_[0] * extent_[1] * sizeof(T); }

  void print(std::ostream &os = std::cout);

  static std::shared_ptr<Buffer<T *>>
  wrap(const size_t n0, const size_t n1, T **x,
       enum MemorySpace mem_space = MEMORY_SPACE_INVALID);

  static std::shared_ptr<Buffer<T *>>
  own(const size_t n0, const size_t n1, T **x,
      std::function<void(int n, void **)> destroy,
      enum MemorySpace mem_space = MEMORY_SPACE_INVALID);

  int to_files(const Path &path);

  void release();

private:
  size_t extent_[2];
  T **ptr_{nullptr};
  std::function<void(int n, void **)> destroy_;
  MemorySpace mem_space_;
};

template <typename T>
std::shared_ptr<Buffer<T>> create_host_buffer(const size_t n);

template <typename T>
std::shared_ptr<Buffer<T *>> create_host_buffer(const size_t n0,
                                                const size_t n1);

template <typename T>
std::shared_ptr<Buffer<T *>> create_host_buffer_fake_SoA(const size_t n0,
                                                         const size_t n1);

template <typename T>
std::shared_ptr<Buffer<T *>>
convert_host_buffer_to_fake_SoA(const size_t n0,
                                const std::shared_ptr<Buffer<T>> &in);

template <typename T>
static std::shared_ptr<Buffer<T>>
soa_to_aos(const size_t in_stride0, const size_t in_stride1,
           const std::shared_ptr<Buffer<T *>> &in);

template <typename T>
std::shared_ptr<Buffer<T>> manage_host_buffer(const size_t n, T *data);

template <typename T>
std::shared_ptr<Buffer<T *>> manage_host_buffer(const size_t n0,
                                                const size_t n1, T **data);

template <typename T>
std::shared_ptr<Buffer<T>> view(const std::shared_ptr<Buffer<T>> &buffer,
                                const size_t begin, const size_t end);

template <typename T>
std::shared_ptr<Buffer<T *>>
view(const std::shared_ptr<Buffer<T *>> &buffer, const size_t begin0,
     const size_t end0, const size_t begin1, const size_t end1);

template <typename T>
std::shared_ptr<Buffer<T>> sub(const std::shared_ptr<Buffer<T *>> &buffer,
                               const size_t i0);

template <typename O>
std::shared_ptr<Buffer<O>> astype(const std::shared_ptr<Buffer<O>> &in);

template <typename O, typename I>
std::shared_ptr<Buffer<O>> astype(const std::shared_ptr<Buffer<I>> &in);

template <typename O, typename I>
std::shared_ptr<Buffer<O *>> astype(const std::shared_ptr<Buffer<I *>> &in);

template <typename T>
std::shared_ptr<Buffer<T *>> copy(const std::shared_ptr<Buffer<T *>> &buffer);

template <typename T>
std::shared_ptr<Buffer<T *>>
zeros_like(const std::shared_ptr<Buffer<T *>> &buffer);

template <typename T> using SharedBuffer = std::shared_ptr<Buffer<T>>;

} // namespace smesh

#endif // SMESH_BUFFER_HPP
