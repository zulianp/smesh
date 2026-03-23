#ifndef SFEM_RESTRICT_HPP
#define SFEM_RESTRICT_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sshex8_restriction.hpp"
#include "smesh_forward_declarations.hpp"

#ifdef SMESH_ENABLE_CUDA
#include "smesh_restriction.cuh"
#endif

#include <memory>

namespace smesh {
template <typename T> class AbstractRestrict {
public:
  virtual ~AbstractRestrict() = default;
  virtual int apply(const T *const x, T *const y) = 0;
  virtual ptrdiff_t rows() const = 0;
  virtual ptrdiff_t cols() const = 0;
  virtual ExecutionSpace execution_space() const = 0;
};

template <typename T> class Restrict final : public AbstractRestrict<T> {
public:
  Restrict(const std::shared_ptr<Mesh> &from, const std::shared_ptr<Mesh> &to,
           const ExecutionSpace es, const int block_size);

  static std::shared_ptr<Restrict> create(const std::shared_ptr<Mesh> &from,
                                          const std::shared_ptr<Mesh> &to,
                                          const ExecutionSpace es,
                                          const int block_size);

  ~Restrict();
  int apply(const T *const x, T *const y) override;
  ptrdiff_t rows() const override;
  ptrdiff_t cols() const override;
  ExecutionSpace execution_space() const override;

  const SharedBuffer<uint16_t> &element_to_node_incidence_count() const;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

template <typename T> class SurfaceRestrict final : public AbstractRestrict<T> {
public:
  SurfaceRestrict(const int from_level, const smesh::ElemType from_elem_type,
                  const ptrdiff_t from_n_nodes,
                  const SharedBuffer<idx_t *> &from_sides,
                  const SharedBuffer<uint16_t> &from_count, const int to_level,
                  const smesh::ElemType to_elem_type,
                  const ptrdiff_t to_n_nodes,
                  const SharedBuffer<idx_t *> &to_sides,
                  const ExecutionSpace es, const int block_size);

  static std::shared_ptr<SurfaceRestrict<T>>
  create(const int from_level, const smesh::ElemType from_elem_type,
         const ptrdiff_t from_n_nodes, const SharedBuffer<idx_t *> &from_sides,
         const SharedBuffer<uint16_t> &from_count, const int to_level,
         const smesh::ElemType to_elem_type, const ptrdiff_t to_n_nodes,
         const SharedBuffer<idx_t *> &to_sides, const ExecutionSpace es,
         const int block_size);

  ~SurfaceRestrict();

  int apply(const T *const x, T *const y) override;
  ptrdiff_t rows() const override;
  ptrdiff_t cols() const override;
  ExecutionSpace execution_space() const override;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif
