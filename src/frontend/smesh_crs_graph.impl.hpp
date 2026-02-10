#ifndef SMESH_CRS_GRAPH_IMPL_HPP
#define SMESH_CRS_GRAPH_IMPL_HPP

#include "smesh_crs_graph.hpp"

#include "smesh_crs_graph.hpp"
#include "smesh_graph.hpp"

namespace smesh {

template <typename count_t, typename idx_t>
class CRSGraph<count_t, idx_t>::Impl {
public:
  std::shared_ptr<Buffer<count_t>> rowptr;
  std::shared_ptr<Buffer<idx_t>> colidx;
  ~Impl() {}
};

template <typename count_t, typename idx_t>
CRSGraph<count_t, idx_t>::CRSGraph(
    const std::shared_ptr<Buffer<count_t>> &rowptr,
    const std::shared_ptr<Buffer<idx_t>> &colidx)
    : impl_(std::make_unique<Impl>()) {
  impl_->rowptr = rowptr;
  impl_->colidx = colidx;
}

template <typename count_t, typename idx_t>
void CRSGraph<count_t, idx_t>::print(std::ostream &os) const {
  auto rowptr = this->rowptr()->data();
  auto colidx = this->colidx()->data();
  const ptrdiff_t nnodes = this->n_nodes();
  const ptrdiff_t nnz = this->nnz();

  os << "CRSGraph (" << nnodes << " nodes, " << nnz << " nnz)\n";
  for (ptrdiff_t i = 0; i < nnodes; i++) {
    os << i << " (" << (rowptr[i + 1] - rowptr[i]) << "): ";
    for (int j = rowptr[i]; j < rowptr[i + 1]; j++) {
      assert(j < nnz);
      os << colidx[j] << " ";
    }
    os << "\n";
  }
}

template <typename count_t, typename idx_t>
CRSGraph<count_t, idx_t>::CRSGraph() : impl_(std::make_unique<Impl>()) {}

template <typename count_t, typename idx_t>
CRSGraph<count_t, idx_t>::~CRSGraph() = default;

template <typename count_t, typename idx_t>
ptrdiff_t CRSGraph<count_t, idx_t>::n_nodes() const {
  return rowptr()->size() - 1;
}
template <typename count_t, typename idx_t>
ptrdiff_t CRSGraph<count_t, idx_t>::nnz() const {
  return colidx()->size();
}

template <typename count_t, typename idx_t>
std::shared_ptr<Buffer<count_t>> CRSGraph<count_t, idx_t>::rowptr() const {
  return impl_->rowptr;
}
template <typename count_t, typename idx_t>
std::shared_ptr<Buffer<idx_t>> CRSGraph<count_t, idx_t>::colidx() const {
  return impl_->colidx;
}

template <typename count_t, typename idx_t>
std::shared_ptr<CRSGraph<count_t, idx_t>>
CRSGraph<count_t, idx_t>::block_to_scalar(const int block_size) {
  auto rowptr = create_host_buffer<count_t>(this->n_nodes() * block_size + 1);
  auto colidx =
      create_host_buffer<idx_t>(this->nnz() * block_size * block_size);

  crs_graph_block_to_scalar(this->n_nodes(), block_size, this->rowptr()->data(),
                            this->colidx()->data(), rowptr->data(),
                            colidx->data());

  return std::make_shared<CRSGraph<count_t, idx_t>>(rowptr, colidx);
}

} // namespace smesh

#endif // SMESH_CRS_GRAPH_IMPL_HPP
