#include "smesh_exchange.hpp"
#include "smesh_mesh.hpp"
#include "smesh_communicator.hpp"
#include "smesh_types.hpp"
#include "smesh_buffer.hpp"

#if defined(SMESH_ENABLE_MPI)
#include "smesh_distributed_aura.hpp"
#include <mpi.h>
#endif

namespace smesh {

class Exchange::Impl {
public:
  std::shared_ptr<Communicator> comm;
  ptrdiff_t nnodes;
  ptrdiff_t n_owned_nodes;
  SharedBuffer<i64> send_count;
  SharedBuffer<i64> send_displs;
  SharedBuffer<i64> recv_count;
  SharedBuffer<i64> recv_displs;
  SharedBuffer<idx_t> sparse_idx;
  SharedBuffer<char> buffer;
  Impl(const std::shared_ptr<Communicator> &comm) : comm(comm) {}
};

Exchange::Exchange(const std::shared_ptr<Communicator> &comm)
    : impl_(std::make_unique<Impl>(comm)) {}

Exchange::~Exchange() = default;

std::shared_ptr<Exchange>
Exchange::create_nodal(const std::shared_ptr<Mesh> &mesh) {
#if defined(SMESH_ENABLE_MPI)
  auto dist = mesh->distributed();
  return create(mesh->comm(), dist->n_nodes_global(), dist->n_nodes_owned(),
                dist->node_owner()->data(), dist->node_offsets()->data(),
                dist->ghosts()->data());

#else
  return std::make_shared<Exchange>(mesh->comm());
#endif
}

#if defined(SMESH_ENABLE_MPI)
std::shared_ptr<Exchange>
Exchange::create(const std::shared_ptr<Communicator> &comm,
                 const ptrdiff_t nnodes, const ptrdiff_t n_owned_nodes,
                 const int *const node_owner, const ptrdiff_t *const node_offsets,
                 const idx_t *const ghosts) {
  int size = comm->size();
  auto ret = std::make_shared<Exchange>(comm);
  ret->impl_->nnodes = nnodes;
  ret->impl_->n_owned_nodes = n_owned_nodes;
  ret->impl_->send_count = create_host_buffer<i64>(size);
  ret->impl_->send_displs = create_host_buffer<i64>(size + 1);
  ret->impl_->recv_count = create_host_buffer<i64>(size);
  ret->impl_->recv_displs = create_host_buffer<i64>(size + 1);
  ret->impl_->sparse_idx = create_host_buffer<idx_t>(nnodes);

  idx_t *sparse_idx = nullptr;
  exchange_create(comm->get(), nnodes, n_owned_nodes, node_owner, node_offsets,
                  ghosts, ret->impl_->send_count->data(),
                  ret->impl_->send_displs->data(),
                  ret->impl_->recv_count->data(),
                  ret->impl_->recv_displs->data(), &sparse_idx);

  auto recv_displs = ret->impl_->recv_displs->data();
  auto recv_count = ret->impl_->recv_count->data();
  const ptrdiff_t buffer_size = recv_count[size - 1] + recv_displs[size - 1];

  ret->impl_->sparse_idx =
      manage_host_buffer<idx_t>((ptrdiff_t)recv_displs[size], sparse_idx);
  ret->impl_->buffer =
      create_host_buffer<char>(buffer_size * SIZE_LARGEST_TYPE);

  return ret;
}
#endif

template <typename T> int Exchange::exchange_add(T *const inout) {
#if defined(SMESH_ENABLE_MPI)
  return exchange_scatter_add(
      impl_->comm->get(), impl_->n_owned_nodes, impl_->send_count->data(),
      impl_->send_displs->data(), impl_->recv_count->data(),
      impl_->recv_displs->data(), impl_->sparse_idx->data(), inout,
      (T *)impl_->buffer->data());
#else
SMESH_UNUSED(inout);
  return SMESH_SUCCESS;
#endif
}

template int Exchange::exchange_add<f16>(f16 *const inout);
template int Exchange::exchange_add<f32>(f32 *const inout);
template int Exchange::exchange_add<f64>(f64 *const inout);
template int Exchange::exchange_add<i8>(i8 *const inout);
template int Exchange::exchange_add<i16>(i16 *const inout);
template int Exchange::exchange_add<i32>(i32 *const inout);
template int Exchange::exchange_add<i64>(i64 *const inout);
template int Exchange::exchange_add<u8>(u8 *const inout);
template int Exchange::exchange_add<u16>(u16 *const inout);
template int Exchange::exchange_add<u32>(u32 *const inout);
template int Exchange::exchange_add<u64>(u64 *const inout);
template int Exchange::exchange_add<char>(mask_t *const inout);


#if defined(__clang__)
template int Exchange::exchange_add<long>(long *const inout);
#endif

} // namespace smesh
