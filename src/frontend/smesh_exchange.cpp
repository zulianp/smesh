#include "smesh_exchange.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_mesh.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"

#include <algorithm>

#if defined(SMESH_ENABLE_MPI)
#include "smesh_distributed_aura.hpp"
#include <mpi.h>
#endif

namespace smesh {

template <typename T>
static void ensure_exchange_buffer(SharedBuffer<char> &buffer,
                                   const ptrdiff_t entries) {
  const size_t nbytes = (size_t)entries * sizeof(T);
  if (!buffer || buffer->size() < nbytes) {
    buffer = create_host_buffer<char>(nbytes);
  }
}

class Exchange::Impl {
public:
  ExchangeScope exchange_scope = ExchangeScope::GhostsOnly;
  std::shared_ptr<Communicator> comm;
  ptrdiff_t nnodes;
  ptrdiff_t n_owned_nodes;
  SharedBuffer<i64> send_count;
  SharedBuffer<i64> send_displs;
  SharedBuffer<i64> recv_count;
  SharedBuffer<i64> recv_displs;
  SharedBuffer<idx_t> sparse_idx;
  SharedBuffer<idx_t> import_idx;
  SharedBuffer<char> send_buffer;
  SharedBuffer<char> recv_buffer;
  Impl(const std::shared_ptr<Communicator> &comm) : comm(comm) {}
};

Exchange::Exchange(const std::shared_ptr<Communicator> &comm)
    : impl_(std::make_unique<Impl>(comm)) {}

Exchange::~Exchange() = default;

std::shared_ptr<Exchange>
Exchange::create_nodal(const std::shared_ptr<Mesh> &mesh,
                       const ExchangeScope exchange_scope) {
#if defined(SMESH_ENABLE_MPI)
  auto dist = mesh->distributed();
  const bool with_aura = exchange_scope == ExchangeScope::GhostsAndAura;
  const ptrdiff_t n_local_nodes = dist->n_nodes_owned() +
                                  dist->n_nodes_ghosts() +
                                  (with_aura ? dist->n_nodes_aura() : 0);
  const idx_t *const import_idx =
      with_aura ? dist->ghosts_and_aura()->data() : dist->ghosts()->data();
  return create(mesh->comm(), exchange_scope, n_local_nodes,
                dist->n_nodes_owned(), dist->node_owner()->data(),
                dist->node_offsets()->data(), import_idx);

#else
  SMESH_UNUSED(exchange_scope);
  return std::make_shared<Exchange>(mesh->comm());
#endif
}

#if defined(SMESH_ENABLE_MPI)
std::shared_ptr<Exchange> Exchange::create(
    const std::shared_ptr<Communicator> &comm,
    const ExchangeScope exchange_scope, const ptrdiff_t n_local_nodes,
    const ptrdiff_t n_owned_nodes, const int *const node_owner,
    const ptrdiff_t *const node_offsets, const idx_t *const ghosts) {
  int size = comm->size();
  auto ret = std::make_shared<Exchange>(comm);
  ret->impl_->exchange_scope = exchange_scope;
  ret->impl_->nnodes = n_local_nodes;
  ret->impl_->n_owned_nodes = n_owned_nodes;
  ret->impl_->send_count = create_host_buffer<i64>(size);
  ret->impl_->send_displs = create_host_buffer<i64>(size + 1);
  ret->impl_->recv_count = create_host_buffer<i64>(size);
  ret->impl_->recv_displs = create_host_buffer<i64>(size + 1);

  idx_t *sparse_idx = nullptr;
  idx_t *import_idx = nullptr;
  SMESH_ASSERT(n_local_nodes >= n_owned_nodes);
  if (exchange_scope == ExchangeScope::GhostsAndAura) {
    exchange_create(comm->get(), n_local_nodes, n_owned_nodes, node_owner,
                    node_offsets, ghosts, ret->impl_->send_count->data(),
                    ret->impl_->send_displs->data(),
                    ret->impl_->recv_count->data(),
                    ret->impl_->recv_displs->data(), &sparse_idx, &import_idx);
  } else {
    exchange_create_ghosts(
        comm->get(), n_local_nodes, n_owned_nodes, node_owner, node_offsets,
        ghosts, ret->impl_->send_count->data(), ret->impl_->send_displs->data(),
        ret->impl_->recv_count->data(), ret->impl_->recv_displs->data(),
        &sparse_idx);
  }

  auto send_displs = ret->impl_->send_displs->data();
  auto recv_displs = ret->impl_->recv_displs->data();
  const ptrdiff_t send_total = send_displs[size];

  ret->impl_->sparse_idx =
      manage_host_buffer<idx_t>((ptrdiff_t)recv_displs[size], sparse_idx);
  if (exchange_scope == ExchangeScope::GhostsAndAura) {
    ret->impl_->import_idx = manage_host_buffer<idx_t>(send_total, import_idx);
  }

  return ret;
}
#endif

template <typename T> int Exchange::scatter_add(T *const inout) {
  return scatter_add(inout, 1);
}

template <typename T>
int Exchange::scatter_add(T *const inout, const ptrdiff_t block_size) {
#if defined(SMESH_ENABLE_MPI)
  SMESH_TRACE_SCOPE("Exchange::scatter_add");
  const int size = impl_->comm->size();
  const ptrdiff_t send_total = impl_->send_displs->data()[size];
  const ptrdiff_t recv_total = impl_->recv_displs->data()[size];
  if (impl_->exchange_scope == ExchangeScope::GhostsOnly) {
    ensure_exchange_buffer<T>(impl_->recv_buffer, recv_total * block_size);
    return exchange_scatter_add_ghosts(
        impl_->comm->get(), impl_->n_owned_nodes, impl_->send_count->data(),
        impl_->send_displs->data(), impl_->recv_count->data(),
        impl_->recv_displs->data(), impl_->sparse_idx->data(), inout,
        (T *)impl_->recv_buffer->data(), block_size);
  }
  ensure_exchange_buffer<T>(impl_->send_buffer, send_total * block_size);
  ensure_exchange_buffer<T>(impl_->recv_buffer, recv_total * block_size);
  return exchange_scatter_add(
      impl_->comm->get(), impl_->n_owned_nodes, impl_->send_count->data(),
      impl_->send_displs->data(), impl_->recv_count->data(),
      impl_->recv_displs->data(), impl_->sparse_idx->data(),
      impl_->import_idx->data(), inout, (T *)impl_->send_buffer->data(),
      (T *)impl_->recv_buffer->data(), block_size);
#else
  SMESH_UNUSED(inout);
  SMESH_UNUSED(block_size);
  return SMESH_SUCCESS;
#endif
}

template <typename T> int Exchange::gather(T *const inout) {
  return gather(inout, 1);
}

template <typename T>
int Exchange::gather(T *const inout, const ptrdiff_t block_size) {
#if defined(SMESH_ENABLE_MPI)
  SMESH_TRACE_SCOPE("Exchange::gather");
  const int size = impl_->comm->size();
  const ptrdiff_t send_total = impl_->send_displs->data()[size];
  const ptrdiff_t recv_total = impl_->recv_displs->data()[size];
  if (impl_->exchange_scope == ExchangeScope::GhostsOnly) {
    ensure_exchange_buffer<T>(impl_->send_buffer, recv_total * block_size);
    return exchange_gather_ghosts(
        impl_->comm->get(), impl_->n_owned_nodes, impl_->send_count->data(),
        impl_->send_displs->data(), impl_->recv_count->data(),
        impl_->recv_displs->data(), impl_->sparse_idx->data(), inout,
        (T *)impl_->send_buffer->data(), block_size);
  }
  ensure_exchange_buffer<T>(impl_->send_buffer, recv_total * block_size);
  ensure_exchange_buffer<T>(impl_->recv_buffer, send_total * block_size);
  return exchange_gather(impl_->comm->get(), impl_->n_owned_nodes,
                         impl_->send_count->data(), impl_->send_displs->data(),
                         impl_->recv_count->data(), impl_->recv_displs->data(),
                         impl_->sparse_idx->data(), impl_->import_idx->data(),
                         inout, (T *)impl_->send_buffer->data(),
                         (T *)impl_->recv_buffer->data(), block_size);
#else
  SMESH_UNUSED(inout);
  SMESH_UNUSED(block_size);
  return SMESH_SUCCESS;
#endif
}

template int Exchange::gather<f16>(f16 *const inout);
template int Exchange::gather<f16>(f16 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<f32>(f32 *const inout);
template int Exchange::gather<f32>(f32 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<f64>(f64 *const inout);
template int Exchange::gather<f64>(f64 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<i8>(i8 *const inout);
template int Exchange::gather<i8>(i8 *const inout, const ptrdiff_t block_size);
template int Exchange::gather<i16>(i16 *const inout);
template int Exchange::gather<i16>(i16 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<i32>(i32 *const inout);
template int Exchange::gather<i32>(i32 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<i64>(i64 *const inout);
template int Exchange::gather<i64>(i64 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<u8>(u8 *const inout);
template int Exchange::gather<u8>(u8 *const inout, const ptrdiff_t block_size);
template int Exchange::gather<u16>(u16 *const inout);
template int Exchange::gather<u16>(u16 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<u32>(u32 *const inout);
template int Exchange::gather<u32>(u32 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<u64>(u64 *const inout);
template int Exchange::gather<u64>(u64 *const inout,
                                   const ptrdiff_t block_size);
template int Exchange::gather<char>(char *const inout);
template int Exchange::gather<char>(char *const inout,
                                    const ptrdiff_t block_size);

#if defined(__clang__)
template int Exchange::gather<long>(long *const inout);
template int Exchange::gather<long>(long *const inout,
                                    const ptrdiff_t block_size);
#endif

template int Exchange::scatter_add<f16>(f16 *const inout);
template int Exchange::scatter_add<f16>(f16 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<f32>(f32 *const inout);
template int Exchange::scatter_add<f32>(f32 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<f64>(f64 *const inout);
template int Exchange::scatter_add<f64>(f64 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<i8>(i8 *const inout);
template int Exchange::scatter_add<i8>(i8 *const inout,
                                       const ptrdiff_t block_size);
template int Exchange::scatter_add<i16>(i16 *const inout);
template int Exchange::scatter_add<i16>(i16 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<i32>(i32 *const inout);
template int Exchange::scatter_add<i32>(i32 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<i64>(i64 *const inout);
template int Exchange::scatter_add<i64>(i64 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<u8>(u8 *const inout);
template int Exchange::scatter_add<u8>(u8 *const inout,
                                       const ptrdiff_t block_size);
template int Exchange::scatter_add<u16>(u16 *const inout);
template int Exchange::scatter_add<u16>(u16 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<u32>(u32 *const inout);
template int Exchange::scatter_add<u32>(u32 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<u64>(u64 *const inout);
template int Exchange::scatter_add<u64>(u64 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::scatter_add<char>(char *const inout);
template int Exchange::scatter_add<char>(char *const inout,
                                         const ptrdiff_t block_size);

#if defined(__clang__)
template int Exchange::scatter_add<long>(long *const inout);
template int Exchange::scatter_add<long>(long *const inout,
                                         const ptrdiff_t block_size);
#endif

} // namespace smesh
