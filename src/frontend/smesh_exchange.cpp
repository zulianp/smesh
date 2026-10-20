#include "smesh_exchange.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_mesh.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <vector>

#if defined(SMESH_ENABLE_MPI)
#include "smesh_alltoallv.impl.hpp"
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
  ptrdiff_t nnodes{0};
  ptrdiff_t n_owned_nodes{0};
  SharedBuffer<i64> send_count;
  SharedBuffer<i64> send_displs;
  SharedBuffer<i64> recv_count;
  SharedBuffer<i64> recv_displs;
  SharedBuffer<idx_t> sparse_idx;
  SharedBuffer<idx_t> import_idx;
  SharedBuffer<char> send_buffer;
  SharedBuffer<char> recv_buffer;

#if defined(SMESH_ENABLE_MPI)
  // A gather is point-to-point with the ranks that take part, nobody else: a collective over
  // the communicator would touch every rank on every gather. Its messages travel on a
  // duplicate of the communicator, where no other traffic can match them.
  struct Peer {
    int rank;
    i64 offset; // in nodes, into the packed sends or the received values
    i64 count;  // nodes
  };
  std::vector<Peer> gather_sends;
  std::vector<Peer> gather_recvs;
  std::vector<MPI_Request> gather_requests;
  MPI_Comm gather_comm{MPI_COMM_NULL};
  bool gather_pending{false};
  std::function<int()> gather_finish;
#endif

  Impl(const std::shared_ptr<Communicator> &comm) : comm(comm) {}

  ~Impl() {
#if defined(SMESH_ENABLE_MPI)
    if (gather_pending && gather_finish) {
      (void)gather_finish();
    }
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (gather_comm != MPI_COMM_NULL && !finalized) {
      MPI_Comm_free(&gather_comm);
    }
#endif
  }
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
  auto ret = create(mesh->comm(), exchange_scope, n_local_nodes,
                    dist->n_nodes_owned(), dist->node_owner()->data(),
                    dist->node_offsets()->data(), import_idx);

  // The neighbours' requests (sparse_idx: the nodes this rank sends) name its nodes by
  // the positions they were created at; a renumbering of the owned-not-shared prefix
  // (packing) moved those, which the aura reaches. Send from where each is now.
  auto renumbering = dist->owned_renumbering();
  if (ret && renumbering && ret->impl_->sparse_idx) {
    const ptrdiff_t n_moved = static_cast<ptrdiff_t>(renumbering->size());
    const idx_t *const ren = renumbering->data();
    idx_t *const sends = ret->impl_->sparse_idx->data();
    const ptrdiff_t n_sends = static_cast<ptrdiff_t>(ret->impl_->sparse_idx->size());
    for (ptrdiff_t i = 0; i < n_sends; i++) {
      if (sends[i] < n_moved) sends[i] = ren[sends[i]];
    }
  }
  return ret;

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

  // A gather sends what the others asked for (recv_count, from sparse_idx) and receives
  // what this rank asked for (send_count, into ghosts and aura).
  const i64 *const send_count = ret->impl_->send_count->data();
  const i64 *const recv_count = ret->impl_->recv_count->data();
  for (int r = 0; r < size; r++) {
    if (recv_count[r] > 0) {
      ret->impl_->gather_sends.push_back({r, recv_displs[r], recv_count[r]});
    }
    if (send_count[r] > 0) {
      ret->impl_->gather_recvs.push_back({r, send_displs[r], send_count[r]});
    }
  }
  ret->impl_->gather_requests.reserve(ret->impl_->gather_sends.size() +
                                      ret->impl_->gather_recvs.size());
  SMESH_MPI_CATCH(MPI_Comm_dup(comm->get(), &ret->impl_->gather_comm));

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
  if (gather_begin(inout, block_size) != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }
  return gather_wait();
}

template <typename T>
int Exchange::gather_begin(T *const inout, const ptrdiff_t block_size) {
#if defined(SMESH_ENABLE_MPI)
  SMESH_TRACE_SCOPE("Exchange::gather_begin");
  if (impl_->gather_pending) {
    SMESH_ERROR("Exchange::gather_begin: gather already pending\n");
    return SMESH_FAILURE;
  }

  const int size = impl_->comm->size();
  const ptrdiff_t send_total = impl_->send_displs->data()[size];
  const ptrdiff_t recv_total = impl_->recv_displs->data()[size];
  const bool ghosts_only = impl_->exchange_scope == ExchangeScope::GhostsOnly;

  // Pack: recv_count/displs are used as send side for gather (see exchange_gather).
  ensure_exchange_buffer<T>(impl_->send_buffer, recv_total * block_size);
  if (!ghosts_only) {
    ensure_exchange_buffer<T>(impl_->recv_buffer, send_total * block_size);
  }

  T *const send_buf = (T *)impl_->send_buffer->data();
  T *const recv_buf =
      ghosts_only ? nullptr : (T *)impl_->recv_buffer->data();

  const idx_t *const gather_idx = impl_->sparse_idx->data();
  // One thread: a gather packs a few ten thousand values, too few to repay a parallel region
  // (measured on Grace, 4 ranks x 72 threads and 8 x 36 pack faster serially).
  for (ptrdiff_t i = 0; i < recv_total; i++) {
    const ptrdiff_t src = (ptrdiff_t)gather_idx[i] * block_size;
    const ptrdiff_t dst = i * block_size;
    for (ptrdiff_t c = 0; c < block_size; ++c) {
      send_buf[dst + c] = inout[src + c];
    }
  }

  T *const recv_ptr =
      ghosts_only ? &inout[impl_->n_owned_nodes * block_size] : recv_buf;

  // Counted in elements of T, in messages of at most INT_MAX of them.
  auto post = [&](const bool send, T *const buf, const i64 n, const int rank) -> int {
    const i64 chunk = std::numeric_limits<int>::max();
    for (i64 b = 0; b < n; b += chunk) {
      const int m = (int)std::min(chunk, n - b);
      MPI_Request request;
      if (send) {
        SMESH_MPI_CATCH(MPI_Isend(buf + b, m, mpi_type<T>(), rank, 0,
                                  impl_->gather_comm, &request));
      } else {
        SMESH_MPI_CATCH(MPI_Irecv(buf + b, m, mpi_type<T>(), rank, 0,
                                  impl_->gather_comm, &request));
      }
      impl_->gather_requests.push_back(request);
    }
    return SMESH_SUCCESS;
  };
  impl_->gather_requests.clear();
  for (const auto &peer : impl_->gather_recvs) {
    if (post(false, recv_ptr + peer.offset * block_size, peer.count * block_size,
             peer.rank) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }
  for (const auto &peer : impl_->gather_sends) {
    if (post(true, send_buf + peer.offset * block_size, peer.count * block_size,
             peer.rank) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }

  impl_->gather_pending = true;
  const ptrdiff_t bs = block_size;
  const bool use_aura = !ghosts_only;
  impl_->gather_finish = [this, inout, bs, use_aura, recv_buf]() -> int {
    SMESH_MPI_CATCH(MPI_Waitall((int)impl_->gather_requests.size(),
                                impl_->gather_requests.data(),
                                MPI_STATUSES_IGNORE));
    impl_->gather_requests.clear();
    if (use_aura) {
      const int size_local = impl_->comm->size();
      const ptrdiff_t n_import = impl_->send_displs->data()[size_local];
      const idx_t *const import_idx = impl_->import_idx->data();
      for (ptrdiff_t i = 0; i < n_import; ++i) {
        const ptrdiff_t dst = (ptrdiff_t)import_idx[i] * bs;
        const ptrdiff_t src = i * bs;
        for (ptrdiff_t c = 0; c < bs; ++c) {
          inout[dst + c] = recv_buf[src + c];
        }
      }
    }
    impl_->gather_pending = false;
    impl_->gather_finish = nullptr;
    return SMESH_SUCCESS;
  };

  return SMESH_SUCCESS;
#else
  SMESH_UNUSED(inout);
  SMESH_UNUSED(block_size);
  return SMESH_SUCCESS;
#endif
}

int Exchange::gather_wait() {
#if defined(SMESH_ENABLE_MPI)
  SMESH_TRACE_SCOPE("Exchange::gather_wait");
  if (!impl_->gather_pending) {
    return SMESH_SUCCESS;
  }
  if (!impl_->gather_finish) {
    SMESH_ERROR("Exchange::gather_wait: missing finish callback\n");
    return SMESH_FAILURE;
  }
  return impl_->gather_finish();
#else
  return SMESH_SUCCESS;
#endif
}

template int Exchange::gather_begin<f16>(f16 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<f32>(f32 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<f64>(f64 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<i8>(i8 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::gather_begin<i16>(i16 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<i32>(i32 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<i64>(i64 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<u8>(u8 *const inout,
                                        const ptrdiff_t block_size);
template int Exchange::gather_begin<u16>(u16 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<u32>(u32 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<u64>(u64 *const inout,
                                         const ptrdiff_t block_size);
template int Exchange::gather_begin<char>(char *const inout,
                                          const ptrdiff_t block_size);
#if defined(__clang__)
template int Exchange::gather_begin<long>(long *const inout,
                                          const ptrdiff_t block_size);
#endif

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
