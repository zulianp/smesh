#include "smesh_exchange.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_mesh.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#if defined(SMESH_ENABLE_MPI)
#include "smesh_alltoallv.impl.hpp"
#include "smesh_distributed_aura.hpp"
#include <mpi.h>
#if __has_include(<mpi-ext.h>)
#include <mpi-ext.h>
#endif
#endif

#if defined(SMESH_ENABLE_CUDA)
#include "smesh_device_buffer.hpp"
#include "smesh_exchange_kernels.hpp"
#if defined(SMESH_ENABLE_MPI)
#include "smesh_device_communicator.hpp"
#endif
#endif

namespace smesh {

template <typename T>
static void ensure_exchange_buffer(SharedBuffer<char> &buffer,
                                   const ptrdiff_t entries,
                                   const ExecutionSpace es = EXECUTION_SPACE_HOST) {
  const size_t nbytes = (size_t)entries * sizeof(T);
  if (!buffer || buffer->size() < nbytes) {
#if defined(SMESH_ENABLE_CUDA)
    if (es == EXECUTION_SPACE_DEVICE) {
      buffer = create_device_buffer<char>(nbytes);
      return;
    }
#endif
    SMESH_UNUSED(es);
    buffer = create_host_buffer<char>(nbytes);
  }
}

#if defined(SMESH_ENABLE_MPI) && defined(SMESH_ENABLE_CUDA)
// Whether this MPI moves device memory, at run time.
static bool mpi_takes_device_memory() {
#if defined(CRAY_MPICH_VERSION)
  // Cray MPICH's MPIX queries answer 0 whether or not its GPU support is on (8.1.30 on Alps,
  // where device-to-device messages then arrive intact); the support is its environment's.
  const char *const enabled = std::getenv("MPICH_GPU_SUPPORT_ENABLED");
  return enabled && std::string(enabled) == "1";
#elif defined(MPIX_GPU_SUPPORT_CUDA)
  int supported = 0;
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA, &supported);
  return supported != 0;
#elif defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  return MPIX_Query_cuda_support() != 0;
#else
  return false;
#endif
}
#endif

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
  // Where gathers' buffers live; on the device, sparse_idx and import_idx copied there.
  ExecutionSpace execution_space{EXECUTION_SPACE_HOST};
  SharedBuffer<idx_t> device_sparse_idx;
  SharedBuffer<idx_t> device_import_idx;

#if defined(SMESH_ENABLE_MPI)
  // A gather is point-to-point with the ranks that take part, nobody else: a collective over
  // the communicator would touch every rank on every gather. Its messages travel on a
  // duplicate of the communicator, where no other traffic can match them.
  struct Peer {
    int rank;
    i64 offset; // in nodes, into the packed sends or the received values
    i64 count;  // nodes
    bool on_node{false}; // a device exchange's peer on this node: through NCCL
  };
  std::vector<Peer> gather_sends;
  std::vector<Peer> gather_recvs;
  std::vector<MPI_Request> gather_requests;
  MPI_Comm gather_comm{MPI_COMM_NULL};
  bool gather_pending{false};
#if defined(SMESH_ENABLE_CUDA)
  // With the communicator's device communication enabled, a device exchange's messages to peers
  // on this node by NCCL, queued on the stream; to the others (any_remote) by GPU-aware MPI.
  std::shared_ptr<DeviceCommunicator> device_comm;
  bool any_remote{true};
#endif
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
                       const ExchangeScope exchange_scope,
                       const ExecutionSpace execution_space) {
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

  if (ret && execution_space == EXECUTION_SPACE_DEVICE) {
#if defined(SMESH_ENABLE_CUDA)
    // Where the run enabled it, peers on this node through NCCL and the others through
    // GPU-aware MPI: NCCL is the faster within a node and not to be trusted between nodes
    // (smesh_device_communicator.hpp).
    auto &impl = *ret->impl_;
    impl.device_comm = mesh->comm()->device();
    if (impl.device_comm) {
      impl.any_remote = false;
      for (auto *peers : {&impl.gather_sends, &impl.gather_recvs}) {
        for (auto &peer : *peers) {
          peer.on_node = impl.device_comm->node_rank(peer.rank) >= 0;
          impl.any_remote = impl.any_remote || !peer.on_node;
        }
      }
    }
    if (impl.any_remote && !mpi_takes_device_memory()) {
      SMESH_ERROR("Exchange::create_nodal: a device exchange with ranks on other nodes needs "
                  "GPU-aware MPI, which this MPI does not provide (with Cray MPICH: "
                  "MPICH_GPU_SUPPORT_ENABLED=1)\n");
    }
    ret->impl_->execution_space = EXECUTION_SPACE_DEVICE;
    ret->impl_->device_sparse_idx = to_device(ret->impl_->sparse_idx);
    if (ret->impl_->import_idx) {
      ret->impl_->device_import_idx = to_device(ret->impl_->import_idx);
    }
#else
    SMESH_ERROR("Exchange::create_nodal: a device exchange in a build without CUDA\n");
#endif
  }
  return ret;

#else
  SMESH_UNUSED(exchange_scope);
  SMESH_UNUSED(execution_space);
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
  if (impl_->execution_space != EXECUTION_SPACE_HOST) {
    SMESH_ERROR("Exchange::scatter_add: host buffers only\n");
  }
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
int Exchange::gather_begin(T *const inout, const ptrdiff_t block_size, void *const stream) {
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

  const ExecutionSpace es = impl_->execution_space;
  // Pack: recv_count/displs are used as send side for gather (see exchange_gather).
  ensure_exchange_buffer<T>(impl_->send_buffer, recv_total * block_size, es);
  if (!ghosts_only) {
    ensure_exchange_buffer<T>(impl_->recv_buffer, send_total * block_size, es);
  }

  T *const send_buf = (T *)impl_->send_buffer->data();
  T *const recv_buf =
      ghosts_only ? nullptr : (T *)impl_->recv_buffer->data();

  const idx_t *const gather_idx = impl_->sparse_idx->data();
#if defined(SMESH_ENABLE_CUDA)
  if (es == EXECUTION_SPACE_DEVICE) {
    // Packed on the device, and finished before MPI reads the buffer; NCCL's messages follow
    // the packing on the stream.
    if (device::gather_by_index(recv_total, block_size, sizeof(T),
                                impl_->device_sparse_idx->data(), inout,
                                send_buf, stream) != SMESH_SUCCESS ||
        (impl_->any_remote && device::synchronize(stream) != SMESH_SUCCESS)) {
      return SMESH_FAILURE;
    }
  } else
#endif
  {
    // One thread: a gather packs a few ten thousand values, too few to repay a parallel region
    // (measured on Grace, 4 ranks x 72 threads and 8 x 36 pack faster serially).
    for (ptrdiff_t i = 0; i < recv_total; i++) {
      const ptrdiff_t src = (ptrdiff_t)gather_idx[i] * block_size;
      const ptrdiff_t dst = i * block_size;
      for (ptrdiff_t c = 0; c < block_size; ++c) {
        send_buf[dst + c] = inout[src + c];
      }
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
#if defined(SMESH_ENABLE_CUDA)
  if (impl_->device_comm) {
    const auto &device = *impl_->device_comm;
    device.group_begin();
    for (const auto &peer : impl_->gather_recvs) {
      if (peer.on_node) {
        device.recv(recv_ptr + peer.offset * block_size,
                    (size_t)(peer.count * block_size) * sizeof(T), peer.rank, stream);
      }
    }
    for (const auto &peer : impl_->gather_sends) {
      if (peer.on_node) {
        device.send(send_buf + peer.offset * block_size,
                    (size_t)(peer.count * block_size) * sizeof(T), peer.rank, stream);
      }
    }
    device.group_end();
  }
#endif
  impl_->gather_requests.clear();
  for (const auto &peer : impl_->gather_recvs) {
    if (peer.on_node) continue;
    if (post(false, recv_ptr + peer.offset * block_size, peer.count * block_size,
             peer.rank) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }
  for (const auto &peer : impl_->gather_sends) {
    if (peer.on_node) continue;
    if (post(true, send_buf + peer.offset * block_size, peer.count * block_size,
             peer.rank) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }

  impl_->gather_pending = true;
  const ptrdiff_t bs = block_size;
  const bool use_aura = !ghosts_only;
  impl_->gather_finish = [this, inout, bs, use_aura, recv_buf, stream]() -> int {
    SMESH_UNUSED(stream);  // the device unpack's
    SMESH_MPI_CATCH(MPI_Waitall((int)impl_->gather_requests.size(),
                                impl_->gather_requests.data(),
                                MPI_STATUSES_IGNORE));
    impl_->gather_requests.clear();
    if (use_aura) {
      const int size_local = impl_->comm->size();
      const ptrdiff_t n_import = impl_->send_displs->data()[size_local];
#if defined(SMESH_ENABLE_CUDA)
      if (impl_->execution_space == EXECUTION_SPACE_DEVICE) {
        if (device::scatter_by_index(n_import, bs, sizeof(T),
                                     impl_->device_import_idx->data(), recv_buf,
                                     inout, stream) != SMESH_SUCCESS) {
          return SMESH_FAILURE;
        }
        impl_->gather_pending = false;
        impl_->gather_finish = nullptr;
        return SMESH_SUCCESS;
      }
#endif
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
  SMESH_UNUSED(stream);
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
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<f32>(f32 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<f64>(f64 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<i8>(i8 *const inout,
                                        const ptrdiff_t block_size,
                                        void *const stream);
template int Exchange::gather_begin<i16>(i16 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<i32>(i32 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<i64>(i64 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<u8>(u8 *const inout,
                                        const ptrdiff_t block_size,
                                        void *const stream);
template int Exchange::gather_begin<u16>(u16 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<u32>(u32 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<u64>(u64 *const inout,
                                         const ptrdiff_t block_size,
                                         void *const stream);
template int Exchange::gather_begin<char>(char *const inout,
                                          const ptrdiff_t block_size,
                                          void *const stream);
#if defined(__clang__)
template int Exchange::gather_begin<long>(long *const inout,
                                          const ptrdiff_t block_size,
                                          void *const stream);
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
