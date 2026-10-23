#include "smesh_device_communicator.hpp"

#include "smesh_communicator.hpp"

#if defined(SMESH_ENABLE_MPI) && defined(SMESH_ENABLE_NCCL)

#include "smesh_cuda_base.cuh"
#include "smesh_tracer.hpp"

#include <mpi.h>
#include <nccl.h>

#include <pthread.h>
#include <sched.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <thread>
#include <vector>

namespace smesh {

#define SMESH_NCCL_CHECK(call)                                                   \
  do {                                                                           \
    const ncclResult_t smesh_nccl_result_ = (call);                              \
    if (smesh_nccl_result_ != ncclSuccess) {                                     \
      SMESH_ERROR("%s: %s\n", #call, ncclGetErrorString(smesh_nccl_result_));   \
    }                                                                            \
  } while (0)

class DeviceCommunicator::Impl {
public:
  ncclComm_t comm{nullptr};
  // NCCL's set-up, on a thread of its own from Communicator::enable_device() to the first use (ready()):
  // the first in a process takes about 3 s on Alps, which the caller's set-up then hides.
  std::thread init;
  ncclResult_t init_result{ncclSuccess};
  double init_seconds{0};
  int rank{0};
  int size{1};
  bool whole{true};
  // The node rank of each rank of the communicator, -1 off the node.
  std::vector<int> node_rank;
};

DeviceCommunicator::DeviceCommunicator() : impl_(std::make_unique<Impl>()) {}

// The result is not checked: at exit CUDA may have been torn down already.
DeviceCommunicator::~DeviceCommunicator() {
  if (impl_->init.joinable()) impl_->init.join();
  if (impl_->comm) (void)ncclCommDestroy(impl_->comm);
}

namespace {
// The communicator once its set-up has ended, connected for reductions: NCCL connects the ranks
// at the first one (up to 0.5 s on Alps), here rather than in the first solve. Traced: the wait
// (DeviceCommunicator::ready) and the set-up's own time (DeviceCommunicator::init).
ncclComm_t ready(DeviceCommunicator::Impl &impl) {
  if (impl.init.joinable()) {
    SMESH_TRACE_SCOPE("DeviceCommunicator::ready");
    impl.init.join();
    if (impl.init_result != ncclSuccess) {
      SMESH_ERROR("ncclCommInitRank: %s\n", ncclGetErrorString(impl.init_result));
    }
#ifdef SMESH_ENABLE_TRACE
    Tracer::instance().record_event("DeviceCommunicator::init", impl.init_seconds);
#endif
    double *warm = nullptr;
    SMESH_CUDA_CHECK(cudaMalloc(&warm, sizeof(double)));
    SMESH_CUDA_CHECK(cudaMemset(warm, 0, sizeof(double)));
    const ncclResult_t r = ncclAllReduce(warm, warm, 1, ncclFloat64, ncclSum, impl.comm, 0);
    if (r != ncclSuccess) SMESH_ERROR("ncclAllReduce: %s\n", ncclGetErrorString(r));
    SMESH_CUDA_CHECK(cudaStreamSynchronize(0));
    SMESH_CUDA_CHECK(cudaFree(warm));
  }
  return impl.comm;
}
} // namespace

int DeviceCommunicator::rank() const { return impl_->rank; }
int DeviceCommunicator::size() const { return impl_->size; }
bool DeviceCommunicator::whole() const { return impl_->whole; }

int DeviceCommunicator::node_rank(const int comm_rank) const {
  return impl_->node_rank[comm_rank];
}

namespace {
ncclDataType_t nccl_type(const enum PrimitiveType type) {
  switch (type) {
  case SMESH_FLOAT32:
    return ncclFloat32;
  case SMESH_FLOAT64:
    return ncclFloat64;
  default:
    SMESH_ERROR("DeviceCommunicator: reductions of f32 or f64\n");
    return ncclFloat64;
  }
}

size_t nccl_bytes(const enum PrimitiveType type) {
  return type == SMESH_FLOAT32 ? 4 : 8;
}

void require_whole(const bool whole) {
  if (!whole) {
    SMESH_ERROR("DeviceCommunicator: a reduction over ranks on more than one node\n");
  }
}

int peer_on_node(const std::vector<int> &node_rank, const int comm_rank) {
  const int peer = node_rank[comm_rank];
  if (peer < 0) {
    SMESH_ERROR("DeviceCommunicator: rank %d is on another node\n", comm_rank);
  }
  return peer;
}

// The CPUs next to `device` (its PCI device's local_cpulist), as NCCL places its own threads; false
// where the system does not say.
bool cpus_near(const int device, cpu_set_t *const cpus) {
  char bus[32];
  if (cudaDeviceGetPCIBusId(bus, (int)sizeof(bus), device) != cudaSuccess) return false;
  std::string id(bus);
  for (char &c : id) c = (char)std::tolower((unsigned char)c);
  FILE *const f = std::fopen(("/sys/bus/pci/devices/" + id + "/local_cpulist").c_str(), "r");
  if (!f) return false;
  CPU_ZERO(cpus);
  int first = 0, last = 0, n = 0;
  char sep = 0;
  while (std::fscanf(f, "%d", &first) == 1) {
    last = first;
    if (std::fscanf(f, "%c", &sep) == 1 && sep == '-') {
      if (std::fscanf(f, "%d", &last) != 1) break;
      if (std::fscanf(f, "%c", &sep) != 1) sep = 0;
    }
    for (int c = first; c <= last && c < CPU_SETSIZE; c++, n++) CPU_SET(c, cpus);
    if (sep != ',') break;
  }
  std::fclose(f);
  return n > 0;
}

std::map<MPI_Comm, std::shared_ptr<DeviceCommunicator>> &registry() {
  static std::map<MPI_Comm, std::shared_ptr<DeviceCommunicator>> communicators;
  return communicators;
}
} // namespace

void DeviceCommunicator::sum(void *const buffer, const int count,
                             const enum PrimitiveType type,
                             void *const stream) const {
  require_whole(impl_->whole);
  if (count <= 0) return;
  SMESH_NCCL_CHECK(ncclAllReduce(buffer, buffer, (size_t)count, nccl_type(type),
                                 ncclSum, ready(*impl_),
                                 static_cast<cudaStream_t>(stream)));
}

void DeviceCommunicator::sum_and_min_last(void *const buffer, const int count,
                                          const enum PrimitiveType type,
                                          void *const stream) const {
  require_whole(impl_->whole);
  if (count <= 0) return;
  char *const b = static_cast<char *>(buffer);
  const cudaStream_t s = static_cast<cudaStream_t>(stream);
  ncclComm_t comm = ready(*impl_);
  // Two reductions launched as one.
  SMESH_NCCL_CHECK(ncclGroupStart());
  if (count > 1) {
    SMESH_NCCL_CHECK(ncclAllReduce(b, b, (size_t)count - 1, nccl_type(type),
                                   ncclSum, comm, s));
  }
  char *const last = b + (size_t)(count - 1) * nccl_bytes(type);
  SMESH_NCCL_CHECK(ncclAllReduce(last, last, 1, nccl_type(type), ncclMin,
                                 comm, s));
  SMESH_NCCL_CHECK(ncclGroupEnd());
}

void DeviceCommunicator::group_begin() const {
  ready(*impl_); // outside the group
  SMESH_NCCL_CHECK(ncclGroupStart());
}
void DeviceCommunicator::group_end() const { SMESH_NCCL_CHECK(ncclGroupEnd()); }

void DeviceCommunicator::send(const void *const buffer, const std::size_t bytes,
                              const int comm_rank, void *const stream) const {
  SMESH_NCCL_CHECK(ncclSend(buffer, bytes, ncclChar,
                            peer_on_node(impl_->node_rank, comm_rank), ready(*impl_),
                            static_cast<cudaStream_t>(stream)));
}

void DeviceCommunicator::recv(void *const buffer, const std::size_t bytes,
                              const int comm_rank, void *const stream) const {
  SMESH_NCCL_CHECK(ncclRecv(buffer, bytes, ncclChar,
                            peer_on_node(impl_->node_rank, comm_rank), ready(*impl_),
                            static_cast<cudaStream_t>(stream)));
}

std::shared_ptr<DeviceCommunicator> Communicator::device() const {
  auto &communicators = registry();
  auto it = communicators.find(const_cast<Communicator *>(this)->get());
  return it == communicators.end() ? nullptr : it->second;
}

void Communicator::enable_device() const {
  MPI_Comm comm = const_cast<Communicator *>(this)->get();
  auto &communicators = registry();
  if (communicators.count(comm)) return;

  std::shared_ptr<DeviceCommunicator> ret(new DeviceCommunicator());
  auto &impl = *ret->impl_;
  int comm_rank = 0, comm_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm node;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, comm_rank, MPI_INFO_NULL, &node);
  MPI_Comm_rank(node, &impl.rank);
  MPI_Comm_size(node, &impl.size);
  impl.whole = impl.size == comm_size;

  // Every rank of the communicator, by its rank there, to its rank on this node.
  MPI_Group comm_group, node_group;
  MPI_Comm_group(comm, &comm_group);
  MPI_Comm_group(node, &node_group);
  std::vector<int> all(comm_size);
  for (int r = 0; r < comm_size; r++) all[r] = r;
  impl.node_rank.resize(comm_size);
  MPI_Group_translate_ranks(comm_group, comm_size, all.data(), node_group,
                            impl.node_rank.data());
  for (int &r : impl.node_rank) {
    if (r == MPI_UNDEFINED) r = -1;
  }
  MPI_Group_free(&comm_group);
  MPI_Group_free(&node_group);

  // The communicator spans one node, so NCCL's network plugin is never used; loading the
  // Slingshot one (aws-ofi-nccl) costs set-up time and crashed a probe at ncclCommDestroy.
  setenv("NCCL_NET_PLUGIN", "none", 0);
  ncclUniqueId id;
  if (impl.rank == 0) SMESH_NCCL_CHECK(ncclGetUniqueId(&id));
  MPI_Bcast(&id, (int)sizeof(id), MPI_BYTE, 0, node);
  MPI_Comm_free(&node);
  int device = 0;
  SMESH_CUDA_CHECK(cudaGetDevice(&device));
  DeviceCommunicator::Impl *const p = &impl;
  // The thread inherits this one's CPUs, under OpenMP's binding a single core, which the two
  // would share: the caller's set-up then took the 3 s longer (mesh input 1.0 -> 4.3 s, FF200 on 4
  // GPUs). It runs next to the GPU instead.
  cpu_set_t cpus;
  const bool near = cpus_near(device, &cpus);
  impl.init = std::thread([p, id, device, near, cpus]() {
    const double start = time_seconds();
    if (near) pthread_setaffinity_np(pthread_self(), sizeof(cpus), &cpus);
    cudaSetDevice(device);
    p->init_result = ncclCommInitRank(&p->comm, p->size, id, p->rank);
    p->init_seconds = time_seconds() - start;
  });

  communicators[comm] = ret;
}

void release_device_communicators() { registry().clear(); }

} // namespace smesh

#endif
