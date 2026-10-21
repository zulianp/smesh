// DeviceCommunicator on device buffers against the values every rank can compute, many times over
// with different values: its map of the ranks on this node; messages between every two of them;
// and, where the node holds every rank, the reductions, f32 and f64 (whole numbers, so every sum
// is exact): the sums, and the sums with the last value the minimum. A build with MPI and NCCL; on
// any number of ranks and nodes.

#include "smesh_communicator.hpp"
#include "smesh_context.hpp"
#include "smesh_test.hpp"

#if defined(SMESH_ENABLE_MPI) && defined(SMESH_ENABLE_NCCL)
#include "smesh_device_communicator.hpp"

#include <cuda_runtime.h>
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <vector>

using namespace smesh;

namespace {
constexpr int count = 5;
constexpr int repeats = 200;

template <typename T> int check(const std::shared_ptr<Communicator> &comm) {
  auto device = comm->device();
  const int rank = comm->rank(), size = comm->size();

  T *d = nullptr;
  SMESH_TEST_ASSERT(cudaMalloc(&d, count * sizeof(T)) == cudaSuccess);
  long wrong = 0;
  for (int rep = 0; rep < repeats; rep++) {
    for (const bool min_last : {false, true}) {
      std::vector<T> v(count);
      for (int k = 0; k < count; k++) v[k] = T((rank + 1) * (k + 1) + rep);
      // The minimum's rank changes with the repeat.
      if (min_last) v[count - 1] = T(1000 - (rank + rep) % size);
      cudaMemcpy(d, v.data(), count * sizeof(T), cudaMemcpyHostToDevice);
      if (min_last) {
        device->sum_and_min_last(d, count, TypeToEnum<T>::value());
      } else {
        device->sum(d, count, TypeToEnum<T>::value());
      }
      cudaMemcpy(v.data(), d, count * sizeof(T), cudaMemcpyDeviceToHost);
      const int sums = min_last ? count - 1 : count;
      for (int k = 0; k < sums; k++) {
        wrong += v[k] != T((k + 1) * size * (size + 1) / 2 + size * rep);
      }
      if (min_last) wrong += v[count - 1] != T(1000 - (size - 1));
    }
  }
  cudaFree(d);
  const long all = comm->sum(wrong);
  if (!rank) {
    printf("  %d ranks, %s: %ld wrong values in %d reductions\n", size,
           sizeof(T) == 4 ? "f32" : "f64", all, 2 * repeats);
  }
  SMESH_TEST_ASSERT(all == 0);
  return SMESH_TEST_SUCCESS;
}
} // namespace

int test_node_ranks() {
  auto comm = Communicator::world();
  // None until a run enables it.
  SMESH_TEST_ASSERT(comm->device() == nullptr);
  comm->enable_device();
  comm->enable_device();
  auto device = comm->device();
  SMESH_TEST_ASSERT(device != nullptr);
  const int size = comm->size();
  MPI_Comm node;
  MPI_Comm_split_type(comm->get(), MPI_COMM_TYPE_SHARED, comm->rank(), MPI_INFO_NULL, &node);
  int node_size = 0;
  MPI_Comm_size(node, &node_size);
  MPI_Comm_free(&node);
  int on_node = 0;
  for (int r = 0; r < size; r++) on_node += device->node_rank(r) >= 0;
  SMESH_TEST_ASSERT(device->size() == node_size && on_node == node_size);
  SMESH_TEST_ASSERT(device->node_rank(comm->rank()) == device->rank());
  SMESH_TEST_ASSERT(device->whole() == (node_size == size));
  // One device communicator per MPI communicator.
  SMESH_TEST_ASSERT(device == Communicator::world()->device());
  return SMESH_TEST_SUCCESS;
}

// Each rank sends every other rank on its node `n` values of its own, all in one group.
int test_messages() {
  auto comm = Communicator::world();
  auto device = comm->device();
  const int rank = comm->rank(), size = comm->size();
  constexpr int n = 6000;
  std::vector<int> peers;
  for (int r = 0; r < size; r++) {
    if (r != rank && device->node_rank(r) >= 0) peers.push_back(r);
  }
  const int m = (int)peers.size();
  double *send = nullptr, *recv = nullptr;
  SMESH_TEST_ASSERT(cudaMalloc(&send, std::max(1, m) * n * sizeof(double)) == cudaSuccess);
  SMESH_TEST_ASSERT(cudaMalloc(&recv, std::max(1, m) * n * sizeof(double)) == cudaSuccess);
  std::vector<double> h(m * n);
  long wrong = 0;
  for (int rep = 0; rep < repeats; rep++) {
    for (int p = 0; p < m; p++) {
      for (int i = 0; i < n; i++) h[p * n + i] = rank * 1e6 + peers[p] * 1e4 + i + rep * 1e9;
    }
    cudaMemcpy(send, h.data(), m * n * sizeof(double), cudaMemcpyHostToDevice);
    device->group_begin();
    for (int p = 0; p < m; p++) {
      device->recv(recv + p * n, n * sizeof(double), peers[p]);
      device->send(send + p * n, n * sizeof(double), peers[p]);
    }
    device->group_end();
    cudaMemcpy(h.data(), recv, m * n * sizeof(double), cudaMemcpyDeviceToHost);
    for (int p = 0; p < m; p++) {
      for (int i = 0; i < n; i++) wrong += h[p * n + i] != peers[p] * 1e6 + rank * 1e4 + i + rep * 1e9;
    }
  }
  cudaFree(send);
  cudaFree(recv);
  const long all = comm->sum(wrong);
  if (!rank) {
    printf("  %d ranks, %d on this node: %ld wrong values in %d exchanges\n", size, device->size(), all, repeats);
  }
  SMESH_TEST_ASSERT(all == 0);
  return SMESH_TEST_SUCCESS;
}

int test_reductions() {
  auto comm = Communicator::world();
  if (!comm->device()->whole()) return SMESH_TEST_SUCCESS;
  SMESH_TEST_ASSERT(check<double>(comm) == SMESH_TEST_SUCCESS);
  SMESH_TEST_ASSERT(check<float>(comm) == SMESH_TEST_SUCCESS);
  return SMESH_TEST_SUCCESS;
}
#endif

int main(int argc, char *argv[]) {
  SMESH_UNIT_TEST_INIT(argc, argv);
#if defined(SMESH_ENABLE_MPI) && defined(SMESH_ENABLE_NCCL)
  SMESH_RUN_TEST(test_node_ranks);
  SMESH_RUN_TEST(test_messages);
  SMESH_RUN_TEST(test_reductions);
#else
  SMESH_UNUSED(smesh_test_argc__);
  SMESH_UNUSED(smesh_test_argv__);
#endif
  SMESH_UNIT_TEST_FINALIZE();
  return SMESH_UNIT_TEST_ERR();
}
