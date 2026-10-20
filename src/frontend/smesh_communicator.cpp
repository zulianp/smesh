#include "smesh_communicator.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#endif

#include <algorithm>
#include <cstring>
#include <map>
#include <sstream>
#include <utility>

namespace smesh {

class Communicator::Impl {
public:
  Impl();
  ~Impl();

#ifdef SMESH_ENABLE_MPI
  MPI_Comm comm;
  Impl(MPI_Comm comm) : comm(comm) {}
#endif
};

Communicator::Communicator() : impl_(std::make_unique<Impl>()) {}
Communicator::~Communicator() {}

std::shared_ptr<Communicator> Communicator::world() {
#ifdef SMESH_ENABLE_MPI
  return std::make_shared<Communicator>(MPI_COMM_WORLD);
#else
  return std::make_shared<Communicator>();
#endif
}

std::shared_ptr<Communicator> Communicator::null() {
#ifdef SMESH_ENABLE_MPI
  return std::make_shared<Communicator>(MPI_COMM_NULL);
#else
  return std::make_shared<Communicator>();
#endif
}

std::shared_ptr<Communicator> Communicator::self() {
#ifdef SMESH_ENABLE_MPI
  return std::make_shared<Communicator>(MPI_COMM_SELF);
#else
  return std::make_shared<Communicator>();
#endif
}

#ifdef SMESH_ENABLE_MPI
Communicator::Communicator(MPI_Comm comm)
    : impl_(std::make_unique<Impl>(comm)) {}
MPI_Comm &Communicator::get() { return impl_->comm; }
std::shared_ptr<Communicator> Communicator::wrap(MPI_Comm comm) {
  return std::make_shared<Communicator>(comm);
}
#endif

int Communicator::rank() const {
#ifdef SMESH_ENABLE_MPI
  int rank;
  MPI_Comm_rank(impl_->comm, &rank);
  return rank;
#else
  return 0;
#endif
}

int Communicator::size() const {
#ifdef SMESH_ENABLE_MPI
  int size;
  MPI_Comm_size(impl_->comm, &size);
  return size;
#else
  return 1;
#endif
}

void Communicator::barrier() const {
#ifdef SMESH_ENABLE_MPI
  MPI_Barrier(impl_->comm);
#endif
}

Communicator::Impl::Impl() = default;
Communicator::Impl::~Impl() = default;

void Communicator::print_callback(
    const std::function<void(std::ostream &)> &callback) const {
  std::cout << std::flush;
  std::stringstream ss;
  callback(ss);
  int rank = this->rank();
  int size = this->size();

  barrier();

  for (int r = 0; r < size; r++) {
    if (r == rank) {
      std::cout << "----------------\n";
      std::cout << "[" << rank << "]\n" << ss.str() << "\n";
      std::cout << "----------------\n";
      std::cout << std::flush;
    }
    barrier();
  }
}

void Communicator::broadcast(void *buffer, int count, enum PrimitiveType type,
                             int root) const {
#ifdef SMESH_ENABLE_MPI
  MPI_Bcast(buffer, count, mpi_type_from_primitive_type(type), root, impl_->comm);
#else
  SMESH_UNUSED(buffer);
  SMESH_UNUSED(count);
  SMESH_UNUSED(type);
  SMESH_UNUSED(root);
#endif
}

void Communicator::sum(void *buffer, int count, enum PrimitiveType type) const {
#ifdef SMESH_ENABLE_MPI
  MPI_Allreduce(MPI_IN_PLACE, buffer, count, mpi_type_from_primitive_type(type), MPI_SUM, impl_->comm);
#else
  SMESH_UNUSED(buffer);
  SMESH_UNUSED(count);
  SMESH_UNUSED(type);
#endif
}

void Communicator::max(void *buffer, int count, enum PrimitiveType type) const {
#ifdef SMESH_ENABLE_MPI
  MPI_Allreduce(MPI_IN_PLACE, buffer, count, mpi_type_from_primitive_type(type), MPI_MAX, impl_->comm);
#else
  SMESH_UNUSED(buffer);
  SMESH_UNUSED(count);
  SMESH_UNUSED(type);
#endif
}

#ifdef SMESH_ENABLE_MPI
namespace {
// One element is a whole tuple of n values, the first n - 1 summed and the last minimised: MPI
// applies an operation to whole elements, never to part of one.
template <typename T>
void sum_and_min_last_op(void *in, void *inout, int *len, MPI_Datatype *type) {
  int bytes = 0;
  MPI_Type_size(*type, &bytes);
  const int n = bytes / (int)sizeof(T);
  const T *const a = static_cast<const T *>(in);
  T *const b = static_cast<T *>(inout);
  for (int e = 0; e < *len; e++) {
    for (int i = 0; i < n - 1; i++) b[e * n + i] += a[e * n + i];
    b[e * n + n - 1] = std::min(b[e * n + n - 1], a[e * n + n - 1]);
  }
}

// The operation and the tuple type, made once each and kept for the process.
void sum_and_min_last_handles(const enum PrimitiveType type, const int count,
                              MPI_Op *op, MPI_Datatype *tuple) {
  static std::map<std::pair<int, int>, MPI_Datatype> tuples;
  static MPI_Op ops[2] = {MPI_OP_NULL, MPI_OP_NULL};
  const bool is_f64 = type == SMESH_FLOAT64;
  if (!is_f64 && type != SMESH_FLOAT32) {
    SMESH_ERROR("Communicator::sum_and_min_last: f32 or f64\n");
  }
  MPI_Op &o = ops[is_f64 ? 1 : 0];
  if (o == MPI_OP_NULL) {
    MPI_Op_create(is_f64 ? &sum_and_min_last_op<f64> : &sum_and_min_last_op<f32>, 1, &o);
  }
  auto &t = tuples[{(int)type, count}];
  if (!t) {
    MPI_Type_contiguous(count, mpi_type_from_primitive_type(type), &t);
    MPI_Type_commit(&t);
  }
  *op = o;
  *tuple = t;
}
} // namespace
#endif

void Communicator::sum_and_min_last(void *buffer, int count, enum PrimitiveType type) const {
#ifdef SMESH_ENABLE_MPI
  MPI_Op op;
  MPI_Datatype tuple;
  sum_and_min_last_handles(type, count, &op, &tuple);
  MPI_Allreduce(MPI_IN_PLACE, buffer, 1, tuple, op, impl_->comm);
#else
  SMESH_UNUSED(buffer);
  SMESH_UNUSED(count);
  SMESH_UNUSED(type);
#endif
}

void Communicator::allgather(const void *sendbuf, void *recvbuf, int count,
                             enum PrimitiveType type) const {
#ifdef SMESH_ENABLE_MPI
  MPI_Allgather(const_cast<void *>(sendbuf), count,
                mpi_type_from_primitive_type(type), recvbuf, count,
                mpi_type_from_primitive_type(type), impl_->comm);
#else
  // Not a no-op, unlike the reductions above: with one rank the gathered result IS this
  // rank's contribution, and leaving recvbuf alone would return uninitialised memory.
  std::memcpy(recvbuf, sendbuf, (size_t)count * num_bytes(type));
#endif
}

void Communicator::allgatherv(const void *sendbuf, int sendcount, void *recvbuf,
                              const int *recvcounts, const int *displs,
                              enum PrimitiveType type) const {
#ifdef SMESH_ENABLE_MPI
  MPI_Allgatherv(const_cast<void *>(sendbuf), sendcount,
                 mpi_type_from_primitive_type(type), recvbuf,
                 const_cast<int *>(recvcounts), const_cast<int *>(displs),
                 mpi_type_from_primitive_type(type), impl_->comm);
#else
  // One rank: its block sits at displs[0], which a caller is free to make non-zero.
  const size_t bytes = num_bytes(type);
  std::memcpy(static_cast<char *>(recvbuf) + (size_t)displs[0] * bytes, sendbuf,
              (size_t)sendcount * bytes);
  SMESH_UNUSED(recvcounts);
#endif
}

} // namespace smesh
