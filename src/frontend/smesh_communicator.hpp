#ifndef SMESH_DISTRIBUTED_COMMUNICATOR_HPP
#define SMESH_DISTRIBUTED_COMMUNICATOR_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#ifdef SMESH_ENABLE_MPI
#include <mpi.h>
#endif

#include <functional>
#include <iostream>
#include <memory>

namespace smesh {

class Communicator {
public:
  Communicator();
  ~Communicator();
  static std::shared_ptr<Communicator> world();
  static std::shared_ptr<Communicator> null();
  static std::shared_ptr<Communicator> self();

  int rank() const;
  int size() const;
  void barrier() const;

#ifdef SMESH_ENABLE_MPI
  static std::shared_ptr<Communicator> wrap(MPI_Comm comm);
  Communicator(MPI_Comm comm);
  MPI_Comm &get();
#endif

  void
  print_callback(const std::function<void(std::ostream &)> &callback) const;

  template <typename T>
  void broadcast(T *const value, int count, int root) const {
    broadcast(value, count, TypeToEnum<T>::value(), root);
  }

  void broadcast(void *buffer, int count, enum PrimitiveType type,
                 int root) const;

  template <typename T> T sum(const T &value) const {
    T result = value;
    sum(&result, 1, TypeToEnum<T>::value());
    return result;
  }

  void sum(void *buffer, int count, enum PrimitiveType type) const;

  template <typename T> T max(const T &value) const {
    T result = value;
    max(&result, 1, TypeToEnum<T>::value());
    return result;
  }

  void max(void *buffer, int count, enum PrimitiveType type) const;

  /// Every rank contributes @p count elements from @p sendbuf; @p recvbuf receives
  /// size() * count elements, rank-major. Unlike sum/max/broadcast this is NOT in place:
  /// send and receive buffers are distinct, and at one rank -- or in a build without MPI --
  /// the contribution is copied straight across, because a gather that left the receive
  /// buffer untouched would hand the caller uninitialised memory rather than its own data.
  template <typename T>
  void allgather(const T *const sendbuf, T *const recvbuf, int count) const {
    allgather(sendbuf, recvbuf, count, TypeToEnum<T>::value());
  }

  void allgather(const void *sendbuf, void *recvbuf, int count,
                 enum PrimitiveType type) const;

  /// The variable-length form: rank r contributes @p sendcount elements and its block lands
  /// at @p displs[r] in @p recvbuf. @p recvcounts and @p displs are size() long and must
  /// agree on every rank -- the usual way to build them is an allgather of the local counts.
  template <typename T>
  void allgatherv(const T *const sendbuf, int sendcount, T *const recvbuf,
                  const int *const recvcounts, const int *const displs) const {
    allgatherv(sendbuf, sendcount, recvbuf, recvcounts, displs,
               TypeToEnum<T>::value());
  }

  void allgatherv(const void *sendbuf, int sendcount, void *recvbuf,
                  const int *recvcounts, const int *displs,
                  enum PrimitiveType type) const;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif // SMESH_DISTRIBUTED_COMMUNICATOR_HPP