#ifndef SMESH_DEVICE_COMMUNICATOR_HPP
#define SMESH_DEVICE_COMMUNICATOR_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <cstddef>
#include <memory>

namespace smesh {

class Communicator;

/// Communication of device buffers among the ranks of a Communicator that share this rank's node,
/// queued on a CUDA stream (null: the legacy default stream) like a kernel: a call returns at once
/// and its result is there for the work queued after it, without a pass through the host. NCCL,
/// in a build with MPI and NCCL; one per MPI communicator, for a run that enables it
/// (Communicator::enable_device(), which starts NCCL's set-up and returns; the first use waits).
///
/// Opt-in: on Alps NCCL's set-up costs 3.1 to 3.7 s once per process, most of it in GPU driver
/// calls of ~25 ms each, against 11 to 39 us saved per QP iteration on one node (FF200, 4 GPUs:
/// 142 -> 106 us); it pays from about 92,000 QP iterations a process, and the LDP cases take
/// 5,000 to 30,000 (measured 2026-10-21, wip/ldp_containment).
///
/// Only within a node: on Alps (GH200, Slingshot, NCCL 2.22 with aws-ofi-nccl 1.9) NCCL's
/// point-to-point messages between nodes delivered wrong values (4 ranks on 2 nodes, 16 KB per
/// peer: 544 to 1,376 in 50 exchanges, in a probe that had first sent the same device buffers
/// with GPU-aware MPI; measured 2026-10-21) and took twice as long as GPU-aware MPI there. Within
/// a node its exchange took a quarter of MPI's time and every value arrived.
class DeviceCommunicator {
public:
  ~DeviceCommunicator();

  /// This rank among the node's ranks, and their number.
  int rank() const;
  int size() const;

  /// Whether the node holds every rank of the communicator: then the reductions below are the
  /// communicator's.
  bool whole() const;

  /// The node rank of the communicator's rank `comm_rank`, -1 if it is on another node.
  int node_rank(const int comm_rank) const;

  /// In place, the sum over the ranks. f32 and f64; whole() communicators only.
  void sum(void *const buffer, const int count, const enum PrimitiveType type,
           void *const stream = nullptr) const;

  /// In place, the first count - 1 values summed and the last the minimum over the ranks
  /// (Communicator::sum_and_min_last). f32 and f64; whole() communicators only.
  void sum_and_min_last(void *const buffer, const int count,
                        const enum PrimitiveType type,
                        void *const stream = nullptr) const;

  /// Messages to and from the communicator's rank `comm_rank`, on this node, posted together
  /// between group_begin() and group_end(): a send matches the peer's receive of the same size.
  void group_begin() const;
  void group_end() const;
  void send(const void *const buffer, const std::size_t bytes, const int comm_rank,
            void *const stream = nullptr) const;
  void recv(void *const buffer, const std::size_t bytes, const int comm_rank,
            void *const stream = nullptr) const;

  class Impl;

private:
  friend class Communicator;
  DeviceCommunicator();
  std::unique_ptr<Impl> impl_;
};

/// Releases the device communicators Communicator::enable_device() created; Context's destructor
/// calls it before MPI_Finalize. A no-op without MPI and NCCL.
void release_device_communicators();

} // namespace smesh

#endif // SMESH_DEVICE_COMMUNICATOR_HPP
