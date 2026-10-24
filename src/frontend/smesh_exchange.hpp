#ifndef SMESH_EXCHANGE_CPP
#define SMESH_EXCHANGE_CPP

#include "smesh_exchange.hpp"

#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_spaces.hpp"
#include "smesh_types.hpp"

#include <memory>

namespace smesh {

class Exchange final {
public:
  enum class ExchangeScope {
    GhostsOnly,
    GhostsAndAura,
  };

  /// The nodal exchange of a distributed mesh. With EXECUTION_SPACE_DEVICE its gathers take
  /// device buffers: packed and unpacked on the device, sent by GPU-aware MPI, or, where the
  /// communicator's device communication is enabled (Communicator::enable_device), to ranks on
  /// this node by NCCL, queued on the stream.
  static std::shared_ptr<Exchange>
  create_nodal(const std::shared_ptr<Mesh> &mesh,
               const ExchangeScope exchange_scope = ExchangeScope::GhostsOnly,
               const ExecutionSpace execution_space = EXECUTION_SPACE_HOST);

#if defined(SMESH_ENABLE_MPI)
  static std::shared_ptr<Exchange>
  create(const std::shared_ptr<Communicator> &comm,
         const ExchangeScope exchange_scope,
         const ptrdiff_t n_local_nodes, const ptrdiff_t n_owned_nodes,
         const int *const node_owner,
         const ptrdiff_t *const node_offsets, const idx_t *const ghosts);
#endif
  Exchange(const std::shared_ptr<Communicator> &comm);
  ~Exchange();

  template <typename T> int scatter_add(T *const inout);
  template <typename T>
  int scatter_add(T *const inout, const ptrdiff_t block_size);
  template <typename T> int gather(T *const inout);
  template <typename T> int gather(T *const inout, const ptrdiff_t block_size);

  /// Pack and post non-blocking gather (or complete synchronously if counts do not fit Ialltoallv).
  /// On the device the pack, NCCL's messages and the unpack at gather_wait are queued on `stream`
  /// (a cudaStream_t; null: the legacy default stream), and only that stream is waited for before
  /// MPI reads the packed buffer.
  template <typename T>
  int gather_begin(T *const inout, const ptrdiff_t block_size = 1, void *const stream = nullptr);
  /// Wait for a pending gather_begin and unpack into ghost/aura slots, on gather_begin's stream.
  int gather_wait();

  class Impl;

private:
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif // SMESH_EXCHANGE_CPP
