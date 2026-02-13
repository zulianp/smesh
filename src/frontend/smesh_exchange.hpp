#ifndef SMESH_EXCHANGE_CPP
#define SMESH_EXCHANGE_CPP

#include "smesh_exchange.hpp"

#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_types.hpp"

#include <memory>

namespace smesh {

class Exchange final {
public:
  static std::shared_ptr<Exchange> create_nodal(const std::shared_ptr<Mesh> &mesh);

#if defined(SMESH_ENABLE_MPI)
  static std::shared_ptr<Exchange>
  create(const std::shared_ptr<Communicator> &comm, const ptrdiff_t nnodes,
         const ptrdiff_t n_owned_nodes, int *const node_owner,
         const idx_t *const node_offsets, const idx_t *const ghosts);
#endif
  Exchange(const std::shared_ptr<Communicator> &comm);
  ~Exchange();

  template <typename T> int exchange_add(T *const inout);

  class Impl;

private:
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif // SMESH_EXCHANGE_CPP
