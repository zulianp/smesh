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
  enum class ExchangeScope {
    GhostsOnly,
    GhostsAndAura,
  };

  static std::shared_ptr<Exchange>
  create_nodal(const std::shared_ptr<Mesh> &mesh,
               const ExchangeScope exchange_scope = ExchangeScope::GhostsOnly);

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
  template <typename T> int gather(T* const inout);

  class Impl;

private:
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif // SMESH_EXCHANGE_CPP
