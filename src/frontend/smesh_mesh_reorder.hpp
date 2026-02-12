#ifndef SMESH_MESH_REORDER_HPP
#define SMESH_MESH_REORDER_HPP

#include "smesh_base.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_types.hpp"

#include <memory>
namespace smesh {

class SFC {
public:
  SFC();
  ~SFC();

  int reorder(Mesh &mesh);

  static std::shared_ptr<SFC> create_from_env();

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif // SMESH_MESH_REORDER_HPP
