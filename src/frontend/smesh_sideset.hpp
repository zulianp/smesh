#ifndef SMESH_SIDESET_HPP
#define SMESH_SIDESET_HPP

#include "smesh_base.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_buffer.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <mpi.h>

namespace smesh {

class Sideset final {
public:
  int read(const std::shared_ptr<Communicator> &comm, const Path &path,
           block_idx_t block_id = 0);
  SharedBuffer<element_idx_t> parent();
  SharedBuffer<i16> lfi();
  block_idx_t block_id() const;
  static std::shared_ptr<Sideset>
  create_from_file(const std::shared_ptr<Communicator> &comm, const Path &path,
                   block_idx_t block_id = 0);
  ptrdiff_t size() const;
  std::shared_ptr<Communicator> comm() const;
  int write(const Path &path) const;

  Sideset(const std::shared_ptr<Communicator> &comm,
          const SharedBuffer<element_idx_t> &parent,
          const SharedBuffer<i16> &lfi, block_idx_t block_id = 0);
  Sideset();
  ~Sideset();

  static std::shared_ptr<Sideset>
  create(const std::shared_ptr<Communicator> &comm,
         const SharedBuffer<element_idx_t> &parent,
         const SharedBuffer<i16> &lfi, block_idx_t block_id = 0);


  static std::vector<std::shared_ptr<Sideset>> create_from_selector(
      const std::shared_ptr<Mesh> &mesh,
      const std::function<bool(const geom_t, const geom_t, const geom_t)>
          &selector,
      const std::vector<std::string> &block_names = {});

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>>
create_surface_from_sideset(
    const std::shared_ptr<Mesh> &space,
    const std::shared_ptr<Sideset> &sideset);


std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>>
create_surface_from_sidesets(
    const std::shared_ptr<Mesh> &space,
    const std::vector<std::shared_ptr<Sideset>> &sideset);

SharedBuffer<idx_t> create_nodeset_from_sideset(
    const std::shared_ptr<Mesh> &space,
    const std::shared_ptr<Sideset> &sideset);

} // namespace smesh

#endif // SMESH_SIDESET_HPP