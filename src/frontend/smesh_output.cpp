#include "smesh_output.hpp"
#include "smesh_env.hpp"
#include "smesh_common.hpp"
#include "smesh_glob.hpp"

#if defined(SMESH_ENABLE_MPI)
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#endif

namespace smesh {

class Output::Impl {
public:
  std::shared_ptr<Mesh> mesh;
  Path path;
  Impl(const std::shared_ptr<Mesh> &mesh, const Path &path)
      : mesh(mesh), path(path) {}
  ~Impl() {}
};

std::shared_ptr<Mesh> Output::mesh() const { return impl_->mesh; }

Output::Output(const std::shared_ptr<Mesh> &mesh, const Path &path)
    : impl_(std::make_unique<Impl>(mesh, path)) {}

Output::~Output() = default;

std::shared_ptr<Output> Output::create(const std::shared_ptr<Mesh> &mesh,
                                       const Path &path) {
  return std::make_shared<Output>(mesh, path);
}

int Output::write_nodal(const std::string &field_name,
                        const enum PrimitiveType type, const void *const data,
                        const int block_size) {
  if (block_size != 1) {
    SMESH_ERROR("write_nodal: Block size must be 1 for nodal data");
    return SMESH_FAILURE;
  }

  auto mesh = this->mesh();
  if(mesh->comm()->rank() == 0) {
    create_directory(impl_->path);
  }

  auto path = Path(impl_->path) / (field_name + "." + str(type));
  
#if defined(SMESH_ENABLE_MPI)
  if (mesh->comm()->size() > 1) {
    auto dist = mesh->distributed();
    auto data_type = mpi_type_from_primitive_type(type);
    return write_mapped_field(
      mesh->comm()->get(), path, dist->n_nodes_owned(), dist->n_nodes_global(),
      dist->node_mapping()->data(), data_type, data);
  }
#endif
  return array_write(path, type, data, mesh->n_nodes());
}

int Output::write_elemental(const std::string &field_name,
                            const enum PrimitiveType type,
                            const void *const data, const int block_size) {
  if (block_size != 1) {
    SMESH_ERROR("write_elemental: Block size must be 1 for elemental data");
    return SMESH_FAILURE;
  }

  auto mesh = this->mesh();
  if(mesh->comm()->rank() == 0) {
    create_directory(impl_->path);
  }

  auto path = Path(impl_->path) / (field_name + "." + str(type));

#if defined(SMESH_ENABLE_MPI)
  if (mesh->comm()->size() > 1) {
    auto dist = mesh->distributed();
    auto data_type = mpi_type_from_primitive_type(type);
    return write_mapped_field(mesh->comm()->get(), path, dist->n_elements_owned(), dist->n_elements_global(),
                              dist->element_mapping()->data(), data_type, data);
  }
#endif
  return array_write(path, type, data, mesh->n_elements());
}

} // namespace smesh
