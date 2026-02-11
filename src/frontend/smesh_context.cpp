#include "smesh_context.hpp"
#include "smesh_communicator.hpp"

#include "smesh_base.hpp"

#ifdef SMESH_ENABLE_CUDA
// TODO: include CUDA init stuff
#endif

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#endif

namespace smesh {

class Context::Impl {
public:
  bool owns_mpi_context{false};
  std::shared_ptr<Communicator> communicator;
  Impl() {}
};

Context::Context(int argc, char *argv[]) : impl_(std::make_unique<Impl>()) {
#ifdef SMESH_ENABLE_MPI

  MPI_Init(&argc, &argv);
  register_mpi_datatypes();
#endif
  impl_->owns_mpi_context = true;
  impl_->communicator = Communicator::world();
}

Context::~Context() {
#ifdef SMESH_ENABLE_MPI
  unregister_mpi_datatypes();
  if (impl_->owns_mpi_context) {
    MPI_Finalize();
  }
#endif
}

std::shared_ptr<Communicator> Context::communicator() {
  return impl_->communicator;
}

std::shared_ptr<Context> initialize(int argc, char *argv[]) {
  return std::make_shared<Context>(argc, argv);
}

#ifdef SMESH_ENABLE_MPI
Context::Context(int argc, char *argv[], MPI_Comm comm)
    : impl_(std::make_unique<Impl>()) {
  SMESH_UNUSED(argc);
  SMESH_UNUSED(argv);
  impl_->communicator = Communicator::wrap(comm);
}

std::shared_ptr<Context> initialize(int argc, char *argv[], MPI_Comm comm) {
  return std::make_shared<Context>(argc, argv, comm);
}

std::shared_ptr<Context> initialize_serial(int argc, char *argv[]) {

  auto context = std::make_shared<Context>(argc, argv);

  if (context->communicator()->size() > 1) {
    SMESH_ERROR("Context initialization does not support parallel runs!\n");
  }

  return context;
}
#else
std::shared_ptr<Context> initialize_serial(int argc, char *argv[]) {
  return std::make_shared<Context>(argc, argv);
}

#endif

} // namespace smesh
