#include "smesh_context.hpp"
#include "smesh_communicator.hpp"
#include "smesh_tracer.hpp"

#include "smesh_base.hpp"

#if defined(SMESH_ENABLE_CUDA) && defined(SMESH_ENABLE_MPI)
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>
#endif

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace smesh {

#if defined(SMESH_ENABLE_CUDA) && defined(SMESH_ENABLE_MPI)
// Each rank to one GPU of its node: its node-local rank modulo the GPUs it may use, made the only
// one it sees (CUDA_VISIBLE_DEVICES), as CSCS binds ranks to GH200 modules. Nothing touches CUDA
// here: its start-up costs seconds per rank, which a run on the host never pays.
static void see_only_the_local_gpu(MPI_Comm comm) {
  MPI_Comm node;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
  int local_rank = 0;
  MPI_Comm_rank(node, &local_rank);
  MPI_Comm_free(&node);

  std::vector<std::string> gpus;
  if (const char *visible = std::getenv("CUDA_VISIBLE_DEVICES")) {
    std::string list(visible);
    for (size_t begin = 0; begin <= list.size();) {
      const size_t end = std::min(list.find(',', begin), list.size());
      if (end > begin) gpus.push_back(list.substr(begin, end - begin));
      begin = end + 1;
    }
  } else {
    std::error_code error;
    const std::filesystem::path driver("/proc/driver/nvidia/gpus");
    if (std::filesystem::is_directory(driver, error)) {
      ptrdiff_t n = 0;
      for (const auto &entry : std::filesystem::directory_iterator(driver, error)) {
        SMESH_UNUSED(entry);
        n++;
      }
      for (ptrdiff_t i = 0; i < n; i++) gpus.push_back(std::to_string(i));
    }
  }
  if (gpus.size() > 1) {
    setenv("CUDA_VISIBLE_DEVICES", gpus[local_rank % gpus.size()].c_str(), 1);
  }
}
#endif

class Context::Impl {
public:
  bool owns_mpi_context{false};
  std::shared_ptr<Communicator> communicator;
  Impl() {}
};

Context::Context(int argc, char *argv[]) : impl_(std::make_unique<Impl>()) {
#ifdef SMESH_ENABLE_MPI
  int need_funneled = 0;
#ifdef _OPENMP
  need_funneled = omp_get_max_threads() > 1 ? 1 : 0;
#endif
  if (need_funneled) {
    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
      SMESH_ERROR(
          "MPI_Init_thread did not provide MPI_THREAD_FUNNELED (got %d)\n",
          provided);
    }
  } else {
    MPI_Init(&argc, &argv);
  }
  register_mpi_datatypes();
#ifdef SMESH_ENABLE_CUDA
  see_only_the_local_gpu(MPI_COMM_WORLD);
#endif
#endif
  SMESH_UNUSED(argc);
  SMESH_UNUSED(argv);
  impl_->owns_mpi_context = true;
  impl_->communicator = Communicator::world();

  Tracer::instance().set_rank(impl_->communicator->rank());
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

  Tracer::instance().set_rank(impl_->communicator->rank());
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
