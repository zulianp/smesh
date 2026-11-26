#include "smesh_context.hpp"
#include "smesh_communicator.hpp"
#include "smesh_device_communicator.hpp"
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

#if defined(__HIP_PLATFORM_SPIRV__)
#include <algorithm>
#include <mutex>
#include <vector>
#endif

namespace smesh {

#if defined(__HIP_PLATFORM_SPIRV__)
static std::mutex &teardown_mutex() {
  static std::mutex m;
  return m;
}

static std::vector<void (*)()> &teardown_releases() {
  static std::vector<void (*)()> releases;
  return releases;
}

void at_context_teardown(void (*release)()) {
  std::lock_guard<std::mutex> lock(teardown_mutex());
  auto &releases = teardown_releases();
  if (std::find(releases.begin(), releases.end(), release) == releases.end()) releases.push_back(release);
}

// The caches' device memory, released in the reverse order of registration. A release may register
// (a cache first touched there), so they run outside the lock, and those registered meanwhile run
// after them, until none is new. The list is kept: a cache registers once for the process (a
// static's initialisation), and a release only frees what its cache holds, which the next use
// allocates again, so a later owning context (one after another, without MPI) releases again.
static void release_registered_caches() {
  size_t n_run = 0;
  for (;;) {
    std::vector<void (*)()> releases;
    {
      std::lock_guard<std::mutex> lock(teardown_mutex());
      const auto &registered = teardown_releases();
      if (registered.size() == n_run) break;
      releases.assign(registered.begin() + n_run, registered.end());
      n_run = registered.size();
    }
    for (auto it = releases.rbegin(); it != releases.rend(); ++it) (*it)();
  }
}
#endif

#if defined(SMESH_ENABLE_CUDA) && defined(SMESH_ENABLE_MPI)
// The node-local rank the launcher gives the process, or -1. It is read before MPI_Init because
// Cray MPICH starts CUDA there when the job spans nodes, after which CUDA_VISIBLE_DEVICES is
// ignored and every rank of a node shares its first GPU. An MPI launcher's own variable is read
// before SLURM_LOCALID: processes that mpiexec starts inside a batch job inherit the batch step's
// SLURM_LOCALID, 0 for every rank, while the launcher sets its own per process.
// variable is set to the name of the variable read, or to "none".
static int launcher_local_rank(const char **variable) {
  for (const char *name : {"PALS_LOCAL_RANKID", "OMPI_COMM_WORLD_LOCAL_RANK", "MPI_LOCALRANKID",
                           "PMI_LOCAL_RANK", "MV2_COMM_WORLD_LOCAL_RANK", "SLURM_LOCALID"}) {
    if (const char *value = std::getenv(name)) {
      *variable = name;
      return std::atoi(value);
    }
  }
  *variable = "none";
  return -1;
}

// Each rank to one GPU of its node: its node-local rank modulo the GPUs it may use, made the only
// one it sees (CUDA_VISIBLE_DEVICES), as CSCS binds ranks to GH200 modules. Nothing touches CUDA
// here: its start-up costs seconds per rank, which a run on the host never pays.
// Returns the number of GPUs the rank could have used.
static size_t see_only_the_local_gpu(const int local_rank) {
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
  if (gpus.size() > 1 && local_rank >= 0) {
    setenv("CUDA_VISIBLE_DEVICES", gpus[local_rank % gpus.size()].c_str(), 1);
  }
  return gpus.size();
}

// The binding is only right if the launcher's node-local rank is MPI's: fails when ranks that
// share a node and its GPUs were not told apart before MPI_Init.
static void check_the_local_gpu(MPI_Comm comm, const int launcher_rank,
                                const char *launcher_variable, const size_t gpus) {
  MPI_Comm node;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
  int local_rank = 0, local_size = 1;
  MPI_Comm_rank(node, &local_rank);
  MPI_Comm_size(node, &local_size);
  MPI_Comm_free(&node);
  if (gpus > 1 && local_size > 1 && launcher_rank != local_rank) {
    SMESH_ERROR("Rank %d of %d on its node cannot be bound to one of the node's %zu GPUs: the "
                "launcher's node-local rank (%s) is %d. Set CUDA_VISIBLE_DEVICES per rank.\n",
                local_rank, local_size, gpus, launcher_variable, launcher_rank);
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
#ifdef SMESH_ENABLE_CUDA
  const char *launcher_variable = "none";
  const int launcher_rank = launcher_local_rank(&launcher_variable);
  const size_t gpus = see_only_the_local_gpu(launcher_rank);
#endif
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
  check_the_local_gpu(MPI_COMM_WORLD, launcher_rank, launcher_variable, gpus);
#endif
#endif
  SMESH_UNUSED(argc);
  SMESH_UNUSED(argv);
  impl_->owns_mpi_context = true;
  impl_->communicator = Communicator::world();

  Tracer::instance().set_rank(impl_->communicator->rank());
}

Context::~Context() {
#if defined(__HIP_PLATFORM_SPIRV__)
  // Only the context that initialized smesh: one wrapping a communicator leaves the caches to the
  // code that still runs (smesh_context.hpp).
  if (impl_->owns_mpi_context) release_registered_caches();
#endif
#ifdef SMESH_ENABLE_MPI
  release_device_communicators();
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
