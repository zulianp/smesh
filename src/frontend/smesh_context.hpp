#ifndef SMESH_CONTEXT_HPP
#define SMESH_CONTEXT_HPP

#include "smesh_base.hpp"

#include <memory>

#ifdef SMESH_ENABLE_MPI
#include <mpi.h>
#endif

namespace smesh {
class Communicator;

class Context {
public:
  Context(int argc, char *argv[]);
  ~Context();
  std::shared_ptr<Communicator> communicator();

  class Impl;
  std::unique_ptr<Impl> impl_;

#ifdef SMESH_ENABLE_MPI
  Context(int argc, char *argv[], MPI_Comm comm);
#endif
};

std::shared_ptr<Context> initialize(int argc, char *argv[]);
std::shared_ptr<Context> initialize_serial(int argc, char *argv[]);

#if defined(__HIP_PLATFORM_SPIRV__)
// chipStar (HIP on SPIR-V, Intel GPUs) faults when device memory is freed after main returns: exit
// destroys the thread_local state of the main thread, chipStar's among it, before any static's
// destructor or atexit handler runs. A static that keeps device memory for the process (a cache)
// registers its release here, once, and does not free in its destructor. The destructor of the
// context that initialized smesh (initialize, not the wrapper of a communicator) runs the releases
// while the runtime is whole, newest first, with any registered while they run. A cache used after
// that is not released: its memory goes with the process.
void at_context_teardown(void (*release)());
#endif

#ifdef SMESH_ENABLE_MPI
std::shared_ptr<Context> initialize(int argc, char *argv[], MPI_Comm comm);
#endif
} // namespace smesh

#endif // SMESH_CONTEXT_HPP
