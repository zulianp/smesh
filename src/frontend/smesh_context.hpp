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

#ifdef SMESH_ENABLE_MPI
std::shared_ptr<Context> initialize(int argc, char *argv[], MPI_Comm comm);
#endif
} // namespace smesh

#endif // SMESH_CONTEXT_HPP
