#ifndef SMESH_DISTRIBUTED_COMMUNICATOR_HPP
#define SMESH_DISTRIBUTED_COMMUNICATOR_HPP

#include "smesh_base.hpp"

#ifdef SMESH_ENABLE_MPI
#include <mpi.h>
#endif

#include <memory>

namespace smesh {

class Communicator {
public:
    Communicator();
    ~Communicator();
    static std::shared_ptr<Communicator> world();
    static std::shared_ptr<Communicator> null();
    static std::shared_ptr<Communicator> self();

    int rank() const;
    int size() const;

#ifdef SMESH_ENABLE_MPI
    static std::shared_ptr<Communicator> wrap(MPI_Comm comm);
    Communicator(MPI_Comm comm);
    MPI_Comm &get();
#endif

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace smesh

#endif // SMESH_DISTRIBUTED_COMMUNICATOR_HPP