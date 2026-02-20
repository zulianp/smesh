#include "smesh_communicator.hpp"

#include <sstream>

namespace smesh {

class Communicator::Impl {
public:
  Impl();
  ~Impl();

#ifdef SMESH_ENABLE_MPI
  MPI_Comm comm;
  Impl(MPI_Comm comm) : comm(comm) {}
#endif
};

Communicator::Communicator() : impl_(std::make_unique<Impl>()) {}
Communicator::~Communicator() {}

std::shared_ptr<Communicator> Communicator::world() {
#ifdef SMESH_ENABLE_MPI
  return std::make_shared<Communicator>(MPI_COMM_WORLD);
#else
  return std::make_shared<Communicator>();
#endif
}

std::shared_ptr<Communicator> Communicator::null() {
#ifdef SMESH_ENABLE_MPI
  return std::make_shared<Communicator>(MPI_COMM_NULL);
#else
  return std::make_shared<Communicator>();
#endif
}

std::shared_ptr<Communicator> Communicator::self() {
#ifdef SMESH_ENABLE_MPI
  return std::make_shared<Communicator>(MPI_COMM_SELF);
#else
  return std::make_shared<Communicator>();
#endif
}

#ifdef SMESH_ENABLE_MPI
Communicator::Communicator(MPI_Comm comm)
    : impl_(std::make_unique<Impl>(comm)) {}
MPI_Comm &Communicator::get() { return impl_->comm; }
std::shared_ptr<Communicator> Communicator::wrap(MPI_Comm comm) {
  return std::make_shared<Communicator>(comm);
}
#endif

int Communicator::rank() const {
#ifdef SMESH_ENABLE_MPI
  int rank;
  MPI_Comm_rank(impl_->comm, &rank);
  return rank;
#else
  return 0;
#endif
}

int Communicator::size() const {
#ifdef SMESH_ENABLE_MPI
  int size;
  MPI_Comm_size(impl_->comm, &size);
  return size;
#else
  return 1;
#endif
}

void Communicator::barrier() const {
#ifdef SMESH_ENABLE_MPI
  MPI_Barrier(impl_->comm);
#endif
}

Communicator::Impl::Impl() = default;
Communicator::Impl::~Impl() = default;

void Communicator::print_callback(
    const std::function<void(std::ostream &)> &callback) const {
  std::cout << std::flush;
  std::stringstream ss;
  callback(ss);
  int rank = this->rank();
  int size = this->size();

  barrier();

  for (int r = 0; r < size; r++) {
    if (r == rank) {
      std::cout << "----------------\n";
      std::cout << "[" << rank << "]\n" << ss.str() << "\n";
      std::cout << "----------------\n";
      std::cout << std::flush;
    }
    barrier();
  }
}

} // namespace smesh
