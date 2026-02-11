#include "smesh_distributed_base.hpp"

namespace smesh {

int register_mpi_datatypes() {
  MPI_Type_contiguous(2, MPI_BYTE, &SMESH_MPI_F16);
  MPI_Type_commit(&SMESH_MPI_F16);
  return SMESH_SUCCESS;
}

int unregister_mpi_datatypes() {
  MPI_Type_free(&SMESH_MPI_F16);
  return SMESH_SUCCESS;
}
} // namespace smesh
