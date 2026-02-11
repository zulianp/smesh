#include "smesh_distributed_base.hpp"

namespace smesh {

MPI_Datatype SMESH_MPI_F16 = MPI_DATATYPE_NULL;

int register_mpi_datatypes() {
  // f16 is treated as an opaque 16-bit payload for MPI transfers.
  MPI_Type_contiguous(static_cast<int>(sizeof(f16)), MPI_BYTE, &SMESH_MPI_F16);
  MPI_Type_commit(&SMESH_MPI_F16);
  return SMESH_SUCCESS;
}

int unregister_mpi_datatypes() {
  MPI_Type_free(&SMESH_MPI_F16);
  return SMESH_SUCCESS;
}
} // namespace smesh
