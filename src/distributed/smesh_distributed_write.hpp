#ifndef SMESH_DISTRIBUTED_WRITE_HPP
#define SMESH_DISTRIBUTED_WRITE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"

#include <mpi.h>

namespace smesh {

int write_mapped_field(MPI_Comm comm, const Path &output_path,
                       const ptrdiff_t n_local, const ptrdiff_t n_global,
                       const idx_t *const SMESH_RESTRICT mapping,
                       MPI_Datatype data_type,
                       const void *const SMESH_RESTRICT data);
}

#endif // SMESH_DISTRIBUTED_WRITE_HPP