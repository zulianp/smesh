#ifndef SMESH_DISTRIBUTED_READ_HPP
#define SMESH_DISTRIBUTED_READ_HPP

#include "smesh_base.hpp"

#include <mpi.h>

namespace smesh {

template <typename idx_t>
int read_mapped_field(MPI_Comm comm, const char *input_path,
                      const ptrdiff_t n_local, const ptrdiff_t n_global,
                      const idx_t *SMESH_RESTRICT const mapping, MPI_Datatype data_type,
                      void *SMESH_RESTRICT const data_out);
}

#endif // SMESH_DISTRIBUTED_READ_HPP
