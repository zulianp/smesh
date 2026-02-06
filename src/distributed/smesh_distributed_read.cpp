#include "smesh_distributed_read.impl.hpp"


namespace smesh {
// Explicit instantiation
template int read_mapped_field<int32_t>(MPI_Comm comm, const char *input_path,
                                        const ptrdiff_t n_local, const ptrdiff_t n_global,
                                        const int32_t *SMESH_RESTRICT const mapping, MPI_Datatype data_type,
                                        void *const SMESH_RESTRICT data_out);

template int read_mapped_field<int64_t>(MPI_Comm comm, const char *input_path,
                                        const ptrdiff_t n_local, const ptrdiff_t n_global,
                                        const int64_t *const SMESH_RESTRICT mapping, MPI_Datatype data_type,
                                        void *const SMESH_RESTRICT data_out);

}