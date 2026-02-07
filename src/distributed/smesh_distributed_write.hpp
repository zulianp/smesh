#ifndef SMESH_DISTRIBUTED_WRITE_HPP
#define SMESH_DISTRIBUTED_WRITE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"

#include <mpi.h>

namespace smesh {

template <typename FileType, typename T>
int array_write_convert(MPI_Comm comm, const Path &path,
                        const T *const SMESH_RESTRICT data,
                        const ptrdiff_t n_local_elements,
                        const ptrdiff_t n_global_elements);

template <typename T>
int array_write_convert_from_extension(MPI_Comm comm, const Path &path,
                                       const T *const SMESH_RESTRICT data,
                                       const ptrdiff_t n_local_elements,
                                       const ptrdiff_t n_global_elements);

template <typename idx_t>
int write_mapped_field(MPI_Comm comm, const Path &output_path,
                       const ptrdiff_t n_local, const ptrdiff_t n_global,
                       const idx_t *const SMESH_RESTRICT mapping,
                       MPI_Datatype data_type,
                       const void *const SMESH_RESTRICT data);
} // namespace smesh

#endif // SMESH_DISTRIBUTED_WRITE_HPP