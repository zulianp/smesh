#ifndef SMESH_DISTRIBUTED_WRITE_HPP
#define SMESH_DISTRIBUTED_WRITE_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"
#include "smesh_types.hpp"

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

template <typename large_idx_t>
int write_mapped_field(MPI_Comm comm, const Path &output_path,
                       const ptrdiff_t n_local, const ptrdiff_t n_global,
                       const large_idx_t *const SMESH_RESTRICT mapping,
                       MPI_Datatype data_type,
                       const void *const SMESH_RESTRICT data);

// Distributed mesh topology writer for single-block meshes.
// Writes connectivity (i*.ext) and coordinates (x/y/z.ext) in the same
// on-disk format as the serial mesh_to_folder helpers, using a mapping-based
// redistribution across MPI ranks.
int write_distributed_mesh_topology(
    MPI_Comm comm, const Path &path, enum ElemType element_type,
    int spatial_dim, ptrdiff_t n_global_elements, ptrdiff_t n_owned_elements,
    const large_idx_t *element_mapping, int nnodesxelem, idx_t **local_elements,
    ptrdiff_t n_global_nodes, ptrdiff_t n_owned_nodes,
    const large_idx_t *node_mapping, geom_t **local_points);
} // namespace smesh

#endif // SMESH_DISTRIBUTED_WRITE_HPP