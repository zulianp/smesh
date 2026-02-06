#ifndef SMESH_DISTRIBUTED_READ_HPP
#define SMESH_DISTRIBUTED_READ_HPP

#include "smesh_base.hpp"

#include "smesh_path.hpp"

#include <mpi.h>

namespace smesh {

template <typename idx_t>
int read_mapped_field(MPI_Comm comm, const char *input_path,
                      const ptrdiff_t n_local, const ptrdiff_t n_global,
                      const idx_t *SMESH_RESTRICT const mapping,
                      MPI_Datatype data_type,
                      void *SMESH_RESTRICT const data_out);

template <typename geom_t>
int mesh_coordinates_from_folder(MPI_Comm comm, const Path &folder,
                                 int *spatial_dim_out, geom_t ***points_out,
                                 ptrdiff_t *n_local_nodes_out,
                                 ptrdiff_t *n_global_nodes_out);

template <typename idx_t>
int mesh_block_from_folder(MPI_Comm comm, const Path &folder,
                           int *nnodesxelem_out, idx_t ***const elems,
                           ptrdiff_t *const n_local_elements_out,
                           ptrdiff_t *const n_global_elements_out);

template <typename idx_t, typename geom_t, typename element_idx_t>
int mesh_from_folder(
    const MPI_Comm comm, const Path &folder, int *nnodesxelem_out,
    ptrdiff_t *nelements_out, idx_t ***elements_out, int *spatial_dim_out,
    ptrdiff_t *nnodes_out, geom_t ***points_out, ptrdiff_t *n_owned_nodes_out,
    ptrdiff_t *n_owned_elements_out, element_idx_t **element_mapping_out,
    idx_t **node_mapping_out, int **node_owner_out, idx_t **node_offsets_out,
    idx_t **ghosts_out, ptrdiff_t *n_owned_nodes_with_ghosts_out,
    ptrdiff_t *n_shared_elements_out,
    ptrdiff_t *n_owned_elements_with_ghosts_out);
} // namespace smesh

#endif // SMESH_DISTRIBUTED_READ_HPP
