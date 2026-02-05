#ifndef SMESH_READ_HPP
#define SMESH_READ_HPP

#include "smesh_base.hpp"
#include <string_view>

namespace smesh {

template <typename T>
int array_read(const Path &path, T **data, ptrdiff_t *n_elements);

template <typename FileType, typename TargetType>
int array_read_convert(const Path &path, TargetType **data,
                       ptrdiff_t *n_elements);

template <typename T>
int array_read_convert_from_extension(const Path &path, T **data,
                                ptrdiff_t *n_elements);

template <typename idx_t, typename geom_t>
int mesh_from_folder(const Path &folder, int *nnodesxelem_out,
                     ptrdiff_t *nelements_out, idx_t ***elems_out,
                     int *spatial_dim_out, ptrdiff_t *nnodes_out,
                     geom_t ***points_out);

} // namespace smesh

#endif // SMESH_READ_HPP