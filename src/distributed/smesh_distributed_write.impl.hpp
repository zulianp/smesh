#define SMESH_DISTRIBUTED_WRITE_IMPL_HPP

#include "matrixio_array.h"
#include "smesh_alltoallv.impl.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#include "smesh_types.hpp"

#include <chrono>
#include <fstream>
#include <limits>

namespace smesh {

template <typename FileType, typename T>
int array_write_convert(MPI_Comm comm, const Path &path,
                        const T *const SMESH_RESTRICT data,
                        const ptrdiff_t n_local_elements,
                        const ptrdiff_t n_global_elements) {
  if (std::is_same_v<FileType, T>) {
    return array_write(comm, path.c_str(), smesh::mpi_type<T>(), data,
                       n_local_elements, n_global_elements);
  }

  FileType *buffer = (FileType *)malloc(n_local_elements * sizeof(FileType));
  for (ptrdiff_t i = 0; i < n_local_elements; ++i) {
    buffer[i] = static_cast<FileType>(data[i]);
  }
  int ret = array_write(comm, path.c_str(), smesh::mpi_type<FileType>(), buffer,
                        n_local_elements, n_global_elements);

  free(buffer);
  return ret;
}

template <typename T>
int array_write_convert_from_extension(MPI_Comm comm, const Path &path,
                                       const T *const SMESH_RESTRICT data,
                                       const ptrdiff_t n_local_elements,
                                       const ptrdiff_t n_global_elements) {
  auto ext = path.extension();
  if (ext == "raw") {
    return array_write(comm, path.c_str(), smesh::mpi_type<T>(), data,
                       n_local_elements, n_global_elements);
  } else if (ext == "float16") {
    return array_write_convert<f16, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "float32") {
    return array_write_convert<f32, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "float64") {
    return array_write_convert<f64, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "int16") {
    return array_write_convert<i16, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "int32") {
    return array_write_convert<i32, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  } else if (ext == "int64") {
    return array_write_convert<i64, T>(comm, path, data, n_local_elements,
                                       n_global_elements);
  }
  // else  if (ext == "txt") {
  //   return SMESH_SUCCESS;
  // }
  else {
    SMESH_ERROR("Unsupported file extension %s for file %s\n", ext.c_str(),
                path.c_str());
    return SMESH_FAILURE;
  }
}

template <typename large_idx_t>
int write_mapped_field(MPI_Comm comm, const Path &output_path,
                       const ptrdiff_t n_local, const ptrdiff_t n_global,
                       const large_idx_t *const mapping, MPI_Datatype data_type,
                       const void *const data_in) {
  using byte_t = uint8_t;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  const byte_t *const data = (const byte_t *const)data_in;

  int type_size;
  SMESH_MPI_CATCH(MPI_Type_size(data_type, &type_size));

  const ptrdiff_t local_output_size_no_remainder = n_global / size;
  const ptrdiff_t begin = (n_global / size) * rank;

  ptrdiff_t local_output_size = local_output_size_no_remainder;
  if (rank == size - 1) {
    local_output_size = n_global - begin;
  }

  i64 *send_count = (i64 *)malloc((size) * sizeof(i64));
  memset(send_count, 0, (size) * sizeof(i64));

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const large_idx_t idx = mapping[i];
    int dest_rank = std::min(size - 1, (int)(idx / local_output_size_no_remainder));
    send_count[dest_rank]++;
  }

  i64 *recv_count = (i64 *)malloc((size) * sizeof(i64));
  SMESH_MPI_CATCH(MPI_Alltoall(send_count, 1, mpi_type<i64>(), recv_count, 1,
                               mpi_type<i64>(), comm));

  i64 *send_displs = (i64 *)malloc(size * sizeof(i64));
  i64 *recv_displs = (i64 *)malloc(size * sizeof(i64));
  i64 *book_keeping = (i64 *)calloc(size, sizeof(i64));

  send_displs[0] = 0;
  recv_displs[0] = 0;

  // Create data displacements for sending
  for (int i = 0; i < size - 1; ++i) {
    send_displs[i + 1] = send_displs[i] + send_count[i];
  }

  // Create data displacements for receiving
  for (int i = 0; i < size - 1; ++i) {
    recv_displs[i + 1] = recv_displs[i] + recv_count[i];
  }

  large_idx_t *send_list = (large_idx_t *)malloc(n_local * sizeof(large_idx_t));

  ptrdiff_t n_buff = std::max(n_local, local_output_size);
  uint8_t *send_data_and_final_storage = (uint8_t *)malloc(n_buff * type_size);

  // Pack data and indices
  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const large_idx_t idx = mapping[i];
    int dest_rank = std::min(size - 1, (int)(idx / local_output_size_no_remainder));
    SMESH_ASSERT(dest_rank < size);

    // Put index and data into buffers
    const ptrdiff_t offset = send_displs[dest_rank] + book_keeping[dest_rank];
    send_list[offset] = idx;
    memcpy((void *)&send_data_and_final_storage[offset * type_size],
           (void *)&data[i * type_size], type_size);

    book_keeping[dest_rank]++;
  }

  large_idx_t *recv_list = (large_idx_t *)malloc(local_output_size * sizeof(large_idx_t));
  uint8_t *recv_data = (uint8_t *)malloc(local_output_size * type_size);

  ///////////////////////////////////
  // Send indices
  ///////////////////////////////////

  ptrdiff_t max_chunk_size = std::numeric_limits<ptrdiff_t>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64_b(send_list, send_count, send_displs,
                                mpi_type<large_idx_t>(), recv_list, recv_count,
                                recv_displs, mpi_type<large_idx_t>(), comm, max_chunk_size));

  ///////////////////////////////////
  // Send data
  ///////////////////////////////////

  SMESH_MPI_CATCH(all_to_allv_64_b(send_data_and_final_storage, send_count,
                                send_displs, data_type, recv_data, recv_count,
                                recv_displs, data_type, comm, max_chunk_size));

  ///////////////////////////////////
  // Unpack indexed data
  ///////////////////////////////////

  for (ptrdiff_t i = 0; i < local_output_size; ++i) {
    ptrdiff_t dest = recv_list[i] - begin;
    SMESH_ASSERT(dest >= 0);
    SMESH_ASSERT(dest < local_output_size);
    memcpy((void *)&send_data_and_final_storage[dest * type_size],
           (void *)&recv_data[i * type_size], type_size);
  }

  array_write(comm, output_path.c_str(), data_type, (void *)send_data_and_final_storage,
              local_output_size, n_global);

  ///////////////////////////////////
  // Clean-up
  ///////////////////////////////////
  free(send_count);
  free(send_displs);
  free(recv_count);
  free(recv_displs);
  free(book_keeping);
  free(send_list);
  free(recv_list);
  free(recv_data);
  free(send_data_and_final_storage);
  return 0;
}

int write_distributed_mesh_topology(
    MPI_Comm comm, const Path &path, enum ElemType /*element_type*/,
    int spatial_dim, ptrdiff_t n_global_elements, ptrdiff_t n_owned_elements,
    const large_idx_t *element_mapping, int nnodesxelem, idx_t **local_elements,
    ptrdiff_t n_global_nodes, ptrdiff_t n_owned_nodes,
    const large_idx_t *node_mapping, geom_t **local_points) {
  SMESH_TRACE_SCOPE("write_distributed_mesh_topology");

  int err = SMESH_SUCCESS;

  // Write coordinates (x/y/z.*) using ownership-based mapping.
  static constexpr char xyz[3] = {'x', 'y', 'z'};
  for (int d = 0; d < spatial_dim; ++d) {
    std::string fname =
        std::string(1, xyz[d]) + "." + std::string(TypeToString<geom_t>::value());
    Path coord_path = path / fname;

    // Only owned nodes participate in mapped write.
    err |= write_mapped_field(comm, coord_path, n_owned_nodes, n_global_nodes,
                              node_mapping, smesh::mpi_type<geom_t>(),
                              local_points[d]);
  }

  if (err != SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  // Write connectivity i*.*
  for (int v = 0; v < nnodesxelem; ++v) {
    std::string fname =
        "i" + std::to_string(v) + "." + std::string(TypeToString<idx_t>::value());
    Path conn_path = path / fname;

    idx_t *buffer =
        (idx_t *)malloc((size_t)n_owned_elements * sizeof(idx_t));
    if (!buffer) {
      return SMESH_FAILURE;
    }

    idx_t *local_col = local_elements[v];
    for (ptrdiff_t e = 0; e < n_owned_elements; ++e) {
      const idx_t local_node = local_col[e];
      buffer[e] = static_cast<idx_t>(node_mapping[local_node]);
    }

    err |= write_mapped_field(comm, conn_path, n_owned_elements,
                              n_global_elements, element_mapping,
                              smesh::mpi_type<idx_t>(), buffer);

    free(buffer);
  }

  return err == SMESH_SUCCESS ? SMESH_SUCCESS : SMESH_FAILURE;
}

} // namespace smesh
