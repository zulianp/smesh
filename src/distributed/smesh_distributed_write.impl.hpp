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

template <typename idx_t>
int write_mapped_field(MPI_Comm comm, const Path &output_path,
                       const ptrdiff_t n_local, const ptrdiff_t n_global,
                       const idx_t *const mapping, MPI_Datatype data_type,
                       const void *const data_in) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  const uint8_t *const data = (const uint8_t *const)data_in;

  int type_size;
  SMESH_MPI_CATCH(MPI_Type_size(data_type, &type_size));

  const ptrdiff_t local_output_size_no_remainder = n_global / size;
  const ptrdiff_t begin = (n_global / size) * rank;

  ptrdiff_t local_output_size = local_output_size_no_remainder;
  if (rank == size - 1) {
    local_output_size = n_global - begin;
  }

  i64 *send_count = (i64 *)malloc((size_t)size * sizeof(i64));
  memset(send_count, 0, (size_t)size * sizeof(i64));

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];
    int dest_rank = std::min(
        size - 1, static_cast<int>(idx / local_output_size_no_remainder));
    send_count[dest_rank]++;
  }

  i64 *recv_count = (i64 *)malloc((size_t)size * sizeof(i64));
  SMESH_MPI_CATCH(MPI_Alltoall(send_count, 1, smesh::mpi_type<i64>(),
                               recv_count, 1, smesh::mpi_type<i64>(), comm));

  i64 *send_displs = (i64 *)malloc(((size_t)size + 1) * sizeof(i64));
  i64 *recv_displs = (i64 *)malloc(((size_t)size + 1) * sizeof(i64));
  ptrdiff_t *book_keeping = (ptrdiff_t *)calloc(size, sizeof(ptrdiff_t));

  send_displs[0] = 0;
  recv_displs[0] = 0;

  // Create data displacements for sending
  for (int i = 0; i < size; ++i) {
    send_displs[i + 1] = send_displs[i] + send_count[i];
  }

  // Create data displacements for receiving
  for (int i = 0; i < size; ++i) {
    recv_displs[i + 1] = recv_displs[i] + recv_count[i];
  }

  const ptrdiff_t total_recv = (ptrdiff_t)recv_displs[size];
  large_idx_t *send_list = (large_idx_t *)malloc(n_local * sizeof(large_idx_t));

  ptrdiff_t n_buff = std::max(n_local, local_output_size);
  uint8_t *send_data_and_final_storage = (uint8_t *)malloc(n_buff * type_size);

  // Pack data and indices
  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const large_idx_t idx = mapping[i];
    int dest_rank = std::min(
        size - 1, static_cast<int>(idx / local_output_size_no_remainder));
    SMESH_ASSERT(dest_rank < size);

    // Put index and data into buffers
    const ptrdiff_t offset = send_displs[dest_rank] + book_keeping[dest_rank];
    send_list[offset] = idx;
    memcpy((void *)&send_data_and_final_storage[offset * type_size],
           (void *)&data[i * type_size], type_size);

    book_keeping[dest_rank]++;
  }

  large_idx_t *recv_list = (large_idx_t *)malloc((size_t)total_recv * sizeof(large_idx_t));
  uint8_t *recv_data =
      (uint8_t *)malloc((size_t)total_recv * (size_t)type_size);

  ///////////////////////////////////
  // Send indices
  ///////////////////////////////////

  const i64 max_chunk_size = (i64)std::numeric_limits<i32>::max() / size;
  SMESH_MPI_CATCH(all_to_allv_64(send_list, send_count, send_displs, recv_list,
                                 recv_count, recv_displs, comm,
                                 max_chunk_size));

  ///////////////////////////////////
  // Send data
  ///////////////////////////////////

  SMESH_MPI_CATCH(all_to_allv_64_b(
      send_data_and_final_storage, send_count, send_displs, data_type,
      recv_data, recv_count, recv_displs, data_type, comm, max_chunk_size));

  int ret =
      array_write(comm, output_path.c_str(), data_type,
                  send_data_and_final_storage, local_output_size, n_global);

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
  return ret;
}

} // namespace smesh
