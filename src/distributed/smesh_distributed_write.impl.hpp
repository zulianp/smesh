#define SMESH_DISTRIBUTED_WRITE_IMPL_HPP

#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#include "smesh_types.hpp"
#include "matrixio_array.h"

#include <chrono>
#include <fstream>

namespace smesh {

// #region agent log
static inline void smesh_dbglog_write(const char *location, const char *message,
                                      const std::string &data_json) {
  using namespace std::chrono;
  const auto ts =
      duration_cast<milliseconds>(system_clock::now().time_since_epoch())
          .count();
  std::ofstream os("/Users/patrickzulian/Desktop/code/smesh/.cursor/debug.log",
                   std::ios::app);
  if (!os.good())
    return;
  os << "{\"id\":\"log_" << ts << "_" << location << "\","
     << "\"timestamp\":" << ts << ","
     << "\"runId\":\"pre\","
     << "\"hypothesisId\":\"W\","
     << "\"location\":\"" << location << "\","
     << "\"message\":\"" << message << "\","
     << "\"data\":" << data_json << "}\n";
}
// #endregion

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
  if (ext == ".raw") {
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
  } else {
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

  int *send_count = (int *)malloc((size) * sizeof(int));
  memset(send_count, 0, (size) * sizeof(int));

  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];
    int dest_rank = std::min(size - 1, static_cast<int>(idx / local_output_size_no_remainder));
    send_count[dest_rank]++;
  }

  int *recv_count = (int *)malloc((size) * sizeof(int));
  SMESH_MPI_CATCH(MPI_Alltoall(send_count, 1, smesh::mpi_type<idx_t>(),
                               recv_count, 1, smesh::mpi_type<idx_t>(), comm));

  int *send_displs = (int *)malloc(size * sizeof(int));
  int *recv_displs = (int *)malloc(size * sizeof(int));
  ptrdiff_t *book_keeping = (ptrdiff_t *)calloc(size, sizeof(ptrdiff_t));

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

//   const ptrdiff_t total_recv = recv_displs[size - 1] + recv_count[size - 1];
  const ptrdiff_t total_recv = (ptrdiff_t)recv_displs[size - 1] + (ptrdiff_t)recv_count[size - 1];
  // #region agent log
  {
    if (total_recv > local_output_size) {
      smesh_dbglog_write(
          "smesh_distributed_write.impl.hpp:write_mapped_field",
          "recv larger than output segment",
          std::string("{\"rank\":") + std::to_string(rank) +
              ",\"size\":" + std::to_string(size) +
              ",\"n_local\":" + std::to_string((long long)n_local) +
              ",\"n_global\":" + std::to_string((long long)n_global) +
              ",\"begin\":" + std::to_string((long long)begin) +
              ",\"local_output_size\":" +
              std::to_string((long long)local_output_size) +
              ",\"total_recv\":" + std::to_string((long long)total_recv) + "}");
    }
  }
  // #endregion

  idx_t *send_list = (idx_t *)malloc(n_local * sizeof(idx_t));

  ptrdiff_t n_buff = std::max(n_local, local_output_size);
  uint8_t *send_data_and_final_storage = (uint8_t *)malloc(n_buff * type_size);

  // Pack data and indices
  for (ptrdiff_t i = 0; i < n_local; ++i) {
    const idx_t idx = mapping[i];
    int dest_rank = std::min(size - 1, static_cast<int>(idx / local_output_size_no_remainder));
    SMESH_ASSERT(dest_rank < size);

    // Put index and data into buffers
    const ptrdiff_t offset = send_displs[dest_rank] + book_keeping[dest_rank];
    send_list[offset] = idx;
    memcpy((void *)&send_data_and_final_storage[offset * type_size],
           (void *)&data[i * type_size], type_size);

    book_keeping[dest_rank]++;
  }

  idx_t *recv_list = (idx_t *)malloc(local_output_size * sizeof(idx_t));
  uint8_t *recv_data = (uint8_t *)malloc(local_output_size * type_size);

  ///////////////////////////////////
  // Send indices
  ///////////////////////////////////

  SMESH_MPI_CATCH(MPI_Alltoallv(send_list, send_count, send_displs,
                                smesh::mpi_type<idx_t>(), recv_list, recv_count,
                                recv_displs, smesh::mpi_type<idx_t>(), comm));

  ///////////////////////////////////
  // Send data
  ///////////////////////////////////

  SMESH_MPI_CATCH(MPI_Alltoallv(send_data_and_final_storage, send_count,
                                send_displs, data_type, recv_data, recv_count,
                                recv_displs, data_type, comm));

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

  int ret = array_write_convert_from_extension(comm, output_path,
                                               send_data_and_final_storage,
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
  return ret;
}

} // namespace smesh