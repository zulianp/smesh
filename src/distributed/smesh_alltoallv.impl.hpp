#ifndef SMESH_ALLTOALLV_IMPL_HPP
#define SMESH_ALLTOALLV_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_communicator.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"
#include "smesh_alloc.hpp"

#include <string>
#include <limits>
#include <mpi.h>
#include <stddef.h>
#include <algorithm>
#include <cstring>

namespace smesh {

inline int
all_to_allv_64_b(const void *send_elements, const i64 *large_send_count,
                 const i64 *large_send_displs, MPI_Datatype send_datatype,
                 void *recv_elements, const i64 *large_recv_count,
                 const i64 *large_recv_displs, MPI_Datatype recv_datatype,
                 MPI_Comm comm, i64 max_chunk_size) {

  using byte_t = unsigned char;

  SMESH_TRACE_SCOPE("all_to_allv_64");
  SMESH_ASSERT(max_chunk_size > 0);

  int size;
  MPI_Comm_size(comm, &size);
  const i64 i32_max = (i64)std::numeric_limits<int>::max();
  // const i64 i32_max  = max_chunk_size * size;

  i64 local_max_peer_count = 0;
  i64 local_send_buffer_size = 0;
  i64 local_recv_buffer_size = 0;
  for (int r = 0; r < size; r++) {
    local_max_peer_count = std::max(local_max_peer_count, large_send_count[r]);
    local_max_peer_count = std::max(local_max_peer_count, large_recv_count[r]);
    local_send_buffer_size += std::min(max_chunk_size, large_send_count[r]);
    local_recv_buffer_size += std::min(max_chunk_size, large_recv_count[r]);
  }

  i64 all_count = large_send_displs[size - 1] + large_send_count[size- 1];
  all_count = std::max(all_count, large_recv_displs[size - 1] + large_recv_count[size - 1]);

  i64 global_max_peer_count = 0;
  SMESH_MPI_CATCH(MPI_Allreduce(&local_max_peer_count, &global_max_peer_count,
                                1, mpi_type<i64>(), MPI_MAX, comm));

  const int local_fits_i32 = all_count <= i32_max;
  int global_fits_i32 = 0;
  SMESH_MPI_CATCH(MPI_Allreduce(&local_fits_i32, &global_fits_i32, 1, MPI_INT,
                                MPI_MIN, comm));

  const i64 n_rounds = div_round_up(global_max_peer_count, max_chunk_size);

  int *send_displs = (int *)SMESH_CALLOC(size + 1, sizeof(int));
  int *send_count = (int *)SMESH_CALLOC(size, sizeof(int));
  int *recv_displs = (int *)SMESH_CALLOC((size + 1), sizeof(int));
  int *recv_count = (int *)SMESH_CALLOC(size, sizeof(int));

  if (global_fits_i32) {
    for (int r = 0; r < size; r++) {
      send_count[r] = (int)large_send_count[r];
      recv_count[r] = (int)large_recv_count[r];
      send_displs[r] = (int)large_send_displs[r];
      recv_displs[r] = (int)large_recv_displs[r];
    }

    SMESH_MPI_CATCH(MPI_Alltoallv(send_elements, send_count, send_displs,
                                  send_datatype, recv_elements, recv_count,
                                  recv_displs, recv_datatype, comm));
  } else {
    SMESH_ASSERT(max_chunk_size <= i32_max);
    SMESH_ASSERT((i64)size * max_chunk_size <= i32_max);

    i64 *send_offsets = (i64 *)SMESH_CALLOC(size, sizeof(i64));
    i64 *recv_offsets = (i64 *)SMESH_CALLOC(size, sizeof(i64));

    int send_type_size;
    int recv_type_size;
    SMESH_MPI_CATCH(MPI_Type_size(send_datatype, &send_type_size));
    SMESH_MPI_CATCH(MPI_Type_size(recv_datatype, &recv_type_size));

    byte_t *send_buffer =
        (byte_t *)SMESH_ALLOC((size_t)local_send_buffer_size * send_type_size);
    byte_t *recv_buffer =
        (byte_t *)SMESH_ALLOC((size_t)local_recv_buffer_size * recv_type_size);

    for (i64 round = 0; round < n_rounds; round++) {
      send_displs[0] = 0;
      recv_displs[0] = 0;

      for (int r = 0; r < size; r++) {
        const i64 send_remaining = large_send_count[r] - send_offsets[r];
        const int send_chunk =
            (send_remaining > 0) ? (int)std::min(max_chunk_size, send_remaining)
                                 : 0;
        send_count[r] = send_chunk;
        send_displs[r + 1] = send_displs[r] + send_chunk;

        if (send_chunk) {
          const i64 send_begin = large_send_displs[r] + send_offsets[r];
          memcpy(&send_buffer[send_displs[r] * send_type_size],
                 &((byte_t *)send_elements)[send_begin * send_type_size],
                 (size_t)send_chunk * send_type_size);
          send_offsets[r] += send_chunk;
        }

        const i64 recv_remaining = large_recv_count[r] - recv_offsets[r];
        const int recv_chunk =
            (recv_remaining > 0) ? (int)std::min(max_chunk_size, recv_remaining)
                                 : 0;
        recv_count[r] = recv_chunk;
        recv_displs[r + 1] = recv_displs[r] + recv_chunk;
      }

      SMESH_MPI_CATCH(MPI_Alltoallv(send_buffer, send_count, send_displs,
                                    send_datatype, recv_buffer, recv_count,
                                    recv_displs, recv_datatype, comm));

      for (int r = 0; r < size; r++) {
        const int recv_chunk = recv_count[r];
        if (!recv_chunk)
          continue;

        const i64 recv_begin = large_recv_displs[r] + recv_offsets[r];
        memcpy(&((byte_t *)recv_elements)[recv_begin * recv_type_size],
               &recv_buffer[recv_displs[r] * recv_type_size],
               (size_t)recv_chunk * recv_type_size);
        recv_offsets[r] += recv_chunk;
      }
    }

    SMESH_FREE(send_buffer);
    SMESH_FREE(recv_buffer);
    SMESH_FREE(send_offsets);
    SMESH_FREE(recv_offsets);
  }

  SMESH_FREE(send_displs);
  SMESH_FREE(send_count);
  SMESH_FREE(recv_displs);
  SMESH_FREE(recv_count);
  return SMESH_SUCCESS;
}

template <typename T>
int all_to_allv_64(const T *send_elements, const i64 *large_send_count,
                   const i64 *large_send_displs, T *recv_elements,
                   const i64 *large_recv_count, const i64 *large_recv_displs,
                   MPI_Comm comm, i64 max_chunk_size) {

  return all_to_allv_64_b(send_elements, large_send_count, large_send_displs,
                          mpi_type<T>(), recv_elements, large_recv_count,
                          large_recv_displs, mpi_type<T>(), comm,
                          max_chunk_size);
}

inline int all_to_allv_64_bv(const void *send_elements,
                             const i64 *large_send_count,
                             const i64 *large_send_displs,
                             MPI_Datatype scalar_send_datatype,
                             void *recv_elements,
                             const i64 *large_recv_count,
                             const i64 *large_recv_displs,
                             MPI_Datatype scalar_recv_datatype,
                             const ptrdiff_t block_size, MPI_Comm comm,
                             i64 max_chunk_size) {
  SMESH_ASSERT(block_size > 0);

  if (block_size == 1) {
    return all_to_allv_64_b(send_elements, large_send_count, large_send_displs,
                            scalar_send_datatype, recv_elements,
                            large_recv_count, large_recv_displs,
                            scalar_recv_datatype, comm, max_chunk_size);
  }

  SMESH_ASSERT(block_size <= std::numeric_limits<int>::max());

  MPI_Datatype send_datatype = MPI_DATATYPE_NULL;
  MPI_Datatype recv_datatype = MPI_DATATYPE_NULL;

  SMESH_MPI_CATCH(
      MPI_Type_contiguous((int)block_size, scalar_send_datatype, &send_datatype));
  SMESH_MPI_CATCH(MPI_Type_commit(&send_datatype));

  if (scalar_send_datatype == scalar_recv_datatype) {
    recv_datatype = send_datatype;
  } else {
    SMESH_MPI_CATCH(MPI_Type_contiguous((int)block_size, scalar_recv_datatype,
                                        &recv_datatype));
    SMESH_MPI_CATCH(MPI_Type_commit(&recv_datatype));
  }

  const int err =
      all_to_allv_64_b(send_elements, large_send_count, large_send_displs,
                       send_datatype, recv_elements, large_recv_count,
                       large_recv_displs, recv_datatype, comm, max_chunk_size);

  if (recv_datatype != send_datatype) {
    SMESH_MPI_CATCH(MPI_Type_free(&recv_datatype));
  }
  SMESH_MPI_CATCH(MPI_Type_free(&send_datatype));
  return err;
}

template <typename T>
int all_to_allv_64v(const T *send_elements, const i64 *large_send_count,
                    const i64 *large_send_displs, T *recv_elements,
                    const i64 *large_recv_count,
                    const i64 *large_recv_displs, const ptrdiff_t block_size,
                    MPI_Comm comm, i64 max_chunk_size) {
  return all_to_allv_64_bv(send_elements, large_send_count, large_send_displs,
                           mpi_type<T>(), recv_elements, large_recv_count,
                           large_recv_displs, mpi_type<T>(), block_size, comm,
                           max_chunk_size);
}

} // namespace smesh
#endif // SMESH_ALLTOALLV_IMPL_HPP
