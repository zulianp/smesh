#ifndef SMESH_DISTRIBUTED_BASE_HPP
#define SMESH_DISTRIBUTED_BASE_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <mpi.h>

#define SMESH_MPI_CATCH(err)                                                   \
  {                                                                            \
    if (err != MPI_SUCCESS) {                                                  \
      char string_buff[4096];                                                  \
      int resultlen = 4096;                                                    \
      MPI_Error_string(err, string_buff, &resultlen);                          \
      fprintf(stderr, "MPI error: %s\n", string_buff);                         \
      fflush(stderr);                                                          \
      assert(0);                                                               \
      MPI_Abort(MPI_COMM_WORLD, -1);                                           \
    }                                                                          \
  }

namespace smesh {
    template <typename T>
    MPI_Datatype mpi_type();

    template <>
    MPI_Datatype mpi_type<float>() {
        return MPI_FLOAT;
    }
    template <>
    MPI_Datatype mpi_type<double>() {
        return MPI_DOUBLE;
    }
    template <>
    MPI_Datatype mpi_type<int>() {
        return MPI_INT;
    }
    template <>
    MPI_Datatype mpi_type<long long>() {
        return MPI_LONG_LONG;
    }

    template <>
    MPI_Datatype mpi_type<short>() {
        return MPI_SHORT;
    }

    // template <>
    // MPI_Datatype mpi_type<f16>() {
    //     // No standard MPI half type; treat as 16-bit payload.
    //     return MPI_UNSIGNED_SHORT;
    // }
}

#endif // SMESH_DISTRIBUTED_BASE_HPP