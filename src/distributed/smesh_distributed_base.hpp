#ifndef SMESH_DISTRIBUTED_BASE_HPP
#define SMESH_DISTRIBUTED_BASE_HPP

#include "smesh_base.hpp"

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
}

#endif // SMESH_DISTRIBUTED_BASE_HPP