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
template <typename T> MPI_Datatype mpi_type();

template <> inline MPI_Datatype mpi_type<char>() { return MPI_CHAR; }
template <> inline MPI_Datatype mpi_type<f64>() { return MPI_DOUBLE; }
template <> inline MPI_Datatype mpi_type<f32>() { return MPI_FLOAT; }
template <> inline MPI_Datatype mpi_type<i16>() { return MPI_SHORT; }
template <> inline MPI_Datatype mpi_type<i32>() { return MPI_INT32_T; }
template <> inline MPI_Datatype mpi_type<i8>() { return MPI_CHAR; }
template <> inline MPI_Datatype mpi_type<i64>() { return MPI_INT64_T; }
// template <> inline MPI_Datatype mpi_type<long long>() { return MPI_LONG_LONG; }
template <> inline MPI_Datatype mpi_type<u16>() { return MPI_UNSIGNED_SHORT; }
template <> inline MPI_Datatype mpi_type<u64>() { return MPI_UNSIGNED_LONG_LONG; }
template <> inline MPI_Datatype mpi_type<u8>() { return MPI_UNSIGNED_CHAR; }
template <> inline MPI_Datatype mpi_type<u32>() { return MPI_UINT32_T; }


extern MPI_Datatype SMESH_MPI_F16;

int register_mpi_datatypes();
int unregister_mpi_datatypes();

template <> inline MPI_Datatype mpi_type<f16>() {
  // No standard MPI float16 type; use a registered 16-bit payload.
  return SMESH_MPI_F16;
}

inline MPI_Datatype
mpi_type_from_primitive_type(const enum PrimitiveType type) {
  switch (type) {
  case PrimitiveType::SMESH_FLOAT16:
    return mpi_type<f16>();
  case PrimitiveType::SMESH_FLOAT32:
    return MPI_FLOAT;
  case PrimitiveType::SMESH_FLOAT64:
    return MPI_DOUBLE;
  case PrimitiveType::SMESH_INT16:
    return MPI_SHORT;
  case PrimitiveType::SMESH_INT32:
    return MPI_INT;
  case PrimitiveType::SMESH_INT64:
    return MPI_LONG_LONG;
  case PrimitiveType::SMESH_UINT8:
    return MPI_UNSIGNED_CHAR;
  case PrimitiveType::SMESH_UINT16:
    return MPI_UNSIGNED_SHORT;
    //   case PrimitiveType::SMESH_UINT32:
    //     return MPI_UNSIGNED_INT;
  case PrimitiveType::SMESH_UINT64:
    return MPI_UNSIGNED_LONG_LONG;
  case PrimitiveType::SMESH_CHAR:
    return MPI_CHAR;
  default:
    SMESH_ERROR("Invalid primitive type: %d", type);
    return MPI_DATATYPE_NULL;
  }
}
} // namespace smesh

#endif // SMESH_DISTRIBUTED_BASE_HPP