#ifndef SMESH_DISTRIBUTED_BASE_HPP
#define SMESH_DISTRIBUTED_BASE_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

#include <mpi.h>
#include <stddef.h>
#include <stdint.h>

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
template <typename T> inline MPI_Datatype mpi_type() {
  using U = std::remove_cv_t<T>;

  if constexpr (std::is_same_v<U, char>)
    return MPI_CHAR;
  else if constexpr (std::is_same_v<U, std::int8_t>)
    return MPI_INT8_T;
  else if constexpr (std::is_same_v<U, std::uint8_t>)
    return MPI_UINT8_T;
  else if constexpr (std::is_same_v<U, std::int16_t>)
    return MPI_INT16_T;
  else if constexpr (std::is_same_v<U, std::uint16_t>)
    return MPI_UINT16_T;
  else if constexpr (std::is_same_v<U, std::int32_t>)
    return MPI_INT32_T;
  else if constexpr (std::is_same_v<U, std::uint32_t>)
    return MPI_UINT32_T;
  else if constexpr (std::is_same_v<U, std::int64_t>)
    return MPI_INT64_T;
  else if constexpr (std::is_same_v<U, std::uint64_t>)
    return MPI_UINT64_T;
  else if constexpr (std::is_same_v<U, float>)
    return MPI_FLOAT;
  else if constexpr (std::is_same_v<U, double>)
    return MPI_DOUBLE;
  else if constexpr (std::is_same_v<U, long double>)
    return MPI_LONG_DOUBLE;
    else if constexpr (std::is_same_v<U, ptrdiff_t>){
      if constexpr (sizeof(ptrdiff_t) == 8)
        return MPI_INT64_T;
      else if constexpr (sizeof(ptrdiff_t) == 4)
        return MPI_INT32_T;
      else if constexpr (sizeof(ptrdiff_t) == 2)
        return MPI_INT16_T;
      else if constexpr (sizeof(ptrdiff_t) == 1)
        return MPI_INT8_T;
      else {
        static_assert(!sizeof(ptrdiff_t), "Unsupported type in smesh::mpi_type<ptrdiff_t>()");
      }
    }
  else {
    static_assert(!sizeof(T), "Unsupported type in smesh::mpi_type<T>()");
  }
}

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