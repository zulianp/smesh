#include "smesh_distributed_write.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(FILE_TYPE, TYPE)        \
  template int array_write_convert<FILE_TYPE, TYPE>(                           \
      MPI_Comm comm, const Path &path, const TYPE *const SMESH_RESTRICT data,  \
      const ptrdiff_t n_local_elements, const ptrdiff_t n_global_elements);

#define SMESH_EXPLICIT_INSTANTIATE(TYPE)                                       \
  template int array_write_convert_from_extension<TYPE>(                       \
      MPI_Comm comm, const Path &path, const TYPE *const SMESH_RESTRICT data,  \
      const ptrdiff_t n_local_elements, const ptrdiff_t n_global_elements);    \
  template int write_mapped_field<TYPE>(                                       \
      MPI_Comm comm, const Path &output_path, const ptrdiff_t n_local,         \
      const ptrdiff_t n_global, const TYPE *const mapping,                     \
      MPI_Datatype data_type, const void *const data_in);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(f16, f32);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(f32, f64);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(f16, f64);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(f64, f64);

SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(i16, i32);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT(i64, i64);

SMESH_EXPLICIT_INSTANTIATE(i8);
SMESH_EXPLICIT_INSTANTIATE(i16);
SMESH_EXPLICIT_INSTANTIATE(u8);
SMESH_EXPLICIT_INSTANTIATE(u16);
SMESH_EXPLICIT_INSTANTIATE(u32);
SMESH_EXPLICIT_INSTANTIATE(u64);
SMESH_EXPLICIT_INSTANTIATE(i32);
SMESH_EXPLICIT_INSTANTIATE(i64);
SMESH_EXPLICIT_INSTANTIATE(mask_t);

SMESH_EXPLICIT_INSTANTIATE(f16);
SMESH_EXPLICIT_INSTANTIATE(f32);
SMESH_EXPLICIT_INSTANTIATE(f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_ARRAY_WRITE_CONVERT
#undef SMESH_EXPLICIT_INSTANTIATE