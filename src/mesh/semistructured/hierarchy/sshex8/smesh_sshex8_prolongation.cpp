#include "smesh_sshex8_prolongation.impl.hpp"


#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATION(IDX_T, REAL_T)          \
  template int sshex8_hierarchical_prolongation<IDX_T, REAL_T>(                  \
      int level, const ptrdiff_t nelements,                                   \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT elements,        \
      const int vec_size, const REAL_T *const SMESH_RESTRICT from,              \
      REAL_T *const SMESH_RESTRICT to)


#define SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATE(IDX_T, REAL_T)          \
  template int sshex8_prolongate<IDX_T, REAL_T>(                              \
      const ptrdiff_t nelements, const int from_level, const int from_level_stride, \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements, const int to_level, const int to_level_stride, \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements, const int vec_size, const REAL_T *const SMESH_RESTRICT from, \
      REAL_T *const SMESH_RESTRICT to)

namespace smesh {
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATION(i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATION(i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATION(i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATION(i64, f64);

    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATE(i32, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATE(i64, f32);
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATE(i32, f64);
    SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATE(i64, f64);
}

#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATION
#undef SMESH_EXPLICIT_INSTANTIATE_SSHEX8_PROLONGATE
