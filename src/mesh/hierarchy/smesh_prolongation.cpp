#include "smesh_prolongation.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_PROLONGATION(IDX_T, T)                       \
  template int hierarchical_prolongation<IDX_T, T>(                              \
      const enum ElemType, const enum ElemType, const ptrdiff_t,                 \
      const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, const int,        \
      const T *const SMESH_RESTRICT, T *const SMESH_RESTRICT);

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_PROLONGATION(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_PROLONGATION(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_PROLONGATION(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_PROLONGATION(i64, f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_PROLONGATION

