#include "smesh_sidesets.impl.hpp"
#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_SIDESETS(T)                                 \
  template int extract_nodeset_from_sideset<T, T>(                             \
      const enum ElemType,                                                     \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT, const ptrdiff_t,    \
      const T *const SMESH_RESTRICT, const i16 *const SMESH_RESTRICT,          \
      ptrdiff_t *, T **SMESH_RESTRICT);                                        \
  template int extract_nodeset_from_sidesets<T, T>(                            \
      ptrdiff_t, const enum ElemType[], T **const SMESH_RESTRICT[],            \
      const ptrdiff_t[], const T *const SMESH_RESTRICT[],                      \
      const i16 *const SMESH_RESTRICT[], ptrdiff_t *, T **SMESH_RESTRICT);     \
  template int extract_surface_from_sideset<T, T>(                             \
      const enum ElemType element_type,                                        \
      const T *const SMESH_RESTRICT *const SMESH_RESTRICT elems,               \
      const ptrdiff_t n_surf_elements,                                         \
      const T *const SMESH_RESTRICT parent_element,                            \
      const i16 *const SMESH_RESTRICT side_idx,                                \
      T *const SMESH_RESTRICT *const SMESH_RESTRICT sides);

      

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_SIDESETS(i32);
SMESH_EXPLICIT_INSTANTIATE_SIDESETS(i64);
} // namespace smesh
