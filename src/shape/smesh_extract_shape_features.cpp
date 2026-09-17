#include "smesh_extract_shape_features.impl.hpp"

#include "smesh_types.hpp"

#define SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(IDX_T, GEOM_T)          \
  template int extract_sharp_edges<IDX_T, GEOM_T, IDX_T>(                      \
      enum ElemType, ptrdiff_t,                                                \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, ptrdiff_t,            \
      GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                      \
      const IDX_T *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,    \
      GEOM_T, ptrdiff_t *, IDX_T **, IDX_T **)

#define SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHAP_FEATURES(IDX_T)                \
  template int extract_sharp_corners<IDX_T>(                            \
      ptrdiff_t, ptrdiff_t, IDX_T *const SMESH_RESTRICT,                       \
      IDX_T *const SMESH_RESTRICT, ptrdiff_t *, IDX_T **, int);                \
  template int extract_disconnected_faces<IDX_T, IDX_T>(                \
      enum ElemType, ptrdiff_t, ptrdiff_t,                                     \
      IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT, ptrdiff_t,            \
      const IDX_T *const SMESH_RESTRICT, const IDX_T *const SMESH_RESTRICT,    \
      ptrdiff_t *, IDX_T **)

namespace smesh {
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHAP_FEATURES(i16);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHAP_FEATURES(i32);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHAP_FEATURES(i64);

SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(i16, f32);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(i32, f32);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(i64, f32);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(i16, f64);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(i32, f64);
SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES(i64, f64);
} // namespace smesh

#undef SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHAP_FEATURES
#undef SMESH_EXPLICIT_INSTANTIATE_EXTRACT_SHARP_EDGES
