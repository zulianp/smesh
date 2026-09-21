#include "smesh_quality.impl.hpp"

namespace smesh {

template f32 mesh_elem_mean_ratio<i32, f32>(const enum ElemType,
                                            const int,
                                            const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const f32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const ptrdiff_t);
template f32 mesh_elem_mean_ratio<i64, f32>(const enum ElemType,
                                            const int,
                                            const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const f32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const ptrdiff_t);
template f64 mesh_elem_mean_ratio<i32, f64>(const enum ElemType,
                                            const int,
                                            const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const f64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const ptrdiff_t);
template f64 mesh_elem_mean_ratio<i64, f64>(const enum ElemType,
                                            const int,
                                            const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const f64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const ptrdiff_t);

template int mesh_element_quality<i32, f32>(const enum ElemType,
                                            const ptrdiff_t,
                                            const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const int,
                                            const f32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            f32 *const SMESH_RESTRICT);
template int mesh_element_quality<i64, f32>(const enum ElemType,
                                            const ptrdiff_t,
                                            const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const int,
                                            const f32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            f32 *const SMESH_RESTRICT);
template int mesh_element_quality<i32, f64>(const enum ElemType,
                                            const ptrdiff_t,
                                            const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const int,
                                            const f64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            f64 *const SMESH_RESTRICT);
template int mesh_element_quality<i64, f64>(const enum ElemType,
                                            const ptrdiff_t,
                                            const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            const int,
                                            const f64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                            f64 *const SMESH_RESTRICT);

template f32 mesh_quality_min<f32>(const ptrdiff_t, const f32 *const SMESH_RESTRICT);
template f64 mesh_quality_min<f64>(const ptrdiff_t, const f64 *const SMESH_RESTRICT);

template void mesh_clamp_to_band<f32>(const int,
                                      f32 *,
                                      f32 *,
                                      f32 *,
                                      const f32,
                                      const f32,
                                      const f32,
                                      const f32,
                                      const f32,
                                      const f32,
                                      const f32,
                                      const f32);
template void mesh_clamp_to_band<f64>(const int,
                                      f64 *,
                                      f64 *,
                                      f64 *,
                                      const f64,
                                      const f64,
                                      const f64,
                                      const f64,
                                      const f64,
                                      const f64,
                                      const f64,
                                      const f64);

#define SMESH_EXPLICIT_INSTANTIATE_SURF_GEOM(IDX_T, GEOM_T)                                          \
    template int mesh_vertex_normals_from_faces<IDX_T, GEOM_T>(                                      \
            const enum ElemType,                                                                     \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                 \
            const int,                                                                               \
            const ptrdiff_t,                                                                         \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            GEOM_T *const SMESH_RESTRICT,                                                            \
            GEOM_T *const SMESH_RESTRICT,                                                            \
            GEOM_T *const SMESH_RESTRICT);                                                           \
    template int mesh_project_to_surface<IDX_T, GEOM_T>(                                             \
            const enum ElemType,                                                                     \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                 \
            const int,                                                                               \
            const ptrdiff_t,                                                                         \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            const ptrdiff_t,                                                                         \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                      \
            const GEOM_T);                                                                           \
    template int mesh_project_to_segments<IDX_T, GEOM_T>(                                            \
            const int,                                                                               \
            const ptrdiff_t,                                                                         \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const IDX_T *const SMESH_RESTRICT,                                                       \
            const ptrdiff_t,                                                                         \
            const GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                \
            const ptrdiff_t,                                                                         \
            const uint8_t *const SMESH_RESTRICT,                                                     \
            GEOM_T *const SMESH_RESTRICT *const SMESH_RESTRICT,                                      \
            const GEOM_T);

SMESH_EXPLICIT_INSTANTIATE_SURF_GEOM(i32, f32)
SMESH_EXPLICIT_INSTANTIATE_SURF_GEOM(i64, f32)
SMESH_EXPLICIT_INSTANTIATE_SURF_GEOM(i32, f64)
SMESH_EXPLICIT_INSTANTIATE_SURF_GEOM(i64, f64)

}  // namespace smesh
