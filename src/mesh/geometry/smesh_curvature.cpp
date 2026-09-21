#include "smesh_curvature.impl.hpp"

namespace smesh {

template int mesh_node_curvature<i32, f32>(const enum ElemType,
                                           const ptrdiff_t,
                                           const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const int,
                                           const ptrdiff_t,
                                           const f32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const uint8_t *const SMESH_RESTRICT,
                                           f32 *const SMESH_RESTRICT);
template int mesh_node_curvature<i64, f32>(const enum ElemType,
                                           const ptrdiff_t,
                                           const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const int,
                                           const ptrdiff_t,
                                           const f32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const uint8_t *const SMESH_RESTRICT,
                                           f32 *const SMESH_RESTRICT);
template int mesh_node_curvature<i32, f64>(const enum ElemType,
                                           const ptrdiff_t,
                                           const i32 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const int,
                                           const ptrdiff_t,
                                           const f64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const uint8_t *const SMESH_RESTRICT,
                                           f64 *const SMESH_RESTRICT);
template int mesh_node_curvature<i64, f64>(const enum ElemType,
                                           const ptrdiff_t,
                                           const i64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const int,
                                           const ptrdiff_t,
                                           const f64 *const SMESH_RESTRICT *const SMESH_RESTRICT,
                                           const uint8_t *const SMESH_RESTRICT,
                                           f64 *const SMESH_RESTRICT);

template int mesh_grade_size_2to1<i32, i32, f32>(const ptrdiff_t,
                                                 const i32 *const SMESH_RESTRICT,
                                                 const i32 *const SMESH_RESTRICT,
                                                 const int,
                                                 f32 *const SMESH_RESTRICT);
template int mesh_grade_size_2to1<i64, i32, f32>(const ptrdiff_t,
                                                 const i32 *const SMESH_RESTRICT,
                                                 const i64 *const SMESH_RESTRICT,
                                                 const int,
                                                 f32 *const SMESH_RESTRICT);
template int mesh_grade_size_2to1<i64, i64, f32>(const ptrdiff_t,
                                                 const i64 *const SMESH_RESTRICT,
                                                 const i64 *const SMESH_RESTRICT,
                                                 const int,
                                                 f32 *const SMESH_RESTRICT);
template int mesh_grade_size_2to1<i32, i32, f64>(const ptrdiff_t,
                                                 const i32 *const SMESH_RESTRICT,
                                                 const i32 *const SMESH_RESTRICT,
                                                 const int,
                                                 f64 *const SMESH_RESTRICT);
template int mesh_grade_size_2to1<i64, i32, f64>(const ptrdiff_t,
                                                 const i32 *const SMESH_RESTRICT,
                                                 const i64 *const SMESH_RESTRICT,
                                                 const int,
                                                 f64 *const SMESH_RESTRICT);
template int mesh_grade_size_2to1<i64, i64, f64>(const ptrdiff_t,
                                                 const i64 *const SMESH_RESTRICT,
                                                 const i64 *const SMESH_RESTRICT,
                                                 const int,
                                                 f64 *const SMESH_RESTRICT);

template int mesh_size_from_curvature<f32>(const ptrdiff_t,
                                           const f32 *const SMESH_RESTRICT,
                                           const f32,
                                           const f32,
                                           const f32,
                                           const f32,
                                           f32 *const SMESH_RESTRICT);
template int mesh_size_from_curvature<f64>(const ptrdiff_t,
                                           const f64 *const SMESH_RESTRICT,
                                           const f64,
                                           const f64,
                                           const f64,
                                           const f64,
                                           f64 *const SMESH_RESTRICT);

template int mesh_size_from_curvature_error<f32>(const ptrdiff_t,
                                                 const f32 *const SMESH_RESTRICT,
                                                 const f32,
                                                 const f32,
                                                 const f32,
                                                 const f32,
                                                 const f32,
                                                 f32 *const SMESH_RESTRICT);
template int mesh_size_from_curvature_error<f64>(const ptrdiff_t,
                                                 const f64 *const SMESH_RESTRICT,
                                                 const f64,
                                                 const f64,
                                                 const f64,
                                                 const f64,
                                                 const f64,
                                                 f64 *const SMESH_RESTRICT);

template void mesh_mark_nodes_from_edges<i32>(const ptrdiff_t,
                                              const ptrdiff_t,
                                              const i32 *const SMESH_RESTRICT,
                                              const i32 *const SMESH_RESTRICT,
                                              uint8_t *const SMESH_RESTRICT);
template void mesh_mark_nodes_from_edges<i64>(const ptrdiff_t,
                                              const ptrdiff_t,
                                              const i64 *const SMESH_RESTRICT,
                                              const i64 *const SMESH_RESTRICT,
                                              uint8_t *const SMESH_RESTRICT);

}  // namespace smesh
