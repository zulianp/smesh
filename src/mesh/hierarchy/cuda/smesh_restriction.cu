namespace smesh {

template <typename From, typename To>
__global__ void cu_macrotet4_to_tet4_restriction_kernel(
    const ptrdiff_t nelements, idx_t **const SMESH_RESTRICT elements,
    const uint16_t *const SMESH_RESTRICT e2n_count, const int vec_size,
    const ptrdiff_t from_stride, const From *const SMESH_RESTRICT from,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to) {
  for (ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x; e < nelements;
       e += blockDim.x * gridDim.x) {
    const idx_t i0 = elements[0][e];
    const idx_t i1 = elements[1][e];
    const idx_t i2 = elements[2][e];
    const idx_t i3 = elements[3][e];

    // P2
    const idx_t i4 = elements[4][e];
    const idx_t i5 = elements[5][e];
    const idx_t i6 = elements[6][e];
    const idx_t i7 = elements[7][e];
    const idx_t i8 = elements[8][e];
    const idx_t i9 = elements[9][e];

    for (int v = 0; v < vec_size; v++) {
      const ptrdiff_t ii0 = i0 * vec_size + v;
      const ptrdiff_t ii1 = i1 * vec_size + v;
      const ptrdiff_t ii2 = i2 * vec_size + v;
      const ptrdiff_t ii3 = i3 * vec_size + v;
      const ptrdiff_t ii4 = i4 * vec_size + v;
      const ptrdiff_t ii5 = i5 * vec_size + v;
      const ptrdiff_t ii6 = i6 * vec_size + v;
      const ptrdiff_t ii7 = i7 * vec_size + v;
      const ptrdiff_t ii8 = i8 * vec_size + v;
      const ptrdiff_t ii9 = i9 * vec_size + v;

      const T to0 = from[ii0 * from_stride] / e2n_count[i0] +
                    from[ii4 * from_stride] * (0.5 / e2n_count[i4]) +
                    from[ii6 * from_stride] * (0.5 / e2n_count[i6]) +
                    from[ii7 * from_stride] * (0.5 / e2n_count[i7]);

      const T to1 = from[ii1 * from_stride] / e2n_count[i1] +
                    from[ii5 * from_stride] * (0.5 / e2n_count[i5]) +
                    from[ii8 * from_stride] * (0.5 / e2n_count[i8]) +
                    from[ii4 * from_stride] * (0.5 / e2n_count[i4]);

      const T to2 = from[ii2 * from_stride] / e2n_count[i2] +
                    from[ii9 * from_stride] * (0.5 / e2n_count[i9]) +
                    from[ii5 * from_stride] * (0.5 / e2n_count[i5]) +
                    from[ii6 * from_stride] * (0.5 / e2n_count[i6]);

      const T to3 = from[ii3 * from_stride] / e2n_count[i3] +
                    from[ii7 * from_stride] * (0.5 / e2n_count[i7]) +
                    from[ii8 * from_stride] * (0.5 / e2n_count[i8]) +
                    from[ii9 * from_stride] * (0.5 / e2n_count[i9]);

      atomicAdd(&to[ii0 * to_stride], to0);
      atomicAdd(&to[ii1 * to_stride], to1);
      atomicAdd(&to[ii2 * to_stride], to2);
      atomicAdd(&to[ii3 * to_stride], to3);
    }
  }
}

template <typename From, typename To>
static int cu_macrotet4_to_tet4_restriction_tpl(
    const ptrdiff_t nelements, idx_t **const SMESH_RESTRICT elements,
    const uint16_t *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const ptrdiff_t from_stride,
    const From *const SMESH_RESTRICT from, const ptrdiff_t to_stride,
    To *const SMESH_RESTRICT to, void *stream) {
  // Hand tuned
  int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
  {
    int min_grid_size;
    cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &block_size,
        cu_macrotet4_to_tet4_restriction_kernel<From, To>, 0, 0);
  }
#endif // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

  ptrdiff_t n_blocks =
      MAX(ptrdiff_t(1), (nelements + block_size - 1) / block_size);

  if (stream) {
    cudaStream_t s = *static_cast<cudaStream_t *>(stream);

    cu_macrotet4_to_tet4_restriction_kernel<From, To>
        <<<n_blocks, block_size, 0, s>>>(
            nelements, elements, element_to_node_incidence_count, vec_size,
            from_stride, from, to_stride, to);
  } else {
    cu_macrotet4_to_tet4_restriction_kernel<From, To>
        <<<n_blocks, block_size, 0>>>(nelements, elements,
                                      element_to_node_incidence_count, vec_size,
                                      from_stride, from, to_stride, to);
  }

  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}

int cu_macrotet4_to_tet4_restriction(
    const ptrdiff_t nelements, idx_t **const SMESH_RESTRICT elements,
    const uint16_t *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const enum RealType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum RealType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream) {
  SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
  if (from_type != to_type) {
    SMESH_ERROR(
        "Unsupported element pair for hierarchical_restriction %d, %d\n",
        from_type, to_type);
    return SMESH_FAILURE;
  }

  switch (from_type) {
  case SMESH_REAL_DEFAULT: {
    return cu_macrotet4_to_tet4_restriction_tpl(
        nelements, elements, element_to_node_incidence_count, vec_size,
        from_stride, (real_t *)from, to_stride, (real_t *)to, stream);
  }
  case SMESH_FLOAT32: {
    return cu_macrotet4_to_tet4_restriction_tpl(
        nelements, elements, element_to_node_incidence_count, vec_size,
        from_stride, (f32 *)from, to_stride, (f32 *)to, stream);
  }
  case SMESH_FLOAT64: {
    return cu_macrotet4_to_tet4_restriction_tpl(
        nelements, elements, element_to_node_incidence_count, vec_size,
        from_stride, (f64 *)from, to_stride, (f64 *)to, stream);
  }
  default: {
    SMESH_ERROR(
        "[Error] cu_tet4_to_macrotet4_prolongation_tpl: not implemented "
        "for type "
        "%s "
        "(code %d)\n",
        real_type_to_string(from_type), from_type);
    return SMESH_FAILURE;
  }
  }
}
} // namespace smesh