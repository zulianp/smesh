
#define PROLONGATE_IN_KERNEL_BASIS_VERSION  // The other version needs debugging
#ifdef PROLONGATE_IN_KERNEL_BASIS_VERSION

// Even TO sub-elements are used to interpolate from FROM sub-elements
template <typename From, typename To>
__global__ void cu_sshex8_prolongate_kernel(const ptrdiff_t                 nelements,
                                            const int                       from_level,
                                            const int                       from_level_stride,
                                            idx_t **const SFEM_RESTRICT     from_elements,
                                            const int                       to_level,
                                            const int                       to_level_stride,
                                            idx_t **const SFEM_RESTRICT     to_elements,
                                            const int                       vec_size,
                                            const enum RealType             from_type,
                                            const ptrdiff_t                 from_stride,
                                            const From *const SFEM_RESTRICT from,
                                            const enum RealType             to_type,
                                            const ptrdiff_t                 to_stride,
                                            To *const SFEM_RESTRICT         to) {
    static_assert(TILE_SIZE == 8, "This only works with tile size 8!");

    // Uunsigned char necessary for multiple instantiations
    extern __shared__ unsigned char cu_buff[];

    const int step_factor = to_level / from_level;

    // Tile number in group
    const int tile    = threadIdx.x >> 3;   // same as threadIdx.x / 8
    const int n_tiles = blockDim.x >> 3;    // same as blockDim.x / 8
    const int sub_idx = threadIdx.x & 0x7;  // same as threadIdx.x % 8

    From *in = (From *)&cu_buff[tile * TILE_SIZE * sizeof(From)];

    // hex8 idx
    const int xi = sub_idx & 0x1;         // equivalent to sub_idx % 2
    const int yi = (sub_idx >> 1) & 0x1;  // equivalent to (sub_idx / 2) % 2
    const int zi = (sub_idx >> 2);        // equivalent to sub_idx / 4
    assert(n_tiles * TILE_SIZE == blockDim.x);

    // 1 macro element per tile
    const ptrdiff_t e = blockIdx.x * n_tiles + tile;

    const To between_h = (From)from_level / (To)to_level;

    const int to_even     = is_even(to_level);
    const int from_nloops = from_level + to_even;
    const int to_nloops   = to_level + to_even;

    // Vector loop
    for (int d = 0; d < vec_size; d++) {
        // loop on all FROM micro elements
        for (int from_zi = 0; from_zi < from_nloops; from_zi++) {
            const int off_from_zi = (from_zi + zi);

            for (int from_yi = 0; from_yi < from_nloops; from_yi++) {
                const int off_from_yi = (from_yi + yi);

                for (int from_xi = 0; from_xi < from_nloops; from_xi++) {
                    const int off_from_xi = (from_xi + xi);

                    const bool from_exists =
                            e < nelements && off_from_zi <= from_level && off_from_yi <= from_level && off_from_xi <= from_level;

                    // Wait for shared memory transactions to be finished
                    __syncwarp();

                    // Gather
                    if (from_exists) {
                        const int idx_from = cu_sshex8_lidx(from_level * from_level_stride,
                                                            off_from_xi * from_level_stride,
                                                            off_from_yi * from_level_stride,
                                                            off_from_zi * from_level_stride);

                        const idx_t     gidx = from_elements[idx_from][e];
                        const ptrdiff_t idx  = (gidx * vec_size + d) * from_stride;
                        in[sub_idx]          = from[idx];
                    } else {
                        in[sub_idx] = 0;
                    }

                    // Wait for in to be filled
                    __syncwarp();

                    int start_zi = from_zi * step_factor;
                    start_zi += is_odd(start_zi);  // Skip odd numbers

                    int start_yi = from_yi * step_factor;
                    start_yi += is_odd(start_yi);  // Skip odd numbers

                    int start_xi = from_xi * step_factor;
                    start_xi += is_odd(start_xi);  // Skip odd numbers

                    const int end_zi = MIN(to_nloops, start_zi + step_factor);
                    const int end_yi = MIN(to_nloops, start_yi + step_factor);
                    const int end_xi = MIN(to_nloops, start_xi + step_factor);

                    // sub-loop on even TO micro-elements
                    for (int to_zi = start_zi; to_zi < end_zi; to_zi += 2) {
                        for (int to_yi = start_yi; to_yi < end_yi; to_yi += 2) {
                            for (int to_xi = start_xi; to_xi < end_xi; to_xi += 2) {
                                const int off_to_zi = (to_zi + zi);
                                const int off_to_yi = (to_yi + yi);
                                const int off_to_xi = (to_xi + xi);

                                const To x = (off_to_xi - from_xi * step_factor) * between_h;
                                const To y = (off_to_yi - from_yi * step_factor) * between_h;
                                const To z = (off_to_zi - from_zi * step_factor) * between_h;

                                assert(x >= 0);
                                assert(x <= 1);
                                assert(y >= 0);
                                assert(y <= 1);
                                assert(z >= 0);
                                assert(z <= 1);

                                // This requires 64 bytes on the stack frame
                                // Cartesian order
                                To f[8] = {// Bottom
                                           (1 - x) * (1 - y) * (1 - z),
                                           x * (1 - y) * (1 - z),
                                           (1 - x) * y * (1 - z),
                                           x * y * (1 - z),
                                           // Top
                                           (1 - x) * (1 - y) * z,
                                           x * (1 - y) * z,
                                           (1 - x) * y * z,
                                           x * y * z};

#ifndef NDEBUG
                                To pou = 0;
                                for (int i = 0; i < 8; i++) {
                                    pou += f[i];
                                }

                                assert(fabs(1 - pou) < 1e-8);
#endif

                                To out = 0;
                                for (int v = 0; v < 8; v++) {
                                    const int round_robin = ROUND_ROBIN(v, sub_idx);
                                    // There should be no bank conflicts due to round robin
                                    out += f[round_robin] * in[round_robin];
                                }

                                // Check if not ghost nodes for scatter assign
                                const bool to_exists =
                                        e < nelements && off_to_zi <= to_level && off_to_yi <= to_level && off_to_xi <= to_level;

                                if (to_exists) {
                                    // set the value to the output
                                    const int idx_to = cu_sshex8_lidx(to_level * to_level_stride,
                                                                      off_to_xi * to_level_stride,
                                                                      off_to_yi * to_level_stride,
                                                                      off_to_zi * to_level_stride);

                                    to[(to_elements[idx_to][e] * vec_size + d) * to_stride] = out;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#else  // PROLONGATE_IN_KERNEL_BASIS_VERSION

template <typename From, typename To>
__global__ void cu_sshex8_prolongate_kernel(const ptrdiff_t                 nelements,
                                            const int                       from_level,
                                            const int                       from_level_stride,
                                            idx_t **const SFEM_RESTRICT     from_elements,
                                            const int                       to_level,
                                            const int                       to_level_stride,
                                            idx_t **const SFEM_RESTRICT     to_elements,
                                            const To *const SFEM_RESTRICT   S,
                                            const int                       vec_size,
                                            const enum RealType             from_type,
                                            const ptrdiff_t                 from_stride,
                                            const From *const SFEM_RESTRICT from,
                                            const enum RealType             to_type,
                                            const ptrdiff_t                 to_stride,
                                            To *const SFEM_RESTRICT         to) {
    static_assert(TILE_SIZE == 8, "This only works with tile size 8!");

    // Uunsigned char necessary for multiple instantiations
    extern __shared__ unsigned char cu_buff[];

    const int step_factor = to_level / from_level;
    const int to_npoints  = to_level + 1;

    // Tile number in group
    const int tile    = threadIdx.x >> 3;   // same as threadIdx.x / 8
    const int n_tiles = blockDim.x >> 3;    // same as blockDim.x / 8
    const int sub_idx = threadIdx.x & 0x7;  // same as threadIdx.x % 8

    // Potential bug ??
    From *in = (From *)&cu_buff[tile * TILE_SIZE * sizeof(From)];

    // hex8 idx
    const int xi = sub_idx & 0x1;         // equivalent to sub_idx % 2
    const int yi = (sub_idx >> 1) & 0x1;  // equivalent to (sub_idx / 2) % 2
    const int zi = (sub_idx >> 2);        // equivalent to sub_idx / 4
    assert(n_tiles * TILE_SIZE == blockDim.x);

    // 1 macro element per tile
    const ptrdiff_t e = blockIdx.x * n_tiles + tile;

    const int to_even     = is_even(to_level);
    const int from_nloops = from_level + to_even;
    const int to_nloops   = to_level + to_even;

    // Vector loop
    for (int d = 0; d < vec_size; d++) {
        // loop on all FROM micro elements
        for (int from_zi = 0; from_zi < from_nloops; from_zi++) {
            for (int from_yi = 0; from_yi < from_nloops; from_yi++) {
                for (int from_xi = 0; from_xi < from_nloops; from_xi++) {
                    const int off_from_zi = (from_zi + zi);
                    const int off_from_yi = (from_yi + yi);
                    const int off_from_xi = (from_xi + xi);

                    const bool from_exists =
                            e < nelements && off_from_zi <= from_level && off_from_yi <= from_level && off_from_xi <= from_level;

                    // Wait for shared memory transactions to be finished
                    __syncwarp();

                    // Gather
                    if (from_exists) {
                        const int idx_from = cu_sshex8_lidx(from_level * from_level_stride,
                                                            off_from_zi * from_level_stride,
                                                            off_from_yi * from_level_stride,
                                                            off_from_xi * from_level_stride);

                        const ptrdiff_t     gidx = from_elements[idx_from][e];
                        const ptrdiff_t idx  = (gidx * vec_size + d) * from_stride;
                        in[sub_idx]          = from[idx];
                    } else {
                        in[sub_idx] = 0;
                    }

                    // Wait for in to be filled
                    __syncwarp();

                    int start_zi = from_zi * step_factor;
                    start_zi += is_odd(start_zi);  // Skip odd numbers

                    int start_yi = from_yi * step_factor;
                    start_yi += is_odd(start_yi);  // Skip odd numbers

                    int start_xi = from_xi * step_factor;
                    start_xi += is_odd(start_xi);  // Skip odd numbers

                    const int end_zi = MIN(to_nloops, start_zi + step_factor);
                    const int end_yi = MIN(to_nloops, start_yi + step_factor);
                    const int end_xi = MIN(to_nloops, start_xi + step_factor);

                    // sub-loop on even TO micro-elements
                    for (int to_zi = start_zi; to_zi < end_zi; to_zi += 2) {
                        for (int to_yi = start_yi; to_yi < end_yi; to_yi += 2) {
                            for (int to_xi = start_xi; to_xi < end_xi; to_xi += 2) {
                                // Tile-level parallelism due to xi, yi, zi
                                const int off_to_zi = (to_zi + zi);
                                const int off_to_yi = (to_yi + yi);
                                const int off_to_xi = (to_xi + xi);

                                const From *const Sx = &S[(off_to_xi - from_xi)];
                                const From *const Sy = &S[(off_to_yi - from_yi)];
                                const From *const Sz = &S[(off_to_zi - from_zi)];

                                To out = 0;
                                // for (int vz = 0; vz < 2; vz++) {
                                //     const int rrvz = ROUND_ROBIN_2(vz, zi);
                                //     for (int vy = 0; vy < 2; vy++) {
                                //         const int rrvy = ROUND_ROBIN_2(vy, yi);
                                //         for (int vx = 0; vx < 2; vx++) {
                                //             const int rrvx = ROUND_ROBIN_2(vx, xi);
                                //             const int idx  = rrvz * 4 + rrvy * 2 + rrvx;
                                //             out += Sx[rrvx * to_npoints] * Sy[rrvy * to_npoints] * Sz[rrvz * to_npoints] *
                                //                    in[idx];
                                //         }
                                //     }
                                // }

                                for (int vz = 0; vz < 2; vz++) {
                                    for (int vy = 0; vy < 2; vy++) {
                                        for (int vx = 0; vx < 2; vx++) {
                                            const int idx = vz * 4 + vy * 2 + vx;
                                            out += Sx[vx * to_npoints] * Sy[vy * to_npoints] * Sz[vz * to_npoints] * in[idx];
                                        }
                                    }
                                }

                                // Check if not ghost nodes for scatter assign
                                const bool to_exists =
                                        e < nelements && off_to_zi <= to_level && off_to_yi <= to_level && off_to_xi <= to_level;

                                if (to_exists) {
                                    // set the value to the output
                                    const int idx_to = cu_sshex8_lidx(to_level * to_level_stride,
                                                                      off_to_xi * to_level_stride,
                                                                      off_to_yi * to_level_stride,
                                                                      off_to_zi * to_level_stride);

                                    to[(to_elements[idx_to][e] * vec_size + d) * to_stride] = out;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif  // PROLONGATE_IN_KERNEL_BASIS_VERSION

template <typename From, typename To>
int cu_sshex8_prolongate_tpl(const ptrdiff_t                 nelements,
                             const int                       from_level,
                             const int                       from_level_stride,
                             idx_t **const SFEM_RESTRICT     from_elements,
                             const int                       to_level,
                             const int                       to_level_stride,
                             idx_t **const SFEM_RESTRICT     to_elements,
                             const int                       vec_size,
                             const enum RealType             from_type,
                             const ptrdiff_t                 from_stride,
                             const From *const SFEM_RESTRICT from,
                             const enum RealType             to_type,
                             const ptrdiff_t                 to_stride,
                             To *const SFEM_RESTRICT         to,
                             void                           *stream) {
    SFEM_DEBUG_SYNCHRONIZE();

    const int block_size      = 128;
    ptrdiff_t n_blocks        = MAX(ptrdiff_t(1), (nelements + block_size / TILE_SIZE - 1) / (block_size / TILE_SIZE));
    size_t    shared_mem_size = block_size * sizeof(From);

#ifndef PROLONGATE_IN_KERNEL_BASIS_VERSION
    ShapeInterpolation<To> S(to_level / from_level);
#endif

    if (stream) {
        cudaStream_t s = *static_cast<cudaStream_t *>(stream);

        cu_sshex8_prolongate_kernel<From, To><<<n_blocks, block_size, shared_mem_size, s>>>(nelements,
                                                                                            from_level,
                                                                                            from_level_stride,
                                                                                            from_elements,
                                                                                            to_level,
                                                                                            to_level_stride,
                                                                                            to_elements,
#ifndef PROLONGATE_IN_KERNEL_BASIS_VERSION
                                                                                            S.data,
#endif
                                                                                            vec_size,
                                                                                            from_type,
                                                                                            from_stride,
                                                                                            from,
                                                                                            to_type,
                                                                                            to_stride,
                                                                                            to);
    } else {
        cu_sshex8_prolongate_kernel<From, To><<<n_blocks, block_size, shared_mem_size>>>(nelements,
                                                                                         from_level,
                                                                                         from_level_stride,
                                                                                         from_elements,
                                                                                         to_level,
                                                                                         to_level_stride,
                                                                                         to_elements,
#ifndef PROLONGATE_IN_KERNEL_BASIS_VERSION
                                                                                         S.data,
#endif
                                                                                         vec_size,
                                                                                         from_type,
                                                                                         from_stride,
                                                                                         from,
                                                                                         to_type,
                                                                                         to_stride,
                                                                                         to);
    }

    SFEM_DEBUG_SYNCHRONIZE();
    return SFEM_SUCCESS;
}

extern int cu_sshex8_prolongate(const ptrdiff_t                 nelements,
                                const int                       from_level,
                                const int                       from_level_stride,
                                idx_t **const SFEM_RESTRICT     from_elements,
                                const int                       to_level,
                                const int                       to_level_stride,
                                idx_t **const SFEM_RESTRICT     to_elements,
                                const int                       vec_size,
                                const enum RealType             from_type,
                                const ptrdiff_t                 from_stride,
                                const void *const SFEM_RESTRICT from,
                                const enum RealType             to_type,
                                const ptrdiff_t                 to_stride,
                                void *const SFEM_RESTRICT       to,
                                void                           *stream) {
    assert(from_type == to_type && "TODO mixed types!");
    if (from_type != to_type) {
        return SFEM_FAILURE;
    }

    switch (from_type) {
        case SFEM_REAL_DEFAULT: {
            return cu_sshex8_prolongate_tpl<real_t, real_t>(nelements,
                                                            from_level,
                                                            from_level_stride,
                                                            from_elements,
                                                            to_level,
                                                            to_level_stride,
                                                            to_elements,
                                                            vec_size,
                                                            from_type,
                                                            from_stride,
                                                            (real_t *)from,
                                                            to_type,
                                                            to_stride,
                                                            (real_t *)to,
                                                            stream);
        }
        case SFEM_FLOAT32: {
            return cu_sshex8_prolongate_tpl<float, float>(nelements,
                                                          from_level,
                                                          from_level_stride,
                                                          from_elements,
                                                          to_level,
                                                          to_level_stride,
                                                          to_elements,
                                                          vec_size,
                                                          from_type,
                                                          from_stride,
                                                          (float *)from,
                                                          to_type,
                                                          to_stride,
                                                          (float *)to,
                                                          stream);
        }
        case SFEM_FLOAT64: {
            return cu_sshex8_prolongate_tpl<double, double>(nelements,
                                                            from_level,
                                                            from_level_stride,
                                                            from_elements,
                                                            to_level,
                                                            to_level_stride,
                                                            to_elements,
                                                            vec_size,
                                                            from_type,
                                                            from_stride,
                                                            (double *)from,
                                                            to_type,
                                                            to_stride,
                                                            (double *)to,
                                                            stream);
        }

        default: {
            SFEM_ERROR(
                    "[Error]  cu_sshex8_prolongate: not implemented for type "
                    "%s "
                    "(code %d)\n",
                    real_type_to_string(from_type),
                    from_type);
            return SFEM_FAILURE;
        }
    }
}
