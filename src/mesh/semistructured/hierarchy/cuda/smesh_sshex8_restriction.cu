#include "smesh_sshex8_restriction.cuh"

#include "smesh_cuda_base.cuh"

#include "smesh_macros.hpp"
#if defined(__HIP_PLATFORM_SPIRV__)
#include "smesh_context.hpp"
#endif
#include "smesh_sshex8.hpp"
#include "smesh_types.hpp"

#include <cstdio>
#include <mutex>
#include <type_traits>
#include <vector>

namespace smesh {

    static const int TILE_SIZE = 8;
#define ROUND_ROBIN(val, shift) ((val + shift) & (TILE_SIZE - 1))
#define ROUND_ROBIN_2(val, shift) ((val + shift) & (2 - 1))

    static inline __device__ int is_even(const int v) { return !(v & 1); }
    static inline __device__ int is_odd(const int v) { return (v & 1); }

    // Each file's own (ssquad4 has one with the same names): an anonymous namespace keeps the two apart.
    namespace {

    template <typename T>
    class ShapeInterpolation {
    public:
        T     *data{nullptr};
        size_t nodes{0};
        int    stride{0};
        int    steps{0};
        int    padding{0};

        int init(const int steps_, const int padding_ = 0) {
            steps   = steps_;
            padding = padding_;
            nodes   = (steps + 1);
            stride  = nodes + padding;
            std::vector<T> S_host(2 * stride, 0);
            f64            h = 1. / steps_;
            for (int i = 0; i < nodes; i++) {
                S_host[0 * stride + i] = (1 - h * i);
                S_host[1 * stride + i] = h * i;
            }

            auto nbytes = S_host.size() * sizeof(T);

            SMESH_CUDA_CHECK(cudaMalloc((void **)&data, nbytes));
            SMESH_CUDA_CHECK(cudaMemcpy(data, S_host.data(), nbytes, cudaMemcpyHostToDevice));

#if 0
        for (int i = 0; i < stride; i++) {
            printf("%g ", (f64)S_host[0 * stride + i]);
        }
        printf("\n");

        for (int i = 0; i < stride; i++) {
            printf("%g ", (f64)S_host[1 * stride + i]);
        }
        printf("\n");
#endif
            return SMESH_SUCCESS;
        }

        bool matches(const int steps_, const int padding_) const { return data && steps == steps_ && padding == padding_; }

        ~ShapeInterpolation() { destroy(); }

        void destroy() {
            if (data) {
                cudaFree(data);
                data = nullptr;
            }
        }
    };

    template <typename T>
    class ShapeInterpolationPool {
    public:
        static ShapeInterpolation<T> *get(const int steps, const int padding) {
            std::lock_guard<std::mutex> lock(mutex());
#if defined(__HIP_PLATFORM_SPIRV__)
            // Released at the context's teardown under chipStar (smesh_context.hpp), not at exit.
            static const bool released_at_teardown = (at_context_teardown(&ShapeInterpolationPool<T>::destroy), true);
            (void)released_at_teardown;
#endif

            ShapeInterpolation<T> *empty_slot = nullptr;
            for (int i = 0; i < capacity(); ++i) {
                auto &entry = entries()[i];
                if (entry.matches(steps, padding)) {
                    return &entry;
                }

                if (!entry.data && !empty_slot) {
                    empty_slot = &entry;
                }
            }

            if (!empty_slot) {
                SMESH_ERROR("ShapeInterpolationPool capacity exceeded for steps=%d padding=%d\n", steps, padding);
                return nullptr;
            }

            if (empty_slot->init(steps, padding) != SMESH_SUCCESS) {
                return nullptr;
            }

            return empty_slot;
        }

        static void destroy() {
            std::lock_guard<std::mutex> lock(mutex());
            for (int i = 0; i < capacity(); ++i) {
                entries()[i].destroy();
            }
        }

    private:
        static constexpr int capacity() { return 64; }

        static ShapeInterpolation<T> *entries() {
            static ShapeInterpolation<T> pool[capacity()];
            return pool;
        }

        static std::mutex &mutex() {
            static std::mutex m;
            return m;
        }
    };

    }  // namespace

    // PROLONGATION

    template <typename From, typename To, typename idx_t>
    __global__ void cu_sshex8_hierarchical_prolongation_kernel(const int                                               level,
                                                               const ptrdiff_t                                         nelements,
                                                               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                                               const int                                               vec_size,
                                                               const ptrdiff_t                  from_stride,
                                                               const From *const SMESH_RESTRICT from,
                                                               const ptrdiff_t                  to_stride,
                                                               To *const SMESH_RESTRICT         to) {
        const int corners[8] = {// Bottom
                                sshex8_lidx(level, 0, 0, 0),
                                sshex8_lidx(level, level, 0, 0),
                                sshex8_lidx(level, level, level, 0),
                                sshex8_lidx(level, 0, level, 0),
                                // Top
                                sshex8_lidx(level, 0, 0, level),
                                sshex8_lidx(level, level, 0, level),
                                sshex8_lidx(level, level, level, level),
                                sshex8_lidx(level, 0, level, level)};

        const scalar_t h = 1. / level;

        for (ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x; e < nelements; e += blockDim.x * gridDim.x) {
            for (int zi = 0; zi < level + 1; zi++) {
                for (int yi = 0; yi < level + 1; yi++) {
                    for (int xi = 0; xi < level + 1; xi++) {
                        idx_t idx = elements[sshex8_lidx(level, xi, yi, zi)][e];

                        const scalar_t x = xi * h;
                        const scalar_t y = yi * h;
                        const scalar_t z = zi * h;

                        // Evaluate Hex8 basis functions at x, y, z
                        const scalar_t xm = (1 - x);
                        const scalar_t ym = (1 - y);
                        const scalar_t zm = (1 - z);

                        scalar_t f[8];
                        f[0] = xm * ym * zm;  // (0, 0, 0)
                        f[1] = x * ym * zm;   // (1, 0, 0)
                        f[2] = x * y * zm;    // (1, 1, 0)
                        f[3] = xm * y * zm;   // (0, 1, 0)
                        f[4] = xm * ym * z;   // (0, 0, 1)
                        f[5] = x * ym * z;    // (1, 0, 1)
                        f[6] = x * y * z;     // (1, 1, 1)
                        f[7] = xm * y * z;    // (0, 1, 1)

                        for (int d = 0; d < vec_size; d++) {
                            scalar_t val = 0;

                            for (int v = 0; v < 8; v++) {
                                const ptrdiff_t global_from_idx = (elements[corners[v]][e] * vec_size + d) * from_stride;

                                SMESH_ASSERT(from[global_from_idx] == from[global_from_idx]);
                                SMESH_ASSERT(f[v] == f[v]);

                                val += f[v] * from[global_from_idx];
                            }

                            SMESH_ASSERT(val == val);

                            const ptrdiff_t global_to_idx = (idx * vec_size + d) * to_stride;
                            to[global_to_idx]             = val;
                        }
                    }
                }
            }
        }
    }

    template <typename From, typename To, typename idx_t>
    static int cu_sshex8_hierarchical_prolongation_tpl(const int                                               level,
                                                       const ptrdiff_t                                         nelements,
                                                       const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                                       const int                                               vec_size,
                                                       const ptrdiff_t                                         from_stride,
                                                       const From *const SMESH_RESTRICT                        from,
                                                       const ptrdiff_t                                         to_stride,
                                                       To *const SMESH_RESTRICT                                to,
                                                       void                                                   *stream) {
        SMESH_DEBUG_SYNCHRONIZE();

        // Hand tuned
        int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
        {
            int min_grid_size;
            cudaOccupancyMaxPotentialBlockSize(
                    &min_grid_size, &block_size, cu_sshex8_hierarchical_prolongation_kernel<From, To>, 0, 0);
        }
#endif  // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

        ptrdiff_t n_blocks = MAX(ptrdiff_t(1), (nelements + block_size - 1) / block_size);

        cu_sshex8_hierarchical_prolongation_kernel<From, To><<<n_blocks, block_size, 0, (cudaStream_t)stream>>>(
                level, nelements, elements, vec_size, from_stride, from, to_stride, to);

        SMESH_DEBUG_SYNCHRONIZE();
        return SMESH_SUCCESS;
    }

    int cu_sshex8_hierarchical_prolongation(const int                        level,
                                            const ptrdiff_t                  nelements,
                                            idx_t **const SMESH_RESTRICT     elements,
                                            const int                        vec_size,
                                            const enum PrimitiveType         from_type,
                                            const ptrdiff_t                  from_stride,
                                            const void *const SMESH_RESTRICT from,
                                            const enum PrimitiveType         to_type,
                                            const ptrdiff_t                  to_stride,
                                            void *const SMESH_RESTRICT       to,
                                            void                            *stream) {
        SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
        if (from_type != to_type) {
            return SMESH_FAILURE;
        }

        switch (from_type) {
            case SMESH_DEFAULT: {
                return cu_sshex8_hierarchical_prolongation_tpl(
                        level, nelements, elements, vec_size, from_stride, (real_t *)from, to_stride, (real_t *)to, stream);
            }
            case SMESH_FLOAT32: {
                return cu_sshex8_hierarchical_prolongation_tpl(
                        level, nelements, elements, vec_size, from_stride, (f32 *)from, to_stride, (f32 *)to, stream);
            }
            case SMESH_FLOAT64: {
                return cu_sshex8_hierarchical_prolongation_tpl(
                        level, nelements, elements, vec_size, from_stride, (f64 *)from, to_stride, (f64 *)to, stream);
            }
            default: {
                SMESH_ERROR(
                        "[Error]  cu_sshex8_prolongation_tpl: not implemented for type %s "
                        "(code %d)\n",
                        to_string(from_type),
                        from_type);
                return SMESH_FAILURE;
            }
        }
    }

    template <typename From, typename To, typename idx_t>
    __global__ void cu_sshex8_hierarchical_restriction_kernel(const int                                               level,
                                                              const ptrdiff_t                                         nelements,
                                                              const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                                              const u16 *const SMESH_RESTRICT                         e2n_count,
                                                              const int                                               vec_size,
                                                              const ptrdiff_t                                         from_stride,
                                                              const From *const SMESH_RESTRICT                        from,
                                                              const ptrdiff_t                                         to_stride,
                                                              To *const SMESH_RESTRICT                                to) {
        const int corners[8] = {// Bottom
                                sshex8_lidx(level, 0, 0, 0),
                                sshex8_lidx(level, level, 0, 0),
                                sshex8_lidx(level, level, level, 0),
                                sshex8_lidx(level, 0, level, 0),
                                // Top
                                sshex8_lidx(level, 0, 0, level),
                                sshex8_lidx(level, level, 0, level),
                                sshex8_lidx(level, level, level, level),
                                sshex8_lidx(level, 0, level, level)};

        const scalar_t h = 1. / level;
        scalar_t       acc[8];

        for (ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x; e < nelements; e += blockDim.x * gridDim.x) {
            for (int d = 0; d < vec_size; d++) {
                for (int i = 0; i < 8; i++) {
                    acc[i] = 0;
                }

                for (int zi = 0; zi < level + 1; zi++) {
                    for (int yi = 0; yi < level + 1; yi++) {
                        for (int xi = 0; xi < level + 1; xi++) {
                            const int       lidx = sshex8_lidx(level, xi, yi, zi);
                            const ptrdiff_t idx  = elements[lidx][e];

                            const scalar_t x = xi * h;
                            const scalar_t y = yi * h;
                            const scalar_t z = zi * h;

                            // Evaluate Hex8 basis functions at x, y, z
                            const scalar_t xm = (1 - x);
                            const scalar_t ym = (1 - y);
                            const scalar_t zm = (1 - z);

                            scalar_t f[8];
                            f[0] = xm * ym * zm;  // (0, 0, 0)
                            f[1] = x * ym * zm;   // (1, 0, 0)
                            f[2] = x * y * zm;    // (1, 1, 0)
                            f[3] = xm * y * zm;   // (0, 1, 0)
                            f[4] = xm * ym * z;   // (0, 0, 1)
                            f[5] = x * ym * z;    // (1, 0, 1)
                            f[6] = x * y * z;     // (1, 1, 1)
                            f[7] = xm * y * z;    // (0, 1, 1)

                            const ptrdiff_t global_from_idx = (idx * vec_size + d) * from_stride;
                            const scalar_t  val             = from[global_from_idx] / e2n_count[idx];

                            SMESH_ASSERT(from[global_from_idx] == from[global_from_idx]);
                            SMESH_ASSERT(e2n_count[idx] > 0);
                            SMESH_ASSERT(val == val);

                            for (int i = 0; i < 8; i++) {
                                acc[i] += f[i] * val;
                            }
                        }
                    }
                }

                for (int v = 0; v < 8; v++) {
                    const ptrdiff_t global_to_idx = (elements[corners[v]][e] * vec_size + d) * to_stride;
                    atomicAdd(&to[global_to_idx], acc[v]);
                }
            }
        }
    }

    template <typename From, typename To, typename idx_t>
    static int cu_sshex8_hierarchical_restriction_tpl(const int                                               level,
                                                      const ptrdiff_t                                         nelements,
                                                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                                      const u16 *const SMESH_RESTRICT  element_to_node_incidence_count,
                                                      const int                        vec_size,
                                                      const ptrdiff_t                  from_stride,
                                                      const From *const SMESH_RESTRICT from,
                                                      const ptrdiff_t                  to_stride,
                                                      To *const SMESH_RESTRICT         to,
                                                      void                            *stream) {
        SMESH_DEBUG_SYNCHRONIZE();

        // Hand tuned
        int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
        {
            int min_grid_size;
            cudaOccupancyMaxPotentialBlockSize(
                    &min_grid_size, &block_size, cu_sshex8_hierarchical_restriction_kernel<From, To, idx_t>, 0, 0);
        }
#endif  // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

        ptrdiff_t n_blocks = MAX(ptrdiff_t(1), (nelements + block_size - 1) / block_size);

        cu_sshex8_hierarchical_restriction_kernel<From, To, idx_t><<<n_blocks, block_size, 0, (cudaStream_t)stream>>>(
                level, nelements, elements, element_to_node_incidence_count, vec_size, from_stride, from, to_stride, to);

        SMESH_DEBUG_SYNCHRONIZE();
        return SMESH_SUCCESS;
    }

    extern int cu_sshex8_hierarchical_restriction(const int                                               level,
                                                  const ptrdiff_t                                         nelements,
                                                  const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                                  const u16 *const SMESH_RESTRICT  element_to_node_incidence_count,
                                                  const int                        vec_size,
                                                  const enum PrimitiveType         from_type,
                                                  const ptrdiff_t                  from_stride,
                                                  const void *const SMESH_RESTRICT from,
                                                  const enum PrimitiveType         to_type,
                                                  const ptrdiff_t                  to_stride,
                                                  void *const SMESH_RESTRICT       to,
                                                  void                            *stream) {
        SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
        if (from_type != to_type) {
            return SMESH_FAILURE;
        }

        switch (from_type) {
            case SMESH_DEFAULT: {
                return cu_sshex8_hierarchical_restriction_tpl(level,
                                                              nelements,
                                                              elements,
                                                              element_to_node_incidence_count,
                                                              vec_size,
                                                              from_stride,
                                                              (real_t *)from,
                                                              to_stride,
                                                              (real_t *)to,
                                                              stream);
            }
            case SMESH_FLOAT32: {
                return cu_sshex8_hierarchical_restriction_tpl(level,
                                                              nelements,
                                                              elements,
                                                              element_to_node_incidence_count,
                                                              vec_size,
                                                              from_stride,
                                                              (f32 *)from,
                                                              to_stride,
                                                              (f32 *)to,
                                                              stream);
            }
            case SMESH_FLOAT64: {
                return cu_sshex8_hierarchical_restriction_tpl(level,
                                                              nelements,
                                                              elements,
                                                              element_to_node_incidence_count,
                                                              vec_size,
                                                              from_stride,
                                                              (f64 *)from,
                                                              to_stride,
                                                              (f64 *)to,
                                                              stream);
            }
            default: {
                SMESH_ERROR(
                        "[Error]  cu_sshex8_prolongation_tpl: not implemented for type "
                        "%s "
                        "(code %d)\n",
                        to_string(from_type),
                        from_type);
                return SMESH_FAILURE;
            }
        }
    }

    template <typename From, typename To>
    __global__ void sshex8_hierarchical_restriction_kernel(const ptrdiff_t nelements,
                                                           const int       from_level,
                                                           const int       from_level_stride,
                                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                                                           const u16 *const SMESH_RESTRICT from_element_to_node_incidence_count,
                                                           const int                       to_level,
                                                           const int                       to_level_stride,
                                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                                                           const To *const SMESH_RESTRICT                          S,
                                                           const int                                               vec_size,
                                                           const enum PrimitiveType                                from_type,
                                                           const ptrdiff_t                                         from_stride,
                                                           const From *const SMESH_RESTRICT                        from,
                                                           const enum PrimitiveType                                to_type,
                                                           const ptrdiff_t                                         to_stride,
                                                           To *const SMESH_RESTRICT                                to) {
        static_assert(TILE_SIZE == 8,
                      "This only works with tile size 8 because the implementation assumes a "
                      "fixed tile size for shared memory "
                      "layout and indexing.");

        // Unsigned char necessary for multiple template instantiations of this kernel
        extern __shared__ unsigned char cu_buff[];

        const int step_factor = from_level / to_level;

        // Tile number in group
        const int tile    = threadIdx.x >> 3;   // same as threadIdx.x / 8
        const int n_tiles = blockDim.x >> 3;    // same as blockDim.x / 8
        const int sub_idx = threadIdx.x & 0x7;  // same as threadIdx.x % 8

        From *in = (From *)&cu_buff[tile * TILE_SIZE * sizeof(From)];

        // hex8 idx
        const int xi = sub_idx & 0x1;         // equivalent to sub_idx % 2
        const int yi = (sub_idx >> 1) & 0x1;  // equivalent to (sub_idx / 2) % 2
        const int zi = (sub_idx >> 2);        // equivalent to sub_idx / 4
        SMESH_ASSERT(n_tiles * TILE_SIZE == blockDim.x);

        // // 1 macro element per tile
        const ptrdiff_t e         = blockIdx.x * n_tiles + tile;
        const int       from_even = is_even(from_level);
        const int       to_nloops = to_level + from_even;

        // Add padding (Make sure that S has the correct padding)
        const int S_stride = step_factor + 1 + from_even;

        // // Vector loop
        for (int d = 0; d < vec_size; d++) {
            // loop on all TO micro elements
            for (int to_zi = 0; to_zi < to_nloops; to_zi++) {
                for (int to_yi = 0; to_yi < to_nloops; to_yi++) {
                    for (int to_xi = 0; to_xi < to_nloops; to_xi++) {
                        // Attention: parallelism in tile using xi, yi, zi
                        To acc = 0;

                        const int z_start = to_zi * step_factor;
                        const int y_start = to_yi * step_factor;
                        const int x_start = to_xi * step_factor;

                        const int z_odd = is_odd(z_start);
                        const int y_odd = is_odd(y_start);
                        const int x_odd = is_odd(x_start);

                        // Figure out local index range
                        const int z_end = (z_start == to_level ? 1 : step_factor) + z_odd;
                        const int y_end = (y_start == to_level ? 1 : step_factor) + y_odd;
                        const int x_end = (x_start == to_level ? 1 : step_factor) + x_odd;

                        for (int from_zi = z_odd; from_zi < z_end; from_zi += 2) {
                            for (int from_yi = y_odd; from_yi < y_end; from_yi += 2) {
                                for (int from_xi = x_odd; from_xi < x_end; from_xi += 2) {
                                    // Parallel read from global mem on 8 fine nodes
                                    {
                                        const int zz = z_start + from_zi;
                                        const int yy = y_start + from_yi;
                                        const int xx = x_start + from_xi;

                                        const int off_from_zi = (zz + zi);
                                        const int off_from_yi = (yy + yi);
                                        const int off_from_xi = (xx + xi);

                                        const bool from_exists = e < nelements &&              // Check element exists
                                                                 off_from_zi <= from_level &&  //
                                                                 off_from_yi <= from_level &&  //
                                                                 off_from_xi <= from_level;

                                        __syncwarp();
                                        if (from_exists) {
                                            const int idx_from = sshex8_lidx(from_level * from_level_stride,
                                                                             off_from_xi * from_level_stride,
                                                                             off_from_yi * from_level_stride,
                                                                             off_from_zi * from_level_stride);

                                            const idx_t     gidx = from_elements[idx_from][e];
                                            const ptrdiff_t idx  = (gidx * vec_size + d) * from_stride;
                                            in[sub_idx]          = from[idx] / from_element_to_node_incidence_count[gidx];
                                        } else {
                                            in[sub_idx] = 0;
                                        }
                                        __syncwarp();
                                    }

                                    const To *const Sz = &S[zi * S_stride + from_zi];
                                    const To *const Sy = &S[yi * S_stride + from_yi];
                                    const To *const Sx = &S[xi * S_stride + from_xi];

                                    for (int dz = 0; dz < 2; dz++) {
                                        const int rrdz = ROUND_ROBIN_2(dz, zi);
                                        for (int dy = 0; dy < 2; dy++) {
                                            const int rrdy = ROUND_ROBIN_2(dy, yi);
                                            for (int dx = 0; dx < 2; dx++) {
                                                const int rrdx = ROUND_ROBIN_2(dx, xi);
                                                // No bank conflicts due to round robin for single
                                                // precision
                                                const To c = in[rrdz * 4 + rrdy * 2 + rrdx];
                                                const To f = Sx[rrdx] * Sy[rrdy] * Sz[rrdz];
                                                SMESH_ASSERT(f >= 0);
                                                SMESH_ASSERT(f <= 1);
                                                acc += c * f;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // Parallel accumulate on 8 coarse nodes
                        const int off_to_zi = to_zi + zi;
                        const int off_to_yi = to_yi + yi;
                        const int off_to_xi = to_xi + xi;

                        const bool exists = e < nelements &&          // Check element exists
                                            off_to_zi <= to_level &&  //
                                            off_to_yi <= to_level &&  //
                                            off_to_xi <= to_level;
                        if (exists) {
                            const int idx_to = sshex8_lidx(to_level * to_level_stride,
                                                           off_to_xi * to_level_stride,
                                                           off_to_yi * to_level_stride,
                                                           off_to_zi * to_level_stride);

                            atomicAdd(&to[((ptrdiff_t)to_elements[idx_to][e] * vec_size + d) * to_stride], acc);
                        }
                    }
                }
            }
        }
    }

    // Node-centric restriction FROM_L -> FROM_L / STEP of EB consecutive macro elements per block:
    // the fine values (divided by their incidence count) are staged in shared memory one component at
    // a time, reduced in x, y, z with the hat weights, and each coarse value is added once. Level 8 ->
    // 4, 44^3 macro elements (132 M dofs), standalone: MI300A 2.37 ms against 8.97 for the tile-of-8
    // kernel above (one element per block of 512: 2.82 ms).
    constexpr int NODE_CENTRIC_BLOCK_SIZE = 1024;

    // Static shared memory of the kernel below: element indices, fine values, and the x-reduced values.
    template <typename T, int FROM_L, int STEP>
    constexpr size_t sshex8_restriction_node_centric_shared_bytes(const int EB) {
        return size_t(EB) * (FROM_L + 1) * (FROM_L + 1) *
               ((FROM_L + 1) * (sizeof(idx_t) + sizeof(T)) + (FROM_L / STEP + 1) * sizeof(T));
    }

    template <typename T, int FROM_L, int STEP, int EB>
    __global__ void __launch_bounds__(NODE_CENTRIC_BLOCK_SIZE) sshex8_restriction_node_centric_kernel(const ptrdiff_t                                         nelements,
                                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                                                           const u16 *const SMESH_RESTRICT                         count,
                                                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                                                           const int                                               vec_size,
                                                           const T *const SMESH_RESTRICT                           from,
                                                           T *const SMESH_RESTRICT                                 to) {
        constexpr int NF  = FROM_L + 1;
        constexpr int NC  = FROM_L / STEP + 1;
        constexpr int NF3 = NF * NF * NF;
        __shared__ idx_t gidx[EB * NF3];           // [f][e]
        __shared__ T     fine[EB * NF3];           // [e][z][y][x]; after the y pass [e][z][Y][X]
        __shared__ T     xred[EB * NF * NF * NC];  // [e][z][y][X]

        for (ptrdiff_t e0 = (ptrdiff_t)blockIdx.x * EB; e0 < nelements; e0 += (ptrdiff_t)gridDim.x * EB) {
            for (int i = threadIdx.x; i < EB * NF3; i += blockDim.x) {
                const int el = i % EB;
                gidx[i]      = e0 + el < nelements ? from_elements[i / EB][e0 + el] : -1;
            }
            __syncthreads();

            for (int d = 0; d < vec_size; d++) {
                for (int i = threadIdx.x; i < EB * NF3; i += blockDim.x) {
                    const idx_t g              = gidx[i];
                    fine[(i % EB) * NF3 + i / EB] = g < 0 ? T(0) : from[(ptrdiff_t)g * vec_size + d] / count[g];
                }
                __syncthreads();

                // Coarse X takes fine STEP X + k, |k| < STEP, with weight 1 - |k| / STEP
                for (int i = threadIdx.x; i < EB * NF * NF * NC; i += blockDim.x) {
                    const int      X   = i % NC;
                    const T *const row = &fine[(i / NC) * NF];
                    T              acc = 0;
#pragma unroll
                    for (int k = 1 - STEP; k < STEP; k++) {
                        const int x = STEP * X + k;
                        if (x >= 0 && x <= FROM_L) acc += (T(1) - T(k < 0 ? -k : k) / STEP) * row[x];
                    }
                    xred[i] = acc;
                }
                __syncthreads();

                for (int i = threadIdx.x; i < EB * NF * NC * NC; i += blockDim.x) {
                    const int X  = i % NC;
                    const int Y  = (i / NC) % NC;
                    const int ez = i / (NC * NC);
                    T         acc = 0;
#pragma unroll
                    for (int k = 1 - STEP; k < STEP; k++) {
                        const int y = STEP * Y + k;
                        if (y >= 0 && y <= FROM_L) acc += (T(1) - T(k < 0 ? -k : k) / STEP) * xred[(ez * NF + y) * NC + X];
                    }
                    fine[i] = acc;
                }
                __syncthreads();

                for (int i = threadIdx.x; i < EB * NC * NC * NC; i += blockDim.x) {
                    const int el   = i % EB;
                    const int node = i / EB;
                    if (e0 + el >= nelements) continue;
                    const int X   = node % NC;
                    const int Y   = (node / NC) % NC;
                    const int Z   = node / (NC * NC);
                    T         acc = 0;
#pragma unroll
                    for (int k = 1 - STEP; k < STEP; k++) {
                        const int z = STEP * Z + k;
                        if (z >= 0 && z <= FROM_L) acc += (T(1) - T(k < 0 ? -k : k) / STEP) * fine[((el * NF + z) * NC + Y) * NC + X];
                    }
                    atomicAdd(&to[(ptrdiff_t)to_elements[node][e0 + el] * vec_size + d], acc);
                }
                __syncthreads();
            }
        }
    }

    template <typename From, typename To, typename idx_t>
    int cu_sshex8_restrict_tpl(const ptrdiff_t                                         nelements,
                               const int                                               from_level,
                               const int                                               from_level_stride,
                               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                               const u16 *const SMESH_RESTRICT from_element_to_node_incidence_count,
                               const int                       to_level,
                               const int                       to_level_stride,
                               const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                               const int                                               vec_size,
                               const enum PrimitiveType                                from_type,
                               const ptrdiff_t                                         from_stride,
                               const From *const SMESH_RESTRICT                        from,
                               const enum PrimitiveType                                to_type,
                               const ptrdiff_t                                         to_stride,
                               To *const SMESH_RESTRICT                                to,
                               void                                                   *stream) {
        SMESH_DEBUG_SYNCHRONIZE();

        if constexpr (std::is_same_v<From, To>) {
            if (from_level == 8 && to_level == 4 && from_level_stride == 1 && to_level_stride == 1 && from_stride == 1 &&
                to_stride == 1) {
                // No blocks to launch on a rank without elements (a grid of 0 is an invalid configuration)
                if (nelements == 0) return SMESH_SUCCESS;

                // Four macro elements per block, three when 64-bit indices would exceed the 48 KiB of
                // static shared memory CUDA allows a block (f64 and i64: 59,616 bytes for four)
                constexpr size_t max_static_shared = 48 * 1024;
                constexpr int    EB =
                        sshex8_restriction_node_centric_shared_bytes<To, 8, 2>(4) <= max_static_shared ? 4 : 3;
                static_assert(sshex8_restriction_node_centric_shared_bytes<To, 8, 2>(EB) <= max_static_shared,
                              "sshex8_restriction_node_centric_kernel exceeds 48 KiB of static shared memory");

                sshex8_restriction_node_centric_kernel<To, 8, 2, EB>
                        <<<(nelements + EB - 1) / EB, NODE_CENTRIC_BLOCK_SIZE, 0, (cudaStream_t)stream>>>(
                                nelements, from_elements, from_element_to_node_incidence_count, to_elements, vec_size, from, to);
                SMESH_DEBUG_SYNCHRONIZE();
                return SMESH_SUCCESS;
            }
        }

        // Hand tuned
        int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
        {
            int min_grid_size;
            cudaOccupancyMaxPotentialBlockSize(
                    &min_grid_size, &block_size, sshex8_hierarchical_restriction_kernel<From, To>, 0, 0);
        }
#endif  // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

        ptrdiff_t n_blocks = MAX(ptrdiff_t(1), (nelements + block_size / TILE_SIZE - 1) / (block_size / TILE_SIZE));

        auto S = ShapeInterpolationPool<To>::get(from_level / to_level, from_level % 2 == 0);
        if (!S) {
            return SMESH_FAILURE;
        }

        size_t shared_mem_size = block_size * sizeof(From);

        sshex8_hierarchical_restriction_kernel<From, To>
                <<<n_blocks, block_size, shared_mem_size, (cudaStream_t)stream>>>(nelements,
                                                                                  from_level,
                                                                                  from_level_stride,
                                                                                  from_elements,
                                                                                  from_element_to_node_incidence_count,
                                                                                  to_level,
                                                                                  to_level_stride,
                                                                                  to_elements,
                                                                                  S->data,
                                                                                  vec_size,
                                                                                  from_type,
                                                                                  from_stride,
                                                                                  from,
                                                                                  to_type,
                                                                                  to_stride,
                                                                                  to);

        SMESH_DEBUG_SYNCHRONIZE();
        return SMESH_SUCCESS;
    }

    int cu_sshex8_restrict(const ptrdiff_t                                         nelements,
                           const int                                               from_level,
                           const int                                               from_level_stride,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
                           const u16 *const SMESH_RESTRICT  from_element_to_node_incidence_count,
                           const int                        to_level,
                           const int                        to_level_stride,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
                           const int                        vec_size,
                           const enum PrimitiveType         from_type,
                           const ptrdiff_t                  from_stride,
                           const void *const SMESH_RESTRICT from,
                           const enum PrimitiveType         to_type,
                           const ptrdiff_t                  to_stride,
                           void *const SMESH_RESTRICT       to,
                           void                            *stream) {
        SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
        if (from_type != to_type) {
            SMESH_ERROR("cu_sshex8_restrict requires matching from/to types (%s != %s)\n", to_string(from_type), to_string(to_type));
            return SMESH_FAILURE;
        }

        if (from_level <= 0 || to_level <= 0 || from_level % to_level != 0) {
            SMESH_ERROR("cu_sshex8_restrict requires positive nested levels with from_level %% to_level == 0 (from=%d, to=%d)\n",
                        from_level,
                        to_level);
            return SMESH_FAILURE;
        }

        if (from_level_stride <= 0 || to_level_stride <= 0 || vec_size <= 0) {
            SMESH_ERROR("cu_sshex8_restrict requires positive strides and vec_size (from_level_stride=%d, to_level_stride=%d, vec_size=%d)\n",
                        from_level_stride,
                        to_level_stride,
                        vec_size);
            return SMESH_FAILURE;
        }

        switch (from_type) {
            case SMESH_DEFAULT: {
                return cu_sshex8_restrict_tpl<real_t, real_t, idx_t>(nelements,
                                                                      from_level,
                                                                      from_level_stride,
                                                                      from_elements,
                                                                      from_element_to_node_incidence_count,
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
            case SMESH_FLOAT32: {
                return cu_sshex8_restrict_tpl<f32, f32, idx_t>(nelements,
                                                                from_level,
                                                                from_level_stride,
                                                                from_elements,
                                                                from_element_to_node_incidence_count,
                                                                to_level,
                                                                to_level_stride,
                                                                to_elements,
                                                                vec_size,
                                                                from_type,
                                                                from_stride,
                                                                (f32 *)from,
                                                                to_type,
                                                                to_stride,
                                                                (f32 *)to,
                                                                stream);
            }
            case SMESH_FLOAT64: {
                return cu_sshex8_restrict_tpl<f64, f64, idx_t>(nelements,
                                                                from_level,
                                                                from_level_stride,
                                                                from_elements,
                                                                from_element_to_node_incidence_count,
                                                                to_level,
                                                                to_level_stride,
                                                                to_elements,
                                                                vec_size,
                                                                from_type,
                                                                from_stride,
                                                                (f64 *)from,
                                                                to_type,
                                                                to_stride,
                                                                (f64 *)to,
                                                                stream);
            }

            default: {
                SMESH_ERROR(
                        "[Error]  cu_sshex8_restrict: not implemented for type "
                        "%s "
                        "(code %d)\n",
                        to_string(from_type),
                        from_type);
                return SMESH_FAILURE;
            }
        }
    }
}  // namespace smesh
