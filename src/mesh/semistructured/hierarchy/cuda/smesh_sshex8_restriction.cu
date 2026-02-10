#include "cu_sshex8_interpolate.h"

#include "sfem_cuda_base.h"

#include "cu_sshex8_inline.hpp"

#include <cSMESH_ASSERT>
#include <cstdio>
#include <vector>

static const int TILE_SIZE = 8;
#define ROUND_ROBIN(val, shift) ((val + shift) & (TILE_SIZE - 1))
#define ROUND_ROBIN_2(val, shift) ((val + shift) & (2 - 1))

static inline __device__ int is_even(const int v) { return !(v & 1); }
static inline __device__ int is_odd(const int v) { return (v & 1); }

template <typename T> class ShapeInterpolation {
public:
  T *data{nullptr};
  size_t nodes{0};
  int stride{0};

  ShapeInterpolation(const int steps, const int padding = 0) {
    nodes = (steps + 1);
    stride = nodes + padding;
    std::vector<T> S_host(2 * stride, 0);
    f64 h = 1. / steps;
    for (int i = 0; i < nodes; i++) {
      S_host[0 * stride + i] = (1 - h * i);
      S_host[1 * stride + i] = h * i;
    }

    auto nbytes = S_host.size() * sizeof(T);

    SMESH_CUDA_CHECK(cudaMalloc((void **)&data, nbytes));
    SMESH_CUDA_CHECK(
        cudaMemcpy(data, S_host.data(), nbytes, cudaMemcpyHostToDevice));

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
  }

  ~ShapeInterpolation() { cudaFree(data); }
};

// PROLONGATION

template <typename From, typename To, typename idx_t>
__global__ void cu_sshex8_hierarchical_prolongation_kernel(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int vec_size, const ptrdiff_t from_stride,
    const From *const SMESH_RESTRICT from, const ptrdiff_t to_stride,
    To *const SMESH_RESTRICT to) {
  const int corners[8] = {// Bottom
                          cu_sshex8_lidx(level, 0, 0, 0),
                          cu_sshex8_lidx(level, level, 0, 0),
                          cu_sshex8_lidx(level, level, level, 0),
                          cu_sshex8_lidx(level, 0, level, 0),
                          // Top
                          cu_sshex8_lidx(level, 0, 0, level),
                          cu_sshex8_lidx(level, level, 0, level),
                          cu_sshex8_lidx(level, level, level, level),
                          cu_sshex8_lidx(level, 0, level, level)};

  const scalar_t h = 1. / level;

  for (ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x; e < nelements;
       e += blockDim.x * gridDim.x) {
    for (int zi = 0; zi < level + 1; zi++) {
      for (int yi = 0; yi < level + 1; yi++) {
        for (int xi = 0; xi < level + 1; xi++) {
          idx_t idx = elements[cu_sshex8_lidx(level, xi, yi, zi)][e];

          const scalar_t x = xi * h;
          const scalar_t y = yi * h;
          const scalar_t z = zi * h;

          // Evaluate Hex8 basis functions at x, y, z
          const scalar_t xm = (1 - x);
          const scalar_t ym = (1 - y);
          const scalar_t zm = (1 - z);

          scalar_t f[8];
          f[0] = xm * ym * zm; // (0, 0, 0)
          f[1] = x * ym * zm;  // (1, 0, 0)
          f[2] = x * y * zm;   // (1, 1, 0)
          f[3] = xm * y * zm;  // (0, 1, 0)
          f[4] = xm * ym * z;  // (0, 0, 1)
          f[5] = x * ym * z;   // (1, 0, 1)
          f[6] = x * y * z;    // (1, 1, 1)
          f[7] = xm * y * z;   // (0, 1, 1)

          for (int d = 0; d < vec_size; d++) {
            scalar_t val = 0;

            for (int v = 0; v < 8; v++) {
              const ptrdiff_t global_from_idx =
                  (elements[corners[v]][e] * vec_size + d) * from_stride;

              SMESH_ASSERT(from[global_from_idx] == from[global_from_idx]);
              SMESH_ASSERT(f[v] == f[v]);

              val += f[v] * from[global_from_idx];
            }

            SMESH_ASSERT(val == val);

            const ptrdiff_t global_to_idx = (idx * vec_size + d) * to_stride;
            to[global_to_idx] = val;
          }
        }
      }
    }
  }
}

template <typename From, typename To, typename idx_t>
static int cu_sshex8_hierarchical_prolongation_tpl(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int vec_size, const ptrdiff_t from_stride,
    const From *const SMESH_RESTRICT from, const ptrdiff_t to_stride,
    To *const SMESH_RESTRICT to, void *stream) {
  SMESH_DEBUG_SYNCHRONIZE();

  // Hand tuned
  int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
  {
    int min_grid_size;
    cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &block_size,
        cu_sshex8_hierarchical_prolongation_kernel<From, To>, 0, 0);
  }
#endif // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

  ptrdiff_t n_blocks =
      MAX(ptrdiff_t(1), (nelements + block_size - 1) / block_size);

  if (stream) {
    cudaStream_t s = *static_cast<cudaStream_t *>(stream);

    cu_sshex8_hierarchical_prolongation_kernel<From, To>
        <<<n_blocks, block_size, 0, s>>>(level, nelements, elements, vec_size,
                                         from_stride, from, to_stride, to);
  } else {
    cu_sshex8_hierarchical_prolongation_kernel<From, To>
        <<<n_blocks, block_size, 0>>>(level, nelements, elements, vec_size,
                                      from_stride, from, to_stride, to);
  }

  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}

template <typename idx_t>
int cu_sshex8_hierarchical_prolongation(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int vec_size, const enum RealType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum RealType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream) {
  SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
  if (from_type != to_type) {
    return SMESH_FAILURE;
  }

  switch (from_type) {
  case SMESH_REAL_DEFAULT: {
    return cu_sshex8_hierarchical_prolongation_tpl(
        level, nelements, elements, vec_size, from_stride, (real_t *)from,
        to_stride, (real_t *)to, stream);
  }
  case SMESH_FLOAT32: {
    return cu_sshex8_hierarchical_prolongation_tpl(
        level, nelements, elements, vec_size, from_stride, (f32 *)from,
        to_stride, (f32 *)to, stream);
  }
  case SMESH_FLOAT64: {
    return cu_sshex8_hierarchical_prolongation_tpl(
        level, nelements, elements, vec_size, from_stride, (f64 *)from,
        to_stride, (f64 *)to, stream);
  }
  default: {
    SMESH_ERROR(
        "[Error]  cu_sshex8_prolongation_tpl: not implemented for type %s "
        "(code %d)\n",
        real_type_to_string(from_type), from_type);
    return SMESH_FAILURE;
  }
  }
}

template <typename From, typename To, typename idx_t>
__global__ void cu_sshex8_hierarchical_restriction_kernel(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT e2n_count, const int vec_size,
    const ptrdiff_t from_stride, const From *const SMESH_RESTRICT from,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to) {
  const int corners[8] = {// Bottom
                          cu_sshex8_lidx(level, 0, 0, 0),
                          cu_sshex8_lidx(level, level, 0, 0),
                          cu_sshex8_lidx(level, level, level, 0),
                          cu_sshex8_lidx(level, 0, level, 0),
                          // Top
                          cu_sshex8_lidx(level, 0, 0, level),
                          cu_sshex8_lidx(level, level, 0, level),
                          cu_sshex8_lidx(level, level, level, level),
                          cu_sshex8_lidx(level, 0, level, level)};

  const scalar_t h = 1. / level;
  scalar_t acc[8];

  for (ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x; e < nelements;
       e += blockDim.x * gridDim.x) {
    for (int d = 0; d < vec_size; d++) {
      for (int i = 0; i < 8; i++) {
        acc[i] = 0;
      }

      for (int zi = 0; zi < level + 1; zi++) {
        for (int yi = 0; yi < level + 1; yi++) {
          for (int xi = 0; xi < level + 1; xi++) {
            const int lidx = cu_sshex8_lidx(level, xi, yi, zi);
            const ptrdiff_t idx = elements[lidx][e];

            const scalar_t x = xi * h;
            const scalar_t y = yi * h;
            const scalar_t z = zi * h;

            // Evaluate Hex8 basis functions at x, y, z
            const scalar_t xm = (1 - x);
            const scalar_t ym = (1 - y);
            const scalar_t zm = (1 - z);

            scalar_t f[8];
            f[0] = xm * ym * zm; // (0, 0, 0)
            f[1] = x * ym * zm;  // (1, 0, 0)
            f[2] = x * y * zm;   // (1, 1, 0)
            f[3] = xm * y * zm;  // (0, 1, 0)
            f[4] = xm * ym * z;  // (0, 0, 1)
            f[5] = x * ym * z;   // (1, 0, 1)
            f[6] = x * y * z;    // (1, 1, 1)
            f[7] = xm * y * z;   // (0, 1, 1)

            const ptrdiff_t global_from_idx =
                (idx * vec_size + d) * from_stride;
            const scalar_t val = from[global_from_idx] / e2n_count[idx];

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
        const ptrdiff_t global_to_idx =
            (elements[corners[v]][e] * vec_size + d) * to_stride;
        atomicAdd(&to[global_to_idx], acc[v]);
      }
    }
  }
}

template <typename From, typename To, typename idx_t>
static int cu_sshex8_hierarchical_restriction_tpl(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const ptrdiff_t from_stride,
    const From *const SMESH_RESTRICT from, const ptrdiff_t to_stride,
    To *const SMESH_RESTRICT to, void *stream) {
  SMESH_DEBUG_SYNCHRONIZE();

  // Hand tuned
  int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
  {
    int min_grid_size;
    cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &block_size,
        cu_sshex8_hierarchical_restriction_kernel<From, To, idx_t>, 0, 0);
  }
#endif // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

  ptrdiff_t n_blocks =
      MAX(ptrdiff_t(1), (nelements + block_size - 1) / block_size);

  if (stream) {
    cudaStream_t s = *static_cast<cudaStream_t *>(stream);

    cu_sshex8_hierarchical_restriction_kernel<From, To, idx_t>
        <<<n_blocks, block_size, 0, s>>>(
            level, nelements, elements, element_to_node_incidence_count,
            vec_size, from_stride, from, to_stride, to);
  } else {
    cu_sshex8_hierarchical_restriction_kernel<From, To, idx_t>
        <<<n_blocks, block_size, 0>>>(level, nelements, elements,
                                      element_to_node_incidence_count, vec_size,
                                      from_stride, from, to_stride, to);
  }

  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}

extern int cu_sshex8_hierarchical_restriction(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const u16 *const SMESH_RESTRICT element_to_node_incidence_count,
    const int vec_size, const enum RealType from_type,
    const ptrdiff_t from_stride, const void *const SMESH_RESTRICT from,
    const enum RealType to_type, const ptrdiff_t to_stride,
    void *const SMESH_RESTRICT to, void *stream) {
  SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
  if (from_type != to_type) {
    return SMESH_FAILURE;
  }

  switch (from_type) {
  case SMESH_REAL_DEFAULT: {
    return cu_sshex8_hierarchical_restriction_tpl(
        level, nelements, elements, element_to_node_incidence_count, vec_size,
        from_stride, (real_t *)from, to_stride, (real_t *)to, stream);
  }
  case SMESH_FLOAT32: {
    return cu_sshex8_hierarchical_restriction_tpl(
        level, nelements, elements, element_to_node_incidence_count, vec_size,
        from_stride, (f32 *)from, to_stride, (f32 *)to, stream);
  }
  case SMESH_FLOAT64: {
    return cu_sshex8_hierarchical_restriction_tpl(
        level, nelements, elements, element_to_node_incidence_count, vec_size,
        from_stride, (f64 *)from, to_stride, (f64 *)to, stream);
  }
  default: {
    SMESH_ERROR("[Error]  cu_sshex8_prolongation_tpl: not implemented for type "
                "%s "
                "(code %d)\n",
                real_type_to_string(from_type), from_type);
    return SMESH_FAILURE;
  }
  }
}

template <typename From, typename To>
__global__ void sshex8_hierarchical_restriction_kernel(
    const ptrdiff_t nelements, const int from_level,
    const int from_level_stride,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
    const u16 *const SMESH_RESTRICT from_element_to_node_incidence_count,
    const int to_level, const int to_level_stride,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT to_elements,
    const To *const SMESH_RESTRICT S, const int vec_size,
    const enum RealType from_type, const ptrdiff_t from_stride,
    const From *const SMESH_RESTRICT from, const enum RealType to_type,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to) {
  static_assert(
      TILE_SIZE == 8,
      "This only works with tile size 8 because the implementation assumes a "
      "fixed tile size for shared memory "
      "layout and indexing.");

  // Unsigned char necessary for multiple template instantiations of this kernel
  extern __shared__ unsigned char cu_buff[];

  const int step_factor = from_level / to_level;

  // Tile number in group
  const int tile = threadIdx.x >> 3;     // same as threadIdx.x / 8
  const int n_tiles = blockDim.x >> 3;   // same as blockDim.x / 8
  const int sub_idx = threadIdx.x & 0x7; // same as threadIdx.x % 8

  From *in = (From *)&cu_buff[tile * TILE_SIZE * sizeof(From)];

  // hex8 idx
  const int xi = sub_idx & 0x1;        // equivalent to sub_idx % 2
  const int yi = (sub_idx >> 1) & 0x1; // equivalent to (sub_idx / 2) % 2
  const int zi = (sub_idx >> 2);       // equivalent to sub_idx / 4
  SMESH_ASSERT(n_tiles * TILE_SIZE == blockDim.x);

  // // 1 macro element per tile
  const ptrdiff_t e = blockIdx.x * n_tiles + tile;
  const int from_even = is_even(from_level);
  const int to_nloops = to_level + from_even;

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

                  const bool from_exists =
                      e < nelements &&             // Check element exists
                      off_from_zi <= from_level && //
                      off_from_yi <= from_level && //
                      off_from_xi <= from_level;

                  __syncwarp();
                  if (from_exists) {
                    const int idx_from =
                        cu_sshex8_lidx(from_level * from_level_stride,
                                       off_from_xi * from_level_stride,
                                       off_from_yi * from_level_stride,
                                       off_from_zi * from_level_stride);

                    const idx_t gidx = from_elements[idx_from][e];
                    const ptrdiff_t idx = (gidx * vec_size + d) * from_stride;
                    in[sub_idx] =
                        from[idx] / from_element_to_node_incidence_count[gidx];
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

          const bool exists = e < nelements &&         // Check element exists
                              off_to_zi <= to_level && //
                              off_to_yi <= to_level && //
                              off_to_xi <= to_level;
          if (exists) {
            const int idx_to = cu_sshex8_lidx(
                to_level * to_level_stride, off_to_xi * to_level_stride,
                off_to_yi * to_level_stride, off_to_zi * to_level_stride);

            atomicAdd(&to[((ptrdiff_t)to_elements[idx_to][e] * vec_size + d) *
                          to_stride],
                      acc);
          }
        }
      }
    }
  }
}

template <typename From, typename To, typename idx_t>
int sshex8_hierarchical_restriction_tpl(
    const ptrdiff_t nelements, const int from_level,
    const int from_level_stride, idx_t **const SMESH_RESTRICT from_elements,
    const u16 *const SMESH_RESTRICT from_element_to_node_incidence_count,
    const int to_level, const int to_level_stride,
    idx_t **const SMESH_RESTRICT to_elements, const int vec_size,
    const enum RealType from_type, const ptrdiff_t from_stride,
    const From *const SMESH_RESTRICT from, const enum RealType to_type,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to, void *stream) {
  SMESH_DEBUG_SYNCHRONIZE();

  // Hand tuned
  int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
  {
    int min_grid_size;
    cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &block_size, sshex8_hierarchical_restriction_kernel<From, To>, 0, 0);
  }
#endif // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

  ptrdiff_t n_blocks =
      MAX(ptrdiff_t(1),
          (nelements + block_size / TILE_SIZE - 1) / (block_size / TILE_SIZE));

  ShapeInterpolation<To> S(from_level / to_level, from_level % 2 == 0);

  size_t shared_mem_size = block_size * sizeof(From);

  if (stream) {
    cudaStream_t s = *static_cast<cudaStream_t *>(stream);

    sshex8_hierarchical_restriction_kernel<From, To>
        <<<n_blocks, block_size, shared_mem_size, s>>>(
            nelements, from_level, from_level_stride, from_elements,
            from_element_to_node_incidence_count, to_level, to_level_stride,
            to_elements, S.data, vec_size, from_type, from_stride, from,
            to_type, to_stride, to);
  } else {
    sshex8_hierarchical_restriction_kernel<From, To>
        <<<n_blocks, block_size, shared_mem_size>>>(
            nelements, from_level, from_level_stride, from_elements,
            from_element_to_node_incidence_count, to_level, to_level_stride,
            to_elements, S.data, vec_size, from_type, from_stride, from,
            to_type, to_stride, to);
  }

  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}

template <typename idx_t>
int sshex8_hierarchical_restriction(
    const ptrdiff_t nelements, const int from_level,
    const int from_level_stride, const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT from_elements,
    const u16 *const SMESH_RESTRICT from_element_to_node_incidence_count,
    const int to_level, const int to_level_stride,
    idx_t **const SMESH_RESTRICT to_elements, const int vec_size,
    const enum RealType from_type, const ptrdiff_t from_stride,
    const void *const SMESH_RESTRICT from, const enum RealType to_type,
    const ptrdiff_t to_stride, void *const SMESH_RESTRICT to, void *stream) {
  SMESH_ASSERT(from_type == to_type && "TODO mixed types!");
  if (from_type != to_type) {
    return SMESH_FAILURE;
  }

  switch (from_type) {
  case SMESH_REAL_DEFAULT: {
    return sshex8_hierarchical_restriction_tpl<real_t, real_t, idx_t>(
        nelements, from_level, from_level_stride, from_elements,
        from_element_to_node_incidence_count, to_level, to_level_stride,
        to_elements, vec_size, from_type, from_stride, (real_t *)from, to_type,
        to_stride, (real_t *)to, stream);
  }
  case SMESH_FLOAT32: {
    return sshex8_hierarchical_restriction_tpl<f32, f32, idx_t>(
        nelements, from_level, from_level_stride, from_elements,
        from_element_to_node_incidence_count, to_level, to_level_stride,
        to_elements, vec_size, from_type, from_stride, (f32 *)from, to_type,
        to_stride, (f32 *)to, stream);
  }
  case SMESH_FLOAT64: {
    return sshex8_hierarchical_restriction_tpl<f64, f64, idx_t>(
        nelements, from_level, from_level_stride, from_elements,
        from_element_to_node_incidence_count, to_level, to_level_stride,
        to_elements, vec_size, from_type, from_stride, (f64 *)from, to_type,
        to_stride, (f64 *)to, stream);
  }

  default: {
    SMESH_ERROR("[Error]  sshex8_hierarchical_restriction: not implemented for type "
                "%s "
                "(code %d)\n",
                real_type_to_string(from_type), from_type);
    return SMESH_FAILURE;
  }
  }
}
