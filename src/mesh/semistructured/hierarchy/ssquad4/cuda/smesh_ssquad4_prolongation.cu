#include "smesh_ssquad4_prolongation.cuh"

#include "smesh_cuda_base.cuh"
#include "smesh_ssquad4_inline.cuh"

#include <cassert>

namespace smesh {

template <typename From, typename To>
__global__ void cu_ssquad4_hierarchical_prolongation_kernel(
    const int level, const ptrdiff_t nelements, const ptrdiff_t stride,
    const idx_t *const SMESH_RESTRICT elements, const int vec_size,
    const ptrdiff_t from_stride, const From *const SMESH_RESTRICT from,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to) {
  const int corners[4] = {
      // Bottom
      cu_ssquad4_lidx(level, 0, 0), cu_ssquad4_lidx(level, level, 0),
      cu_ssquad4_lidx(level, level, level), cu_ssquad4_lidx(level, 0, level)};

  const scalar_t h = 1. / level;

  for (ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x; e < nelements;
       e += blockDim.x * gridDim.x) {
    for (int yi = 0; yi < level + 1; yi++) {
      for (int xi = 0; xi < level + 1; xi++) {
        idx_t idx = elements[cu_ssquad4_lidx(level, xi, yi) * stride + e];

        const scalar_t x = xi * h;
        const scalar_t y = yi * h;

        // Evaluate Hex8 basis functions at x, y, z
        const scalar_t xm = (1 - x);
        const scalar_t ym = (1 - y);

        scalar_t f[4];
        f[0] = xm * ym; // (0, 0)
        f[1] = x * ym;  // (1, 0)
        f[2] = x * y;   // (1, 1)
        f[3] = xm * y;  // (0, 1)

        for (int d = 0; d < vec_size; d++) {
          scalar_t val = 0;

          for (int v = 0; v < 4; v++) {
            const ptrdiff_t global_from_idx =
                (elements[corners[v] * stride + e] * vec_size + d) *
                from_stride;

            assert(from[global_from_idx] == from[global_from_idx]);
            assert(f[v] == f[v]);

            val += f[v] * from[global_from_idx];
          }

          assert(val == val);

          const ptrdiff_t global_to_idx = (idx * vec_size + d) * to_stride;
          to[global_to_idx] = val;
        }
      }
    }
  }
}

template <typename From, typename To>
static int cu_ssquad4_hierarchical_prolongation_tpl(
    const int level, const ptrdiff_t nelements, const ptrdiff_t stride,
    const idx_t *const SMESH_RESTRICT elements, const int vec_size,
    const ptrdiff_t from_stride, const From *const SMESH_RESTRICT from,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to, void *stream) {
  SMESH_DEBUG_SYNCHRONIZE();

  // Hand tuned
  int block_size = 128;
#ifdef SMESH_USE_OCCUPANCY_MAX_POTENTIAL
  {
    int min_grid_size;
    cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &block_size,
        cu_ssquad4_hierarchical_prolongation_kernel<From, To>, 0, 0);
  }
#endif // SMESH_USE_OCCUPANCY_MAX_POTENTIAL

  ptrdiff_t n_blocks =
      MAX(ptrdiff_t(1), (nelements + block_size - 1) / block_size);

  if (stream) {
    cudaStream_t s = *static_cast<cudaStream_t *>(stream);

    cu_ssquad4_hierarchical_prolongation_kernel<From, To>
        <<<n_blocks, block_size, 0, s>>>(level, nelements, stride, elements,
                                         vec_size, from_stride, from, to_stride,
                                         to);
  } else {
    cu_ssquad4_hierarchical_prolongation_kernel<From, To>
        <<<n_blocks, block_size, 0>>>(level, nelements, stride, elements,
                                      vec_size, from_stride, from, to_stride,
                                      to);
  }

  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}

int cu_ssquad4_hierarchical_prolongation(
    const int level, const ptrdiff_t nelements, const ptrdiff_t stride,
    const idx_t *const SMESH_RESTRICT elements, const int vec_size,
    const enum PrimitiveType from_type, const ptrdiff_t from_stride,
    const void *const SMESH_RESTRICT from, const enum PrimitiveType to_type,
    const ptrdiff_t to_stride, void *const SMESH_RESTRICT to, void *stream) {
  assert(from_type == to_type && "TODO mixed types!");
  if (from_type != to_type) {
    return SMESH_FAILURE;
  }

  switch (from_type) {
  case SMESH_DEFAULT: {
    return cu_ssquad4_hierarchical_prolongation_tpl(
        level, nelements, stride, elements, vec_size, from_stride,
        (real_t *)from, to_stride, (real_t *)to, stream);
  }
  case SMESH_FLOAT32: {
    return cu_ssquad4_hierarchical_prolongation_tpl(
        level, nelements, stride, elements, vec_size, from_stride,
        (float *)from, to_stride, (float *)to, stream);
  }
  case SMESH_FLOAT64: {
    return cu_ssquad4_hierarchical_prolongation_tpl(
        level, nelements, stride, elements, vec_size, from_stride,
        (double *)from, to_stride, (double *)to, stream);
  }
  default: {
    SMESH_ERROR(
        "[Error]  cu_ssquad4_prolongation_tpl: not implemented for type %s "
        "(code %d)\n",
        to_string(from_type), from_type);
    return SMESH_FAILURE;
  }
  }
}

static const int TILE_SIZE = 4;
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
    double h = 1. / steps;
    for (int i = 0; i < nodes; i++) {
      S_host[0 * stride + i] = (1 - h * i);
      S_host[1 * stride + i] = h * i;
    }

    auto nbytes = S_host.size() * sizeof(T);

    SMESH_CUDA_CHECK(cudaMalloc((void **)&data, nbytes));
    SMESH_CUDA_CHECK(
        cudaMemcpy(data, S_host.data(), nbytes, cudaMemcpyHostToDevice));
  }

  ~ShapeInterpolation() { cudaFree(data); }
};

// Even TO sub-elements are used to interpolate from FROM sub-elements
template <typename From, typename To>
__global__ void cu_ssquad4_prolongate_kernel(
    const ptrdiff_t nelements, const ptrdiff_t stride, const int from_level,
    const int from_level_stride, idx_t *const SMESH_RESTRICT from_elements,
    const int to_level, const int to_level_stride,
    idx_t *const SMESH_RESTRICT to_elements, const int vec_size,
    const ptrdiff_t from_stride, const From *const SMESH_RESTRICT from,
    const ptrdiff_t to_stride, To *const SMESH_RESTRICT to) {
  static_assert(TILE_SIZE == 4, "This only works with tile size 8!");

  // Uunsigned char necessary for multiple instantiations
  __shared__ unsigned char cu_buff[];

  const int step_factor = to_level / from_level;

  // Tile number in group
  const int tile = threadIdx.x >> 2;
  const int n_tiles = blockDim.x >> 2;
  const int sub_idx = threadIdx.x & 0x3;

  From *in = (From *)&cu_buff[tile * TILE_SIZE * sizeof(From)];

  // quad4 idx
  const int xi = sub_idx & 0x1;        // equivalent to sub_idx % 2
  const int yi = (sub_idx >> 1) & 0x1; // equivalent to (sub_idx / 2) % 2

  assert(n_tiles * TILE_SIZE == blockDim.x);

  // 1 macro element per tile
  const ptrdiff_t e = blockIdx.x * n_tiles + tile;

  const To between_h = (From)from_level / (To)to_level;

  const int to_even = is_even(to_level);
  const int from_nloops = from_level + to_even;
  const int to_nloops = to_level + to_even;

  // Vector loop
  for (int d = 0; d < vec_size; d++) {
    // loop on all FROM micro elements

    for (int from_yi = 0; from_yi < from_nloops; from_yi++) {
      const int off_from_yi = (from_yi + yi);

      for (int from_xi = 0; from_xi < from_nloops; from_xi++) {
        const int off_from_xi = (from_xi + xi);

        const bool from_exists = e < nelements && off_from_yi <= from_level &&
                                 off_from_xi <= from_level;

        // Wait for shared memory transactions to be finished
        __syncwarp();

        // Gather
        if (from_exists) {
          const int idx_from = cu_ssquad4_lidx(from_level * from_level_stride,
                                               off_from_xi * from_level_stride,
                                               off_from_yi * from_level_stride);

          const idx_t gidx = from_elements[idx_from * stride + e];
          const ptrdiff_t idx = (gidx * vec_size + d) * from_stride;
          in[sub_idx] = from[idx];
        } else {
          in[sub_idx] = 0;
        }

        // Wait for in to be filled
        __syncwarp();

        int start_yi = from_yi * step_factor;
        start_yi += is_odd(start_yi); // Skip odd numbers

        int start_xi = from_xi * step_factor;
        start_xi += is_odd(start_xi); // Skip odd numbers

        const int end_yi = MIN(to_nloops, start_yi + step_factor);
        const int end_xi = MIN(to_nloops, start_xi + step_factor);

        // sub-loop on even TO micro-elements

        for (int to_yi = start_yi; to_yi < end_yi; to_yi += 2) {
          for (int to_xi = start_xi; to_xi < end_xi; to_xi += 2) {
            const int off_to_yi = (to_yi + yi);
            const int off_to_xi = (to_xi + xi);

            const To x = (off_to_xi - from_xi * step_factor) * between_h;
            const To y = (off_to_yi - from_yi * step_factor) * between_h;

            assert(x >= 0);
            assert(x <= 1);
            assert(y >= 0);
            assert(y <= 1);

            // This requires 64 bytes on the stack frame
            // Cartesian order
            To f[4] = {// Bottom
                       (1 - x) * (1 - y), x * (1 - y), (1 - x) * y, x * y};

#ifndef NDEBUG
            To pou = 0;
            for (int i = 0; i < 4; i++) {
              pou += f[i];
            }

            assert(fabs(1 - pou) < 1e-8);
#endif

            To out = 0;
            for (int v = 0; v < 4; v++) {
              const int round_robin = ROUND_ROBIN(v, sub_idx);
              // There should be no bank conflicts due to round robin
              out += f[round_robin] * in[round_robin];
            }

            // Check if not ghost nodes for scatter assign
            const bool to_exists =
                e < nelements && off_to_yi <= to_level && off_to_xi <= to_level;

            if (to_exists) {
              // set the value to the output
              const int idx_to = cu_ssquad4_lidx(to_level * to_level_stride,
                                                 off_to_xi * to_level_stride,
                                                 off_to_yi * to_level_stride);

              to[((ptrdiff_t)to_elements[idx_to * stride + e] * vec_size + d) *
                 to_stride] = out;
            }
          }
        }
      }
    }
  }
}

template <typename From, typename To>
int cu_ssquad4_prolongate_tpl(const ptrdiff_t nelements, const ptrdiff_t stride,
                              const int from_level, const int from_level_stride,
                              idx_t *const SMESH_RESTRICT from_elements,
                              const int to_level, const int to_level_stride,
                              idx_t *const SMESH_RESTRICT to_elements,
                              const int vec_size, const ptrdiff_t from_stride,
                              const From *const SMESH_RESTRICT from,
                              const ptrdiff_t to_stride,
                              To *const SMESH_RESTRICT to, void *stream) {
  SMESH_DEBUG_SYNCHRONIZE();

  const int block_size = 128;
  ptrdiff_t n_blocks =
      MAX(ptrdiff_t(1),
          (nelements + block_size / TILE_SIZE - 1) / (block_size / TILE_SIZE));
  size_t shared_mem_size = block_size * sizeof(From);

  if (stream) {
    cudaStream_t s = *static_cast<cudaStream_t *>(stream);

    cu_ssquad4_prolongate_kernel<From, To>
        <<<n_blocks, block_size, shared_mem_size, s>>>(
            nelements, stride, from_level, from_level_stride, from_elements,
            to_level, to_level_stride, to_elements, vec_size, from_stride, from,
            to_stride, to);
  } else {
    cu_ssquad4_prolongate_kernel<From, To>
        <<<n_blocks, block_size, shared_mem_size>>>(
            nelements, stride, from_level, from_level_stride, from_elements,
            to_level, to_level_stride, to_elements, vec_size, from_stride, from,
            to_stride, to);
  }

  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}

int cu_ssquad4_prolongate(
    const ptrdiff_t nelements, const ptrdiff_t stride, const int from_level,
    const int from_level_stride, idx_t *const SMESH_RESTRICT from_elements,
    const int to_level, const int to_level_stride,
    idx_t *const SMESH_RESTRICT to_elements, const int vec_size,
    const enum PrimitiveType from_type, const ptrdiff_t from_stride,
    const void *const SMESH_RESTRICT from, const enum PrimitiveType to_type,
    const ptrdiff_t to_stride, void *const SMESH_RESTRICT to, void *stream) {
  assert(from_type == to_type && "TODO mixed types!");
  if (from_type != to_type) {
    return SMESH_FAILURE;
  }

  switch (from_type) {
  case SMESH_DEFAULT: {
    return cu_ssquad4_prolongate_tpl<real_t, real_t>(
        nelements, stride, from_level, from_level_stride, from_elements,
        to_level, to_level_stride, to_elements, vec_size, from_stride,
        (real_t *)from, to_stride, (real_t *)to, stream);
  }
  case SMESH_FLOAT32: {
    return cu_ssquad4_prolongate_tpl<float, float>(
        nelements, stride, from_level, from_level_stride, from_elements,
        to_level, to_level_stride, to_elements, vec_size, from_stride,
        (float *)from, to_stride, (float *)to, stream);
  }
  case SMESH_FLOAT64: {
    return cu_ssquad4_prolongate_tpl<double, double>(
        nelements, stride, from_level, from_level_stride, from_elements,
        to_level, to_level_stride, to_elements, vec_size, from_stride,
        (double *)from, to_stride, (double *)to, stream);
  }

  default: {
    SMESH_ERROR("[Error]  cu_ssquad4_prolongate: not implemented for type "
                "%s "
                "(code %d)\n",
                to_string(from_type), from_type);
    return SMESH_FAILURE;
  }
  }
}

} // namespace smesh