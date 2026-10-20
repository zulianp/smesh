#include "smesh_exchange_kernels.hpp"

#include "smesh_base.hpp"
#include "smesh_cuda_base.cuh"

#include <algorithm>
#include <cstdint>

namespace smesh {
namespace device {

namespace {
constexpr int N_THREADS = 256;

int n_blocks(const ptrdiff_t n) {
  return (int)std::max<ptrdiff_t>(
      1, std::min<ptrdiff_t>((n + N_THREADS - 1) / N_THREADS, 4096));
}

template <typename W>
__global__ void gather_kernel(const ptrdiff_t n, const ptrdiff_t block_size,
                              const idx_t *const SMESH_RESTRICT idx,
                              const W *const SMESH_RESTRICT in,
                              W *const SMESH_RESTRICT out) {
  for (ptrdiff_t k = blockIdx.x * (ptrdiff_t)blockDim.x + threadIdx.x;
       k < n * block_size; k += (ptrdiff_t)blockDim.x * gridDim.x) {
    const ptrdiff_t i = k / block_size;
    out[k] = in[(ptrdiff_t)idx[i] * block_size + (k - i * block_size)];
  }
}

template <typename W>
__global__ void scatter_kernel(const ptrdiff_t n, const ptrdiff_t block_size,
                               const idx_t *const SMESH_RESTRICT idx,
                               const W *const SMESH_RESTRICT in,
                               W *const SMESH_RESTRICT out) {
  for (ptrdiff_t k = blockIdx.x * (ptrdiff_t)blockDim.x + threadIdx.x;
       k < n * block_size; k += (ptrdiff_t)blockDim.x * gridDim.x) {
    const ptrdiff_t i = k / block_size;
    out[(ptrdiff_t)idx[i] * block_size + (k - i * block_size)] = in[k];
  }
}

template <typename W> struct Gather {
  static void launch(const ptrdiff_t n, const ptrdiff_t bs, const idx_t *idx,
                     const void *in, void *out) {
    gather_kernel<W><<<n_blocks(n * bs), N_THREADS>>>(
        n, bs, idx, static_cast<const W *>(in), static_cast<W *>(out));
  }
};

template <typename W> struct Scatter {
  static void launch(const ptrdiff_t n, const ptrdiff_t bs, const idx_t *idx,
                     const void *in, void *out) {
    scatter_kernel<W><<<n_blocks(n * bs), N_THREADS>>>(
        n, bs, idx, static_cast<const W *>(in), static_cast<W *>(out));
  }
};

// Copies words of the element's size: what the element is does not matter.
template <template <typename> class Copy>
int copy(const ptrdiff_t n, const ptrdiff_t block_size,
         const size_t element_size, const idx_t *const idx,
         const void *const in, void *const out) {
  if (n <= 0) return SMESH_SUCCESS;
  switch (element_size) {
  case 1:
    Copy<std::uint8_t>::launch(n, block_size, idx, in, out);
    break;
  case 2:
    Copy<std::uint16_t>::launch(n, block_size, idx, in, out);
    break;
  case 4:
    Copy<std::uint32_t>::launch(n, block_size, idx, in, out);
    break;
  case 8:
    Copy<std::uint64_t>::launch(n, block_size, idx, in, out);
    break;
  default:
    SMESH_ERROR("device copy by index: elements of %zu bytes\n", element_size);
    return SMESH_FAILURE;
  }
  SMESH_CUDA_CHECK(cudaPeekAtLastError());
  SMESH_DEBUG_SYNCHRONIZE();
  return SMESH_SUCCESS;
}
} // namespace

int gather_by_index(const ptrdiff_t n, const ptrdiff_t block_size,
                    const size_t element_size, const idx_t *const idx,
                    const void *const in, void *const out) {
  return copy<Gather>(n, block_size, element_size, idx, in, out);
}

int scatter_by_index(const ptrdiff_t n, const ptrdiff_t block_size,
                     const size_t element_size, const idx_t *const idx,
                     const void *const in, void *const out) {
  return copy<Scatter>(n, block_size, element_size, idx, in, out);
}

int synchronize() {
  SMESH_CUDA_CHECK(cudaDeviceSynchronize());
  return SMESH_SUCCESS;
}

} // namespace device
} // namespace smesh
