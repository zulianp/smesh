#ifndef SMESH_CUDA_BASE_CUH
#define SMESH_CUDA_BASE_CUH

#include "smesh_gpu_runtime.cuh"

#if !defined(SMESH_ENABLE_HIP)
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_fp16.h>
#endif

#include <stdio.h>
#include <assert.h>

inline void sfem_cuda_check(cudaError_t code, const char* file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "cuda_check: %s %s:%d\n", cudaGetErrorString(code), file, line);
        assert(!code);
        if (abort) exit(code);
    }
}

#define SMESH_CUDA_CHECK(ans) \
    { sfem_cuda_check((ans), __FILE__, __LINE__); }

#ifndef NDEBUG
#define SMESH_DEBUG_SYNCHRONIZE()                \
    do {                                        \
        cudaDeviceSynchronize();                \
        SMESH_CUDA_CHECK(cudaPeekAtLastError()); \
    } while (0)
#else
#define SMESH_DEBUG_SYNCHRONIZE()
#endif

// NVTX ranges are NVIDIA's; a HIP build traces nothing on the device timeline.
#if !defined(SMESH_ENABLE_HIP)
#define SMESH_ENABLE_NVTX
#endif

#ifdef SMESH_ENABLE_NVTX
// NVTX v3, header-only; CUDA 12.9 removed the v2 header (nvToolsExt.h), which
// earlier toolkits also ship.
#if __has_include(<nvtx3/nvToolsExt.h>)
#include <nvtx3/nvToolsExt.h>
#else
#include <nvToolsExt.h>
#endif
namespace sfem {
    namespace details {
        class Tracer {
        public:
            Tracer(const char* name) { nvtxRangePushA(name); }
            ~Tracer() { nvtxRangePop(); }
        };
    }  // namespace details
}  // namespace sfem

#define SMESH_NVTX_SCOPE(name) sfem::details::Tracer uniq_name_using_macros(name);
#define SMESH_RANGE_PUSH(name_) \
    do {                       \
        nvtxRangePushA(name_); \
    } while (0)
#define SMESH_RANGE_POP() \
    do {                 \
        nvtxRangePop();  \
    } while (0)

#else //SMESH_ENABLE_NVTX

#define SMESH_NVTX_SCOPE(name)
#define SMESH_RANGE_PUSH(name_)
#define SMESH_RANGE_POP()
 
#endif //SMESH_ENABLE_NVTX

// Launch bounds configuration for device kernels. HIP's second argument is the
// minimum number of waves per execution unit, not blocks per multiprocessor, so
// a HIP build passes only the block size.
#if defined(SMESH_ENABLE_HIP)
#define SMESH_LAUNCH_BOUNDS(threads_per_block, min_blocks_per_sm) \
    __launch_bounds__(threads_per_block)
#else
#define SMESH_LAUNCH_BOUNDS(threads_per_block, min_blocks_per_sm) \
    __launch_bounds__(threads_per_block, min_blocks_per_sm)
#endif
#endif  // SMESH_CUDA_BASE_CUH
