#ifndef SMESH_GPU_RUNTIME_CUH
#define SMESH_GPU_RUNTIME_CUH

// The GPU runtime the device path is built against: CUDA, or HIP for AMD GPUs.
//
// smesh's device code is written against the CUDA runtime API and compiles
// unchanged under HIP, whose runtime mirrors it name for name. Under
// SMESH_ENABLE_HIP this header includes HIP and maps the cuda* names smesh uses
// onto their hip* counterparts; otherwise it includes the CUDA runtime. Every
// translation unit that calls the runtime includes it instead of
// <cuda_runtime.h>. Only the names smesh calls are mapped.
//
// SCCD has a header like this one, and a translation unit can include both. A
// macro may be redefined only with the same replacement, so a name mapped in
// both headers must map to exactly the same tokens in each.

#include "smesh_config.hpp"

#if defined(SMESH_ENABLE_HIP)

#include <hip/hip_runtime.h>

#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaPeekAtLastError hipPeekAtLastError

#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaMemGetInfo hipMemGetInfo

#define cudaMemcpy hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemset hipMemset

#define cudaStream_t hipStream_t
#define cudaStreamSynchronize hipStreamSynchronize
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaGetDevice hipGetDevice
#define cudaSetDevice hipSetDevice
#define cudaDeviceGetPCIBusId hipDeviceGetPCIBusId
#define cudaOccupancyMaxPotentialBlockSize hipOccupancyMaxPotentialBlockSize

#define cudaPointerAttributes hipPointerAttribute_t
#define cudaPointerGetAttributes hipPointerGetAttributes
#define cudaMemoryTypeDevice hipMemoryTypeDevice

#if defined(__HIPCC__)
// HIP 6.4 has no __syncwarp (ROCm 7 adds one). smesh uses it as the barrier of
// a tile of 4 or 8 threads exchanging values through shared memory, and a tile
// always lies inside one wavefront. The lanes of a wavefront execute together;
// what the barrier must add is the ordering of their shared-memory accesses,
// which a release and an acquire fence at wavefront scope give. This is what
// ROCm 7's __syncwarp does.
__device__ inline void smesh_syncwarp() {
    __builtin_amdgcn_fence(__ATOMIC_RELEASE, "wavefront");
    __builtin_amdgcn_wave_barrier();
    __builtin_amdgcn_fence(__ATOMIC_ACQUIRE, "wavefront");
}
#define __syncwarp() smesh_syncwarp()
#endif

#else

#include <cuda_runtime.h>

#endif

#endif  // SMESH_GPU_RUNTIME_CUH
