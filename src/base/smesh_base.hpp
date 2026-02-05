#ifndef SMESH_BASE_HPP
#define SMESH_BASE_HPP

#include <assert.h>
#include <chrono>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __GNUC__
#include <execinfo.h>
#include <sys/time.h>
#endif
#if defined(__HIPCC__)
#include <hip/driver_types.h>
#include <hip/hip_runtime.h>
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#endif

#define SMESH_SUCCESS 0
#define SMESH_FAILURE 1

#ifdef SMESH_ENABLE_MPI
#include <mpi.h>
#define SMESH_ABORT() MPI_Abort(MPI_COMM_WORLD, SMESH_FAILURE)
#else
#define SMESH_ABORT() abort()
#endif

#ifndef __PRETTY_FUNCTION__
#if defined(__func__)
#define __PRETTY_FUNCTION__ __func__
#else
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif
#endif

#define SMESH_READ_ENV(name, conversion)                                       \
  do {                                                                         \
    char *var = getenv(#name);                                                 \
    if (var) {                                                                 \
      name = conversion(var);                                                  \
    }                                                                          \
  } while (0)

#define SMESH_REQUIRE_ENV(name, conversion)                                    \
  do {                                                                         \
    char *var = getenv(#name);                                                 \
    if (var) {                                                                 \
      name = conversion(var);                                                  \
    } else {                                                                   \
      fprintf(stderr, "[Error] %s is required (%s:%d)", #name, __FILE__,       \
              __LINE__);                                                       \
      assert(0);                                                               \
      SMESH_ABORT();                                                           \
    }                                                                          \
  } while (0)

#define SMESH_ERROR(...)                                                       \
  do {                                                                         \
    fprintf(stderr, __VA_ARGS__);                                              \
    fprintf(stderr, "Aborting at %s:%d\n", __FILE__, __LINE__);                \
    fflush(stderr);                                                            \
    assert(0);                                                                 \
    SMESH_ABORT();                                                             \
  } while (0)

#define SMESH_IMPLEMENT_ME() SMESH_ERROR("Implement me!\n")

#ifdef NDEBUG
#define SMESH_INLINE inline
#define SMESH_FORCE_INLINE inline __attribute__((always_inline))
#else
#define SMESH_INLINE
#define SMESH_FORCE_INLINE
#endif

#define SMESH_UNUSED(var) (void)var
#ifndef _WIN32
#define SMESH_RESTRICT __restrict__
#else
#define SMESH_RESTRICT __restrict
#endif

#if defined(__CUDACC__)
#define SMESH_DEVICE __device__
#define SMESH_HOST __host__
#elif defined(__HIPCC__)
#define SMESH_DEVICE __device__
#define SMESH_HOST __host__
#else
#define SMESH_DEVICE /* ignore */
#define SMESH_HOST   /* ignore */
#endif

#define SMESH_BOTH SMESH_HOST SMESH_DEVICE

namespace smesh {
inline SMESH_BOTH int32_t div_round_up(int32_t a, int32_t b) {
  return (a + b - 1) / b;
}

inline SMESH_BOTH uint32_t div_round_up(uint32_t a, uint32_t b) {
  return (a + b - 1) / b;
}

inline SMESH_BOTH int64_t div_round_up(int64_t a, int64_t b) {
  return (a + b - 1) / b;
}

inline SMESH_BOTH uint64_t div_round_up(uint64_t a, uint64_t b) {
  return (a + b - 1) / b;
}

inline double time_milliseconds() {
  auto now = std::chrono::high_resolution_clock::now();
  auto duration = now.time_since_epoch();
  return std::chrono::duration<double, std::milli>(duration).count();
}

inline double time_seconds() { return time_milliseconds() / 1000.0; }

#if defined(__CUDACC__) || defined(__HIPCC__)
#define SMESH_CUDA_CHECK(call)                                                 \
  {                                                                            \
    cudaError_t rc = call;                                                     \
    if (rc != cudaSuccess) {                                                   \
      fprintf(stderr, "CUDA call (%s) failed with code %d (line %d): %s\n",    \
              #call, rc, __LINE__, cudaGetErrorString(rc));                    \
      SMESH_ERROR("fatal cuda error");                                         \
    }                                                                          \
  }
#endif


} // namespace smesh

#endif // SMESH_BASE_HPP