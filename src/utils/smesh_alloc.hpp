#ifndef SMESH_ALLOC_HPP
#define SMESH_ALLOC_HPP

#include "smesh_config.hpp"

namespace smesh {
namespace internal_ {
void *track_malloc(size_t nbytes, const char *file, int line);
void track_free(void *ptr, const char *file, int line);
void *track_realloc(void *ptr, size_t nbytes, const char *file, int line);
void *track_calloc(size_t n, size_t size, const char *file, int line);
} // namespace internal_
} // namespace smesh

#ifndef SMESH_ENABLE_MEM_DIAGNOSTICS

#define SMESH_ALLOC(_nbytes) malloc(_nbytes)
#define SMESH_FREE(_ptr) free(_ptr)
#define SMESH_REALLOC(_ptr, _nbytes) realloc(_ptr, _nbytes)
#define SMESH_CALLOC(_n, _size) calloc(_n, _size)

#else // SMESH_ENABLE_MEM_DIAGNOSTICS

#define SMESH_ALLOC(_nbytes)                                                   \
  smesh::internal_::track_malloc(_nbytes, __FILE__, __LINE__)
#define SMESH_FREE(_ptr) smesh::internal_::track_free(_ptr, __FILE__, __LINE__)
#define SMESH_REALLOC(_ptr, _nbytes)                                           \
  smesh::internal_::track_realloc(_ptr, _nbytes, __FILE__, __LINE__)
#define SMESH_CALLOC(_n, _size)                                                \
  smesh::internal_::track_calloc(_n, _size, __FILE__, __LINE__)
#endif // SMESH_ENABLE_MEM_DIAGNOSTICS

#endif // SMESH_ALLOC_HPP