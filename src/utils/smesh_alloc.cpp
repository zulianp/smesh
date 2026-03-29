#include "smesh_alloc.hpp"

#include "smesh_base.hpp"

#include <cstdlib>

namespace smesh {
namespace internal_ {

static size_t allocated_bytes = 0;

void *track_malloc(size_t nbytes, const char *file, int line) {
  void *ptr = std::malloc(nbytes);

  if (!ptr) {
    SMESH_ERROR("Failed to allocate %zu bytes (%zu total allocated) at %s:%d\n",
                nbytes, allocated_bytes, file, line);
  }

  allocated_bytes += nbytes;
  return ptr;
}

void track_free(void *ptr, const char *file, int line) {
  SMESH_UNUSED(file);
  SMESH_UNUSED(line);
  std::free(ptr);
}

void *track_realloc(void *ptr, size_t nbytes, const char *file, int line) {
  ptr = std::realloc(ptr, nbytes);

  if (!ptr) {
    SMESH_ERROR(
        "Failed to reallocate %zu bytes (%zu total allocated) at %s:%d\n",
        nbytes, allocated_bytes, file, line);
  }

  allocated_bytes += nbytes; // TODO FIXME: this is not correct
  return ptr;
}

void *track_calloc(size_t n, size_t size, const char *file, int line) {
  void *ptr = std::calloc(n, size);
  if (!ptr) {
    SMESH_ERROR("Failed to calloc %zu bytes (%zu total allocated) at %s:%d\n",
                n * size, allocated_bytes, file, line);
  }

  allocated_bytes += n * size;
  return ptr;
}

} // namespace internal_
} // namespace smesh
