#include "smesh_alloc.hpp"
#include "smesh_base.hpp"
#include "smesh_communicator.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace smesh {
namespace internal_ {

static size_t allocated_bytes = 0;

static const char *path_basename(const char *path) {
  const char *slash = std::strrchr(path, '/');
#if defined(_WIN32)
  const char *bs = std::strrchr(path, '\\');
  if (bs && (!slash || bs > slash)) {
    slash = bs;
  }
#endif
  return slash ? slash + 1 : path;
}

static void format_byte_size(char *buf, size_t buf_len, size_t bytes) {
  if (bytes < 1024) {
    std::snprintf(buf, buf_len, "%zu B", bytes);
    return;
  }
  static const char *const k_units[] = {"KiB", "MiB", "GiB", "TiB"};
  double v = static_cast<double>(bytes);
  unsigned divs = 0;
  while (v >= 1024.0 && divs < 4) {
    v /= 1024.0;
    ++divs;
  }
  std::snprintf(buf, buf_len, "%.2f %s (%zu B)", v, k_units[divs - 1], bytes);
}

static void print_memory_usage(const char *file, int line) {
  if (Communicator::world()->rank() != 0) {
    return;
  }

  int SMESH_MEMORY_USAGE_LOG_LEVEL = 0;
  SMESH_READ_ENV(SMESH_MEMORY_USAGE_LOG_LEVEL, atoi);

  if (!SMESH_MEMORY_USAGE_LOG_LEVEL)
    return;

  char human[80];
  format_byte_size(human, sizeof(human), allocated_bytes);
  const char *base = path_basename(file);
  std::printf("---- smesh memory (rank 0) ----------\n"
              "  total allocated :  %s\n"
              "  last touch      :  %s:%d\n"
              "-------------------------------------\n",
              human, base, line);
  std::fflush(stdout);
}

void *track_malloc(size_t nbytes, const char *file, int line) {
  void *ptr = std::malloc(nbytes);

  if (!ptr) {
    SMESH_ERROR("Failed to allocate %zu bytes (%zu total allocated) at %s:%d\n",
                nbytes, allocated_bytes, file, line);
  }

  allocated_bytes += nbytes;
  print_memory_usage(file, line);
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
  print_memory_usage(file, line);
  return ptr;
}

void *track_calloc(size_t n, size_t size, const char *file, int line) {
  void *ptr = std::calloc(n, size);
  if (!ptr) {
    SMESH_ERROR("Failed to calloc %zu bytes (%zu total allocated) at %s:%d\n",
                n * size, allocated_bytes, file, line);
  }

  allocated_bytes += n * size;
  print_memory_usage(file, line);
  return ptr;
}

void track_external_alloc(void *ptr, size_t nbytes, const char *file,
                          int line) {
  if (!ptr) {
    SMESH_ERROR("Failed to allocate %zu bytes (%zu total allocated) at %s:%d\n",
                nbytes, allocated_bytes, file, line);
  }
  allocated_bytes += nbytes;
  print_memory_usage(file, line);
}

} // namespace internal_
} // namespace smesh
