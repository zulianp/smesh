#include "smesh_tracer.hpp"

#include "smesh_base.hpp"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <map>

#ifdef SMESH_ENABLE_CUDA
#include <cuda_runtime.h>
#include <nvToolsExt.h>
#endif

// #define SMESH_ENABLE_BLOCK_KERNELS

namespace smesh {
class Tracer::Impl {
public:
  std::map<std::string, std::pair<int, double>> events;
  int rank{0};
  bool log_mode{false};

  void dump() {
    if (rank)
      return;

    const char *SMESH_TRACE_FILE = "smesh.trace.csv";
    SMESH_READ_ENV(SMESH_TRACE_FILE, );

    std::ofstream os(SMESH_TRACE_FILE);

    if (!os.good()) {
      SMESH_ERROR("Unable to write trace file!\n");
    }

    os << "name,calls,total,avg\n";
    for (auto &e : events) {
      os << e.first << "," << e.second.first << "," << e.second.second << ","
         << e.second.second / e.second.first << "\n";
    }

    os << std::flush;
    os.close();
  }

  void record_event(const std::string &name, const double duration) {
    auto &e = events[name];
    e.first++;
    e.second += duration;

    if (log_mode) {
#ifdef SMESH_ENABLE_CUDA
      size_t free, total;
      cudaMemGetInfo(&free, &total);

      if (!rank) {
        printf("-- LOG[%d]: %s (%g)\n"
               "   MEMORY: free %g [GB] (total %g [GB])\n",
               rank, name.c_str(), duration, free * 1e-9, total * 1e-9);
      }

#else
      if (!rank) {
        printf("-- LOG: %s (%g)\n", name.c_str(), duration);
      }
#endif
      fflush(stdout);
    }
  }

  Impl() {}
};

Tracer &Tracer::instance() {
  static Tracer instance_;
  return instance_;
}

void Tracer::set_rank(const int rank) { impl_->rank = rank; }

void Tracer::record_event(const char *name, const double duration) {
  impl_->record_event(name, duration);
}

void Tracer::record_event(std::string &&name, const double duration) {
  impl_->record_event(name, duration);
}

Tracer::Tracer() : impl_(std::make_unique<Impl>()) {
  int SMESH_ENABLE_LOG = 0;
  SMESH_READ_ENV(SMESH_ENABLE_LOG, atoi);
  impl_->log_mode = SMESH_ENABLE_LOG;
}

Tracer::~Tracer() {
  impl_->dump();
  impl_ = nullptr;
}

ScopedEvent::ScopedEvent(const char *format, int num) {
  char str[1024];
  int err = snprintf(str, 1024, format, num);
  if (err < 0)
    SMESH_ERROR("UNABLE TO TRACE %s\n", format);

  name = str;

#ifdef SMESH_ENABLE_BLOCK_KERNELS
  smesh::device_synchronize();
#endif
#ifdef SMESH_ENABLE_CUDA
  nvtxRangePushA(name.c_str());
#endif
  elapsed = time_seconds();
}

ScopedEvent::ScopedEvent(const char *name) : name(name) {
#ifdef SMESH_ENABLE_BLOCK_KERNELS
  smesh::device_synchronize();
#endif
#ifdef SMESH_ENABLE_CUDA
  nvtxRangePushA(this->name.c_str());
#endif
  elapsed = time_seconds();
}

ScopedEvent::~ScopedEvent() {
#ifdef SMESH_ENABLE_BLOCK_KERNELS
  smesh::device_synchronize();
#endif

  elapsed = time_seconds() - elapsed;

#ifdef SMESH_ENABLE_CUDA
  nvtxRangePop();
#endif
  Tracer::instance().record_event(std::move(name), elapsed);
}

} // namespace smesh
