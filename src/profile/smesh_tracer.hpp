#ifndef SMESH_TRACER_HPP
#define SMESH_TRACER_HPP

#include <stdio.h>

#include "smesh_base.hpp"

namespace smesh {

class ScopedTimeTracer {
public:
    const char *name;
    double start_time;
    SMESH_INLINE ScopedTimeTracer(const char *name) {
        this->name = name;
        this->start_time = time_seconds();
    }

    SMESH_INLINE ~ScopedTimeTracer() {
        double end_time = time_seconds();
        printf("%s: %g seconds\n", this->name, end_time - this->start_time);
    }
};

#define SMESH_TRACE_SCOPE(name) ScopedTimeTracer _scoped_time_tracer_##__LINE__(name)

class ScopedPerfTracer {
public:
    const char *name;
    double start_time;
    long long int flops;
    
    SMESH_INLINE ScopedPerfTracer(const char *name, long long int flops) {
        this->name = name;
        this->flops = flops;
        this->start_time = time_seconds();
    }
    SMESH_INLINE ~ScopedPerfTracer() {
      double end_time = time_seconds();
      double duration = end_time - this->start_time;
      double throughput = this->flops / duration;
      printf("%s: %g GFlops/s\n", this->name, throughput / 1e9);
    }
};

#define SMESH_TRACE_PERF(name, flops) ScopedPerfTracer _scoped_perf_tracer_##__LINE__(name, flops)

} // namespace smesh

#endif // SMESH_TRACER_HPP