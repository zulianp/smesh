#ifndef SMESH_TRACER_HPP
#define SMESH_TRACER_HPP

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

} // namespace smesh

#endif // SMESH_TRACER_HPP