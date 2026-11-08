#ifndef SMESH_TRACER_HPP
#define SMESH_TRACER_HPP

#include "smesh_base.hpp"

#include <stdio.h>
#include <string>
#include <memory>

namespace smesh {

    class Tracer {
    public:
        static Tracer &instance();
        void           record_event(const char *name, const double duration);
        void           record_event(std::string &&name, const double duration);
        Tracer();
        ~Tracer();

        void set_rank(const int rank);

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    class ScopedEvent {
    public:
        // const char *name{nullptr};
        std::string name;
        double      elapsed{0};
        bool        active{true};  // false inside a suppressed concurrent scope
        ScopedEvent(const char *name);
        ScopedEvent(const char *format, int num);
        ~ScopedEvent();
    };

    /// A scope run concurrently with the rest of the program (a thread of its own), traced only
    /// with SMESH_ENABLE_TRACE_MUTEX (SMESH_TRACE_SCOPE_MUTEX): without it, the scope and every
    /// scope it opens on its thread record nothing, so no event needs a mutex.
    class ScopedTraceSuppression {
    public:
        ScopedTraceSuppression();
        ~ScopedTraceSuppression();
        ScopedTraceSuppression(const ScopedTraceSuppression &)            = delete;
        ScopedTraceSuppression &operator=(const ScopedTraceSuppression &) = delete;
    };

    class ScopedTimeTracer {
    public:
        const char  *name;
        double       start_time;
        SMESH_INLINE ScopedTimeTracer(const char *name) {
            this->name       = name;
            this->start_time = time_seconds();
        }

        SMESH_INLINE ~ScopedTimeTracer() {
            double end_time = time_seconds();
            printf("%s: %g [s]\n", this->name, end_time - this->start_time);
        }
    };



// #define SMESH_TRACE_SCOPE(name) ScopedTimeTracer _scoped_time_tracer_##__LINE__(name)

    class ScopedPerfTracer {
    public:
        const char   *name;
        double        start_time;
        long long int flops;

        SMESH_INLINE ScopedPerfTracer(const char *name, long long int flops) {
            this->name       = name;
            this->flops      = flops;
            this->start_time = time_seconds();
        }
        SMESH_INLINE ~ScopedPerfTracer() {
            double end_time   = time_seconds();
            double duration   = end_time - this->start_time;
            double throughput = this->flops / duration;
            printf("%s: %g [GFlops/s]\n", this->name, throughput / 1e9);
        }
    };
    

#ifdef SMESH_ENABLE_TRACE
#define SMESH_TRACE_SCOPE(name) smesh::ScopedEvent smesh_scoped_trace_event_(name);
#define SMESH_TRACE_SCOPE_VARIANT(format__, num__) smesh::ScopedEvent smesh_scoped_trace_event_(format__, num__);
#define SMESH_TRACE_PERF(name, flops) ScopedPerfTracer _scoped_perf_tracer_##__LINE__(name, flops)
// A scope that runs on a thread of its own: traced, every event under a mutex, with
// SMESH_ENABLE_TRACE_MUTEX; otherwise it and the scopes it opens are not traced.
#ifdef SMESH_ENABLE_TRACE_MUTEX
#define SMESH_TRACE_SCOPE_MUTEX(name) smesh::ScopedEvent smesh_scoped_trace_event_(name);
#else
#define SMESH_TRACE_SCOPE_MUTEX(name) smesh::ScopedTraceSuppression smesh_scoped_trace_suppression_;
#endif
#else
#define SMESH_TRACE_SCOPE(...)
#define SMESH_TRACE_SCOPE_VARIANT(...)
#define SMESH_TRACE_PERF(...)
#define SMESH_TRACE_SCOPE_MUTEX(...)
#endif

}  // namespace smesh

#endif  // SMESH_TRACER_HPP