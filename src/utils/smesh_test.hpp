#ifndef SMESH_TEST_H
#define SMESH_TEST_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <exception>
#include <iostream>


#include "smesh_base.hpp"
#include "smesh_context.hpp"
#include "smesh_tracer.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SMESH_ENABLE_MPI
#include <mpi.h>
#endif

#define SMESH_TEST_SUCCESS 0
#define SMESH_TEST_FAILURE 1
#define SMESH_TEST_SKIPPED 2

static void smesh_print_test_info() {
    printf("=======================\n");
    printf("SMESH_TESTING Info:\n");

    int size = 1, rank = 0;

#ifdef SMESH_ENABLE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (!rank) {
        printf("Comm_size=%d\n", size);
    }

#ifdef _OPENMP
#pragma omp parallel
    {
        size_t start, len;
        int    id  = omp_get_thread_num();
        int    num = omp_get_num_threads();

        if (!id) {
            printf("OMP_NUM_THREADS=%d\n", num);
        }
    }
#endif
    printf("=======================\n");
}

#define SMESH_UNIT_TEST_INIT(argc, argv)  \
    smesh::Context context__(argc, argv); \
    smesh_print_test_info();              \
    int err = 0;

#define SMESH_RUN_TEST(test_)                                                                                    \
    do {                                                                                                        \
        SMESH_TRACE_SCOPE(#test_);                                                                               \
        int this_test = 0;                                                                                      \
        try {                                                                                                   \
            if ((this_test = test_())) fprintf(stderr, "TEST: %s failed! %s:%d\n", #test_, __FILE__, __LINE__); \
        } catch (const std::exception &ex) {                                                                    \
            fprintf(stderr, "Exception: %s, in test %s! %s:%d\n", ex.what(), #test_, __FILE__, __LINE__);       \
            this_test = SMESH_TEST_FAILURE;                                                                      \
        }                                                                                                       \
        err += this_test;                                                                                       \
    } while (0)

#define SMESH_UNIT_TEST_FINALIZE()                                                     \
    do {                                                                              \
        if (err) fprintf(stderr, "The number of tests failed in unit is %d!\n", err); \
    } while (0)
#define SMESH_UNIT_TEST_ERR() (err)

static inline int smesh_test_assert(const int expr, const char *expr_string, const char *file, const int line) {
    if (!expr) {
        fprintf(stderr, "\nAssertion failure: %s\nAt %s:%d\n\n", expr_string, file, line);
        assert(0);
        return SMESH_TEST_FAILURE;
    }

    return 0;
}

static inline int smesh_test_assert_approxeq(double a, double b, double tol, const char *file, const int line) {
    if (fabs(a - b) > tol) {
        fprintf(stderr, "\nAssertion failure: %g != %g (diff %g > %g) \nAt %s:%d\n\n", a, b, fabs(a - b), tol, file, line);
        assert(0);
        return SMESH_TEST_FAILURE;
    }

    return 0;
}

static inline int smesh_assert_array_approx_eq(const ptrdiff_t i,
                                              double          a,
                                              double          b,
                                              double          tol,
                                              const char     *file,
                                              const int       line) {
    if (fabs(a - b) > tol) {
        fprintf(stderr,
                "\nAssertion failure at entry %ld: %g != %g (diff %g > %g) \nAt %s:%d\n\n",
                i,
                a,
                b,
                fabs(a - b),
                tol,
                file,
                line);
        assert(0);
        return SMESH_TEST_FAILURE;
    }

    return 0;
}

#define SMESH_TEST_ASSERT(expr)                                                    \
    if (smesh_test_assert(expr, #expr, __FILE__, __LINE__) == SMESH_TEST_FAILURE) { \
        return SMESH_TEST_FAILURE;                                                 \
    }

#define SMESH_TEST_APPROXEQ(a, b, tol)                                                    \
    if (smesh_test_assert_approxeq(a, b, tol, __FILE__, __LINE__) == SMESH_TEST_FAILURE) { \
        return SMESH_TEST_FAILURE;                                                        \
    }

#define SMESH_ASSERT_ARRAY_APPROX_EQ(n__, a__, b__, tol__)                                                               \
    do {                                                                                                                \
        ptrdiff_t nfails = 0;                                                                                           \
        for (ptrdiff_t i__ = 0; i__ < (n__); i__++) {                                                                     \
            if (smesh_assert_array_approx_eq(i__, (a__)[i__], (b__)[i__], tol__, __FILE__, __LINE__) == SMESH_TEST_FAILURE) { \
                nfails++;                                                                                                \
            }                                                                                                           \
        }                                                                                                               \
        if (nfails) return SMESH_TEST_FAILURE;                                                                           \
    } while (0)


template <typename T>
int smesh_test_assert_equal(T a, T b, const char *file, const int line) {
    if (a != b) {
        std::cerr << "\nAssertion failure: " << a << "!= " << b << "  \nAt" << file << ":" << line << "\n";
        assert(0);
        return SMESH_TEST_FAILURE;
    }

    return 0;
}

#define SMESH_TEST_EQ(a, b)                                                       \
    if (smesh_test_assert_equal(a, b, __FILE__, __LINE__) == SMESH_TEST_FAILURE) { \
        return SMESH_TEST_FAILURE;                                                \
    }

#define SMESH_ASSERT_ARRAY_EQ(n__, a__, b__)     \
    for (ptrdiff_t i__ = 0; i__ < n__; i__++) { \
        SMESH_TEST_EQ(a__[i__], b__[i__]);       \
    }

#endif  // SMESH_TEST
