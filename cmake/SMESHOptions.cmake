# SFEMOptions

# ##############################################################################
option(BUILD_SHARED_LIBS "build shared libraries" OFF)
option(SMESH_ENABLE_AVX2 "Enable AVX2 intrinsics" OFF)
option(SMESH_ENABLE_AVX512 "Enable AVX512 intrinsics" OFF)
option(SMESH_ENABLE_AVX512_SORT "Enable AVX512 sort" OFF)
option(SMESH_ENABLE_BLAS "Enable BLAS" OFF)
option(SMESH_ENABLE_CUDA "Enable CUDA support" OFF)
option(SMESH_ENABLE_CUDA_LINEINFO "Enable cuda line info for profiling" OFF)
option(SMESH_ENABLE_DEV_MODE "Add additional flags for more strict compilation" ON)
option(SMESH_ENABLE_GLIBCXX_DEBUG "uses flags -D_GLIBCXX_DEBUG when compiling in debug mode" OFF)
option(SMESH_ENABLE_LAPACK "Enable Lapck support" OFF)
option(SMESH_ENABLE_METIS "Enable METIS graph-partitioning" OFF)
option(SMESH_ENABLE_MPI "Enable MPI support" ON)
option(SMESH_ENABLE_OPENMP "Enable OpenMP support" OFF)
# option(SMESH_ENABLE_PYTHON "Enable python bindings for SFEM" OFF)
option(SMESH_ENABLE_RYAML "Enable YAML input files with RapidYAML" OFF)
# option(SMESH_ENABLE_AGGRESSIVE_OPT "Enable aggressive optimizations" OFF)

option(SMESH_USE_OCCUPANCY_MAX_POTENTIAL "Enable usage of cudaOccupancyMaxPotentialBlockSize" OFF)
option(SMESH_ENABLE_TRACE "Eneable trace facilities and output sfem.trace.csv (Override with SMESH_TRACE_FILE in the env)" ON)

get_directory_property(HAS_PARENT PARENT_DIRECTORY)

if(HAS_PARENT)
    option(SMESH_ENABLE_TESTING "Build the tests" OFF)
    option(SMESH_ENABLE_BENCHMARK "enable benchmark suite" OFF)
else()
    option(SMESH_ENABLE_TESTING "Build the tests" ON)
    option(SMESH_ENABLE_BENCHMARK "enable benchmark suite" OFF)
endif()

# ##############################################################################
# Handle xSDK defaults
# ##############################################################################

include(CMakeDependentOption)

if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE
            "Release"
            CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

        message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

# ##############################################################################
# ##############################################################################
# ##############################################################################

if(SMESH_ENABLE_DEV_MODE)
    set(SMESH_DEV_FLAGS
        "-Wall -Wextra -pedantic -Werror -Werror=enum-compare -Werror=delete-non-virtual-dtor -Werror=reorder -Werror=return-type" # -Werror=uninitialized
    )
endif()

if(SMESH_ENABLE_GLIBCXX_DEBUG)
    set(SMESH_SPECIAL_DEBUG_FLAGS
        "${SMESH_SPECIAL_DEBUG_FLAGS} -D_GLIBCXX_DEBUG")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SMESH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} ${SMESH_SPECIAL_DEBUG_FLAGS}")

if(SMESH_ENABLE_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=core-avx2 -DSMESH_ENABLE_AVX2_SORT")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx2 -DSMESH_ENABLE_AVX2_SORT")
endif()




