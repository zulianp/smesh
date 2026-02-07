
if(SMESH_ENABLE_RYAML)
    set(RYML_REPO_URL https://github.com/biojppm/rapidyaml CACHE STRING "")
    set(RYML_BRANCH_NAME master CACHE STRING "")
    include(FetchContent)
    FetchContent_Declare(ryml
        GIT_REPOSITORY ${RYML_REPO_URL}
        GIT_TAG ${RYML_BRANCH_NAME}
        GIT_SHALLOW FALSE  # ensure submodules are checked out
    )
    FetchContent_MakeAvailable(ryml)
    list(APPEND SMESH_SUBMODULES ryml::ryml)
endif()

# ##############################################################################
if(WIN32)
    set(GLOB_REPO_URL https://github.com/p-ranav/glob.git CACHE STRING "")
    set(GLOB_BRANCH_NAME master CACHE STRING "")
    include(FetchContent)
    FetchContent_Declare(Glob
        GIT_REPOSITORY ${GLOB_REPO_URL}
        GIT_TAG ${GLOB_BRANCH_NAME}
        GIT_SHALLOW FALSE  # ensure submodules are checked out
    )
    FetchContent_MakeAvailable(Glob)
    list(APPEND SMESH_SUBMODULES Glob)
endif()


# ##############################################################################

if(SMESH_ENABLE_BLAS)
    if(APPLE)
        # Add Accelerate framework for macOS BLAS/LAPACK
        find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
        list(APPEND SMESH_DEP_LIBRARIES ${ACCELERATE_FRAMEWORK})
    else()
        find_package(BLAS REQUIRED)
        list(APPEND SMESH_DEP_LIBRARIES BLAS::BLAS)
    endif()
endif()

# ##############################################################################

if(SMESH_ENABLE_LAPACK)
    find_package(LAPACK REQUIRED)
    list(APPEND SMESH_DEP_LIBRARIES  LAPACK::LAPACK)
endif()

# ##############################################################################

if(SMESH_ENABLE_AVX512_SORT)
	include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/x86-simd-sort/src") 
endif()

# ##############################################################################

if(SMESH_ENABLE_CUDA)
    enable_language(CUDA)
    if(NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()
endif()

# ##############################################################################

if(SMESH_ENABLE_OPENMP)
    if(OPENMP_DIR)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
        list(APPEND SMESH_DEP_INCLUDES "${OPENMP_DIR}/include")
        list(APPEND SMESH_DEP_LIBRARIES "-L${OPENMP_DIR}/lib -lomp" )
    else()
        find_package(OpenMP REQUIRED)
        if(OpenMP_FOUND)
            message(STATUS "OpenMP: ${OpenMP_INCLUDES}")
            set(SMESH_DEP_LIBRARIES "${SMESH_DEP_LIBRARIES};OpenMP::OpenMP_C;OpenMP::OpenMP_CXX")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()
    endif()
endif()

# ##############################################################################

if(SMESH_ENABLE_METIS)
  find_package(METIS REQUIRED)
  if(METIS_FOUND)
    list(APPEND SMESH_BUILD_INCLUDES ${METIS_INCLUDES})
    list(APPEND SMESH_DEP_LIBRARIES ${METIS_LIBRARIES})
    set(SMESH_ENABLE_METIS ON)
  else()
    message(FATAL_ERROR "[Warning] Metis not found")
  endif()
endif()

# ##############################################################################

if(SMESH_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS C CXX)
    list(APPEND SMESH_PUBLIC_SUBMODULES MPI::MPI_CXX MPI::MPI_C)

    set(MATRIXIO_REPO_URL https://github.com/zulianp/matrix.io.git CACHE STRING "")
    set(MATRIXIO_BRANCH_NAME main CACHE STRING "")
    include(FetchContent)
    
    FetchContent_Declare(matrixio
        GIT_REPOSITORY ${MATRIXIO_REPO_URL}
        GIT_TAG ${MATRIXIO_BRANCH_NAME}
        GIT_SHALLOW FALSE  # ensure submodules are checked out
    )
    
    # matrixio is consumed as a build-time dependency; enable its install/export
    # rules so smesh's own export set can be generated without errors.
    set(MATRIXIO_INSTALL ON CACHE BOOL "" FORCE)
    set(MATRIXIO_BUILD_EXTRAS OFF CACHE BOOL "" FORCE)
    set(MATRIXIO_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(matrixio)
    # matrixio exports matrixio::matrixio (alias of target "matrixio")
    list(APPEND SMESH_SUBMODULES matrixio::matrixio)
endif()


# ##############################################################################

find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
    # Generate docs into repo-root api/html (not build/api/html).
    set(SMESH_DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/api")

    # Keep inputs tight: sources + top-level docs (space-separated list).
    set(SMESH_DOXYGEN_INPUT
        "${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/README.md ${PROJECT_SOURCE_DIR}/tetstream_codec.py"
    )

    # Graphviz is optional; enable dot graphs when available.
    find_program(_smesh_dot_exe dot)
    if(_smesh_dot_exe)
        set(SMESH_DOXYGEN_HAVE_DOT YES)
    else()
        set(SMESH_DOXYGEN_HAVE_DOT NO)
    endif()

    set(_smesh_doxyfile_in  "${PROJECT_SOURCE_DIR}/cmake/Doxyfile.in")
    set(_smesh_doxyfile_out "${CMAKE_BINARY_DIR}/Doxyfile.txt")
    configure_file("${_smesh_doxyfile_in}" "${_smesh_doxyfile_out}" @ONLY)

    add_custom_target(docs
        COMMAND "${DOXYGEN_EXECUTABLE}" "${_smesh_doxyfile_out}"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        DEPENDS "${_smesh_doxyfile_out}"
        BYPRODUCTS
            "${SMESH_DOXYGEN_OUTPUT_DIRECTORY}/html/index.html"
        COMMENT "Generating API docs (Doxygen) -> ${SMESH_DOXYGEN_OUTPUT_DIRECTORY}/html"
        VERBATIM
    )

    if(SMESH_ENABLE_INSTALL_DOCS)
        # Pick a reasonable doc install destination even when CMAKE_INSTALL_DOCDIR is empty.
        if("${CMAKE_INSTALL_DOCDIR}" STREQUAL "")
            set(_smesh_install_docdir "${CMAKE_INSTALL_DATAROOTDIR}/doc/${PROJECT_NAME}")
        else()
            set(_smesh_install_docdir "${CMAKE_INSTALL_DOCDIR}")
        endif()

        # Ensure docs exist (and are current) when installing.
        install(CODE
            "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build \"${CMAKE_BINARY_DIR}\" --target docs RESULT_VARIABLE _smesh_docs_res)\n"
            "if(NOT _smesh_docs_res EQUAL 0)\n"
            "  message(FATAL_ERROR \"Failed to build Doxygen docs via target 'docs'.\")\n"
            "endif()\n"
            COMPONENT docs
        )

        # Installs to: <prefix>/${CMAKE_INSTALL_DOCDIR}/api/html/...
        install(
            DIRECTORY "${SMESH_DOXYGEN_OUTPUT_DIRECTORY}/html"
            DESTINATION "${_smesh_install_docdir}/api"
            COMPONENT docs
        )
    endif()
endif()

if(SMESH_ENABLE_INSTALL_DOCS AND NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "SMESH_ENABLE_INSTALL_DOCS requires Doxygen, but it was not found.")
endif()
