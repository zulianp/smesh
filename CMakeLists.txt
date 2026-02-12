cmake_minimum_required(VERSION 3.20)

project(smesh VERSION 0.1.0 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(SMESH_PUBLIC_SUBMODULES "")
include(cmake/SMESHOptions.cmake)
include(cmake/SMESHDependencies.cmake)
include(cmake/SMESHSanitizer.cmake)

if(SMESH_ENABLE_PYTHON)
    include(cmake/SMESHPython.cmake)
endif()

set(SMESH_LIBRARY_SCAN_DIRS
    "src/base"
    "src/drivers"
    "src/frontend"
    "src/graph"
    "src/io"
    "src/mesh"
    "src/packedmesh"
    "src/profile"
    "src/quadrature"
    "src/shape"
    "src/sorting"
    "src/ssmesh"
    "src/utils"
)

if(SMESH_ENABLE_MPI)
    list(APPEND SMESH_LIBRARY_SCAN_DIRS "src/distributed")
endif()

set(_smesh_library_files "")
set(_smesh_library_include_dirs "")
foreach(_dir IN LISTS SMESH_LIBRARY_SCAN_DIRS)
    if(NOT IS_ABSOLUTE "${_dir}")
        set(_dir "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}")
    endif()

    file(GLOB_RECURSE _dir_files
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        "${_dir}/*.c"
        "${_dir}/*.cc"
        "${_dir}/*.cpp"
        "${_dir}/*.cxx"
        "${_dir}/*.h"
        "${_dir}/*.hh"
        "${_dir}/*.hpp"
        "${_dir}/*.hxx"
        "${_dir}/*.inl"
    )

    list(APPEND _smesh_library_include_dirs "${_dir}")
    foreach(_file IN LISTS _dir_files)
        if(_file MATCHES "\\.(h|hh|hpp|hxx|inl)$")
            get_filename_component(_file_dir "${_file}" DIRECTORY)
            list(APPEND _smesh_library_include_dirs "${_file_dir}")
        endif()
    endforeach()
    list(APPEND _smesh_library_files ${_dir_files})
endforeach()

list(REMOVE_DUPLICATES _smesh_library_include_dirs)
set(_smesh_library_build_include_dirs "")
foreach(_inc_dir IN LISTS _smesh_library_include_dirs)
    list(APPEND _smesh_library_build_include_dirs "$<BUILD_INTERFACE:${_inc_dir}>")
endforeach()


list(FILTER _smesh_library_files EXCLUDE REGEX "\\.exe\\.(c|cc|cpp|cxx)$")
list(REMOVE_DUPLICATES _smesh_library_files)

add_library(smesh ${_smesh_library_files})
target_link_libraries(smesh
    PUBLIC
        ${SMESH_PUBLIC_SUBMODULES}
        ${SMESH_DEP_LIBRARIES}
    PRIVATE
        ${SMESH_SUBMODULES}
)

# Dependency include paths (kept PRIVATE unless they leak through public headers).
if(SMESH_DEP_INCLUDES OR SMESH_BUILD_INCLUDES)
    target_include_directories(smesh PRIVATE ${SMESH_DEP_INCLUDES} ${SMESH_BUILD_INCLUDES})
endif()

add_library(smesh::smesh ALIAS smesh)

set_target_properties(smesh PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
    POSITION_INDEPENDENT_CODE ON
)

# Ensure any target that links `smesh` (and thus includes smesh headers) is
# compiled with C++17. This fixes build failures on toolchains defaulting to
# C++14 (e.g. `if constexpr` requires C++17).
target_compile_features(smesh PUBLIC cxx_std_17)

target_include_directories(smesh
    PUBLIC
        ${_smesh_library_build_include_dirs}
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(_smesh_driver_scan_dir "${CMAKE_CURRENT_SOURCE_DIR}/src/drivers")
file(GLOB_RECURSE _smesh_driver_sources
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS
    "${_smesh_driver_scan_dir}/*.exe.cpp"
)

set(_smesh_driver_targets "")
foreach(_driver_source IN LISTS _smesh_driver_sources)
    file(RELATIVE_PATH _driver_rel "${_smesh_driver_scan_dir}" "${_driver_source}")

    get_filename_component(_driver_target "${_driver_rel}" NAME_WE) # => foo.exe
    string(REGEX REPLACE "\\.exe$" "" _driver_target "${_driver_target}")

    if("${_driver_target}" IN_LIST _smesh_driver_targets)
        string(REGEX REPLACE "\\.exe\\.cpp$" "" _driver_target "${_driver_rel}")
        string(REPLACE "/" "_" _driver_target "${_driver_target}")
        string(REPLACE "\\" "_" _driver_target "${_driver_target}")
    endif()

    add_executable("${_driver_target}" "${_driver_source}")
    target_link_libraries("${_driver_target}" PRIVATE smesh)
    list(APPEND _smesh_driver_targets "${_driver_target}")
endforeach()

install(
    TARGETS smesh ${_smesh_driver_targets}
    EXPORT smeshTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)


install(
    EXPORT smeshTargets
    NAMESPACE smesh::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smesh
)

set(_smesh_config_hpp_in "${CMAKE_CURRENT_SOURCE_DIR}/src/base/smesh_config.hpp.in")
set(_smesh_config_hpp_out "${CMAKE_CURRENT_BINARY_DIR}/smesh_config.hpp")
configure_file("${_smesh_config_hpp_in}" "${_smesh_config_hpp_out}" @ONLY)

target_include_directories(smesh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

install(
    FILES "${_smesh_config_hpp_out}"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/smeshConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/smeshConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smesh
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/smeshConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/smeshConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/smeshConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smesh
)
