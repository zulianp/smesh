cmake_minimum_required(VERSION 3.20)

project(smesh VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/SMESHDependencies.cmake)
include(cmake/SMESHOptions.cmake)

if(SMESH_ENABLE_PYTHON)
    include(cmake/SMESHPython.cmake)
endif()


set(SMESH_LIBRARY_SCAN_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    CACHE STRING "Folders to scan for smesh library sources/headers (semicolon-separated)"
)

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

add_library(smesh::smesh ALIAS smesh)

set_target_properties(smesh PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(smesh
    PUBLIC
        ${_smesh_library_build_include_dirs}
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_executable(mesh_to_packed src/drivers/mesh_to_packed.exe.cpp)
target_link_libraries(mesh_to_packed PRIVATE smesh)

install(
    TARGETS smesh mesh_to_packed
    EXPORT smeshTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    FILES src/smesh.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT smeshTargets
    NAMESPACE smesh::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smesh
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
