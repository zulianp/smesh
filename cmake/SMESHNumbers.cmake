# ##############################################################################
# XSDK_PRECISION
# ##############################################################################

option(SMESH_ENABLE_CUSTOM_NUMBERS "Enable custom numbers support. (this option will be removed in the future)" ON)

set(SMESH_REAL_TYPE          "float64"           CACHE STRING "Set SMESH real_t type. Used for solution vectors.")
set(SMESH_SCALAR_TYPE        ${SMESH_REAL_TYPE}   CACHE STRING "Set SMESH scalar_t type. Used for local kernel computations")
set(SMESH_GEOM_TYPE          "float32"           CACHE STRING "Set SMESH geom_t type")
set(SMESH_JACOBIAN_CPU_TYPE  "float32"           CACHE STRING "Set SMESH jacobian_t type")
set(SMESH_JACOBIAN_GPU_TYPE  "float"             CACHE STRING "Set SMESH cu_jacobian_t type (goal is half)")
set(SMESH_ACCUMULATOR_TYPE   ${SMESH_REAL_TYPE}   CACHE STRING "Set SMESH accumulator_t type")
set(SMESH_IDX_TYPE           "int32"             CACHE STRING "Set SMESH idx_t type")
set(SMESH_COUNT_TYPE         ${SMESH_IDX_TYPE}    CACHE STRING "Set SMESH count_t type")
set(SMESH_ELEMENT_IDX_TYPE   ${SMESH_IDX_TYPE}    CACHE STRING "Set SMESH element_idx_t type")
set(SMESH_LOCAL_IDX_TYPE     "int16"             CACHE STRING "Set SMESH local_idx_t type")

if(SMESH_REAL_TYPE STREQUAL "float64" )
    set(SMESH_SIZEOF_REAL_T 8)
elseif(SMESH_REAL_TYPE STREQUAL "float32" )
    set(SMESH_SIZEOF_REAL_T 4)
else()
    message(FATAL_ERROR "Not real number type for `${SMESH_REAL_TYPE}`!")
endif()

function(sfem_simd_vector_size type vec_size)
    if("${type}" STREQUAL "float32")
        set(${vec_size} 8 PARENT_SCOPE)
    elseif("${type}" STREQUAL "float64")
        set(${vec_size} 4 PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Not simd vector size type for `${type}`!")
    endif()
endfunction()

sfem_simd_vector_size(${SMESH_SCALAR_TYPE} SMESH_VEC_SIZE_DEFAUT)
set(SMESH_VEC_SIZE  ${SMESH_VEC_SIZE_DEFAUT}   CACHE STRING "Set number of simd lanes for vec_t type")

function(sfem_num_bytes type num_bytes)
    if("${type}" STREQUAL "int16")
        set(${num_bytes} 2 PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint16")
        set(${num_bytes} 2 PARENT_SCOPE)
    elseif("${type}" STREQUAL "int32")
        set(${num_bytes} 4 PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint32")
        set(${num_bytes} 4 PARENT_SCOPE)
    elseif("${type}" STREQUAL "int64")
        set(${num_bytes} 8 PARENT_SCOPE)
    elseif("${type}" STREQUAL "float16")
        set(${num_bytes} 2 PARENT_SCOPE)
    elseif("${type}" STREQUAL "float32")
        set(${num_bytes} 4 PARENT_SCOPE)
    elseif("${type}" STREQUAL "float64")
        set(${num_bytes} 8 PARENT_SCOPE)
    else()
        message(FATAL_ERROR "No number of bytes for `${type}`! (extend sfem_num_bytes function!)")
    endif()
endfunction()

function(sfem_c_type type c_type)
    if("${type}" STREQUAL "int16")
        set(${c_type} "int16_t " PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint16")
            set(${c_type} "uint16_t " PARENT_SCOPE)
    elseif("${type}" STREQUAL "int32")
        set(${c_type} "int32_t " PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint32")
        set(${c_type} "uint32_t " PARENT_SCOPE)
    elseif("${type}" STREQUAL "int64")
        set(${c_type} "int64_t " PARENT_SCOPE)
    # elseif("${type}" STREQUAL "float16")
    #     set(${c_type} "NO direct MPI Support" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float32")
        set(${c_type} "float" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float64")
        set(${c_type} "double" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Not C type for `${type}`!")
    endif()
endfunction()

function(sfem_mpi_type type mpi_type)
    if("${type}" STREQUAL "int16")
        set(${mpi_type} "MPI_INT16_T " PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint16")
        set(${mpi_type} "MPI_UINT16_T " PARENT_SCOPE)
    elseif("${type}" STREQUAL "int32")
        set(${mpi_type} "MPI_INT32_T " PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint32")
        set(${mpi_type} "MPI_UINT32_T " PARENT_SCOPE)
    elseif("${type}" STREQUAL "int64")
        set(${mpi_type} "MPI_INT64_T " PARENT_SCOPE)
    # elseif("${type}" STREQUAL "float16")
    #     set(${mpi_type} "NO direct MPI Support" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float32")
        set(${mpi_type} "MPI_FLOAT" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float64")
        set(${mpi_type} "MPI_DOUBLE" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Not MPI type for `${type}`!")
    endif()
endfunction()

function(sfem_cusparse_type type cusparse_type)
    # if("${type}" STREQUAL "int16")
    #     set(${cusparse_type} "int16_t " PARENT_SCOPE)
    # else
    if("${type}" STREQUAL "int32")
        set(${cusparse_type} "CUSPARSE_INDEX_32I" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint32")
        set(${cusparse_type} "UINT32_NOT_SUPPORTED" PARENT_SCOPE)
    elseif("${type}" STREQUAL "int64")
        set(${cusparse_type} "CUSPARSE_INDEX_64I" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float16")
        set(${cusparse_type} "CUDA_R_16F" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float32")
        set(${cusparse_type} "CUDA_R_32F" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float64")
        set(${cusparse_type} "CUDA_R_64F" PARENT_SCOPE)
    elseif("${type}" STREQUAL "int16" AND NOT SMESH_ENABLE_CUDA)
        set(${cusparse_type} "INT16_NOT_SUPPORTED" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint16" AND NOT SMESH_ENABLE_CUDA)
        set(${cusparse_type} "UINT16_NOT_SUPPORTED" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Not CuSparse type for `${type}`!")
    endif()
endfunction()

function(sfem_print_type type print_type)
    if("${type}" STREQUAL "int16")
        set(${print_type} "hd" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint16")
        set(${print_type} "hu" PARENT_SCOPE)
    elseif("${type}" STREQUAL "int32")
        set(${print_type} "d" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint32")
        set(${print_type} "u" PARENT_SCOPE)
    elseif("${type}" STREQUAL "int64")
        set(${print_type} "ld" PARENT_SCOPE)
    # elseif("${type}" STREQUAL "float16")
    #     set(${print_type} "NO direct MPI Support" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float32")
        set(${print_type} "f" PARENT_SCOPE)
    elseif("${type}" STREQUAL "float64")
        set(${print_type} "g" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "No print type for `${type}`!")
    endif()
endfunction()

function(sfem_invalid_index type invalid_index)
    if("${type}" STREQUAL "int16")
        set(${invalid_index} "-1" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint16")
        set(${invalid_index} "((uint16_t)65535U)" PARENT_SCOPE)
    elseif("${type}" STREQUAL "int32")
        set(${invalid_index} "-1" PARENT_SCOPE)
    elseif("${type}" STREQUAL "int64")
        set(${invalid_index} "-1" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint32")
        set(${invalid_index} "4294967295" PARENT_SCOPE)
    elseif("${type}" STREQUAL "uint64")
        set(${invalid_index} "18446744073709551615" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "No invalid type for `${type}`!")
    endif()
endfunction()

# C types
sfem_num_bytes(${SMESH_REAL_TYPE} SMESH_REAL_SIZE)
sfem_num_bytes(${SMESH_SCALAR_TYPE} SMESH_SCALAR_SIZE)
sfem_num_bytes(${SMESH_GEOM_TYPE} SMESH_GEOM_SIZE)
sfem_num_bytes(${SMESH_JACOBIAN_CPU_TYPE} SMESH_JACOBIAN_CPU_SIZE)
sfem_num_bytes(${SMESH_ACCUMULATOR_TYPE} SMESH_ACCUMULATOR_SIZE)
sfem_num_bytes(${SMESH_IDX_TYPE} SMESH_IDX_SIZE)
sfem_num_bytes(${SMESH_COUNT_TYPE} SMESH_COUNT_SIZE)
sfem_num_bytes(${SMESH_ELEMENT_IDX_TYPE} SMESH_ELEMENT_IDX_SIZE)
sfem_num_bytes(${SMESH_LOCAL_IDX_TYPE} SMESH_LOCAL_IDX_SIZE)

# Num Bytes
sfem_c_type(${SMESH_REAL_TYPE} SMESH_REAL_C_TYPE)
sfem_c_type(${SMESH_SCALAR_TYPE} SMESH_SCALAR_C_TYPE)
sfem_c_type(${SMESH_GEOM_TYPE} SMESH_GEOM_C_TYPE)
sfem_c_type(${SMESH_JACOBIAN_CPU_TYPE} SMESH_JACOBIAN_CPU_C_TYPE)
sfem_c_type(${SMESH_ACCUMULATOR_TYPE} SMESH_ACCUMULATOR_C_TYPE)
sfem_c_type(${SMESH_IDX_TYPE} SMESH_IDX_C_TYPE)
sfem_c_type(${SMESH_COUNT_TYPE} SMESH_COUNT_C_TYPE)
sfem_c_type(${SMESH_ELEMENT_IDX_TYPE} SMESH_ELEMENT_IDX_C_TYPE)
sfem_c_type(${SMESH_LOCAL_IDX_TYPE} SMESH_LOCAL_IDX_C_TYPE)

# CuSparse types
sfem_cusparse_type(${SMESH_REAL_TYPE} SMESH_REAL_CUSPARSE_TYPE)
sfem_cusparse_type(${SMESH_SCALAR_TYPE} SMESH_SCALAR_CUSPARSE_TYPE)
sfem_cusparse_type(${SMESH_GEOM_TYPE} SMESH_GEOM_CUSPARSE_TYPE)
sfem_cusparse_type(${SMESH_JACOBIAN_CPU_TYPE} SMESH_JACOBIAN_CPU_CUSPARSE_TYPE)
sfem_cusparse_type(${SMESH_IDX_TYPE} SMESH_IDX_CUSPARSE_TYPE)
sfem_cusparse_type(${SMESH_COUNT_TYPE} SMESH_COUNT_CUSPARSE_TYPE)
sfem_cusparse_type(${SMESH_ELEMENT_IDX_TYPE} SMESH_ELEMENT_IDX_CUSPARSE_TYPE)
# sfem_cusparse_type(${SMESH_LOCAL_IDX_TYPE} SMESH_LOCAL_IDX_CUSPARSE_TYPE)

# MPI types
sfem_mpi_type(${SMESH_REAL_TYPE} SMESH_REAL_MPI_TYPE)
sfem_mpi_type(${SMESH_SCALAR_TYPE} SMESH_SCALAR_MPI_TYPE)
sfem_mpi_type(${SMESH_ACCUMULATOR_TYPE} SMESH_ACCUMULATOR_MPI_TYPE)
sfem_mpi_type(${SMESH_GEOM_TYPE} SMESH_GEOM_MPI_TYPE)
sfem_mpi_type(${SMESH_JACOBIAN_CPU_TYPE} SMESH_JACOBIAN_CPU_MPI_TYPE)
sfem_mpi_type(${SMESH_IDX_TYPE} SMESH_IDX_MPI_TYPE)
sfem_mpi_type(${SMESH_COUNT_TYPE} SMESH_COUNT_MPI_TYPE)
sfem_mpi_type(${SMESH_ELEMENT_IDX_TYPE} SMESH_ELEMENT_IDX_MPI_TYPE)
sfem_mpi_type(${SMESH_LOCAL_IDX_TYPE} SMESH_LOCAL_IDX_MPI_TYPE)

# Print format
sfem_print_type(${SMESH_REAL_TYPE} SMESH_REAL_PRINT_TYPE)
sfem_print_type(${SMESH_SCALAR_TYPE} SMESH_SCALAR_PRINT_TYPE)
sfem_print_type(${SMESH_ACCUMULATOR_TYPE} SMESH_ACCUMULATOR_PRINT_TYPE)
sfem_print_type(${SMESH_GEOM_TYPE} SMESH_GEOM_PRINT_TYPE)
sfem_print_type(${SMESH_JACOBIAN_CPU_TYPE} SMESH_JACOBIAN_CPU_PRINT_TYPE)
sfem_print_type(${SMESH_IDX_TYPE} SMESH_IDX_PRINT_TYPE)
sfem_print_type(${SMESH_COUNT_TYPE} SMESH_COUNT_PRINT_TYPE)
sfem_print_type(${SMESH_ELEMENT_IDX_TYPE} SMESH_ELEMENT_IDX_PRINT_TYPE)
sfem_print_type(${SMESH_LOCAL_IDX_TYPE} SMESH_LOCAL_IDX_PRINT_TYPE)

# Invalid index
sfem_invalid_index(${SMESH_IDX_TYPE} SMESH_IDX_INVALID)
sfem_invalid_index(${SMESH_COUNT_TYPE} SMESH_COUNT_INVALID)
sfem_invalid_index(${SMESH_ELEMENT_IDX_TYPE} SMESH_ELEMENT_IDX_INVALID)
sfem_invalid_index(${SMESH_LOCAL_IDX_TYPE} SMESH_LOCAL_IDX_INVALID)

message(STATUS 
    "--------------------------------------------------------------------------------------\n"
    "\nSMESH Numbers\n"
    "--------------------------------------------------------------------------------------\n"
    "type\t\tC\t\tId\tMPI\t\tprintf\t(CMake option)\n"
    "--------------------------------------------------------------------------------------\n"
    "real_t\t\t${SMESH_REAL_C_TYPE}\t\t${SMESH_REAL_TYPE}\t${SMESH_REAL_MPI_TYPE}\t${SMESH_REAL_PRINT_TYPE}\t(SMESH_REAL_TYPE)\n"
    "scalar_t\t${SMESH_SCALAR_C_TYPE}\t\t${SMESH_SCALAR_TYPE}\t${SMESH_SCALAR_MPI_TYPE}\t${SMESH_SCALAR_PRINT_TYPE}\t(SMESH_SCALAR_TYPE)\n"
    "accumulator_t\t${SMESH_ACCUMULATOR_C_TYPE}\t\t${SMESH_ACCUMULATOR_TYPE}\t${SMESH_ACCUMULATOR_MPI_TYPE}\t${SMESH_ACCUMULATOR_PRINT_TYPE}\t(SMESH_ACCUMULATOR_TYPE)\n"
    "geom_t\t\t${SMESH_GEOM_C_TYPE}\t\t${SMESH_GEOM_TYPE}\t${SMESH_GEOM_MPI_TYPE}\t${SMESH_GEOM_PRINT_TYPE}\t(SMESH_GEOM_TYPE)\n"
    "jacobian_t\t${SMESH_JACOBIAN_CPU_C_TYPE}\t\t${SMESH_JACOBIAN_CPU_TYPE}\t${SMESH_JACOBIAN_CPU_MPI_TYPE}\t${SMESH_JACOBIAN_CPU_PRINT_TYPE}\t(SMESH_JACOBIAN_CPU_TYPE)\n"
    "cu_jacobian_t\t${SMESH_JACOBIAN_GPU_TYPE}\t\t-\t-\t\t-\t(SMESH_JACOBIAN_GPU_TYPE)\n"
    "idx_t\t\t${SMESH_IDX_C_TYPE}\t${SMESH_IDX_TYPE}\t${SMESH_IDX_MPI_TYPE}\t${SMESH_IDX_PRINT_TYPE}\t(SMESH_IDX_TYPE)\n"
    "count_t\t\t${SMESH_COUNT_C_TYPE}\t${SMESH_COUNT_TYPE}\t${SMESH_COUNT_MPI_TYPE}\t${SMESH_COUNT_PRINT_TYPE}\t(SMESH_COUNT_TYPE)\n"
    "element_idx_t\t${SMESH_ELEMENT_IDX_C_TYPE}\t${SMESH_ELEMENT_IDX_TYPE}\t${SMESH_ELEMENT_IDX_MPI_TYPE}\t${SMESH_ELEMENT_IDX_PRINT_TYPE}\t(SMESH_ELEMENT_IDX_TYPE)\n"
    "local_idx_t\t${SMESH_LOCAL_IDX_C_TYPE}\t${SMESH_LOCAL_IDX_TYPE}\t${SMESH_LOCAL_IDX_MPI_TYPE}\t${SMESH_LOCAL_IDX_PRINT_TYPE}\t(SMESH_LOCAL_IDX_TYPE)\n"
    "--------------------------------------------------------------------------------------\n"
    "SMESH_VEC_SIZE=${SMESH_VEC_SIZE}\n"
    "--------------------------------------------------------------------------------------\n"

)

# if(NOT XSDK_PRECISION OR USE_XSDK_DEFAULTS)
#     set(XSDK_PRECISION "DOUBLE")
# endif()

# string(COMPARE EQUAL ${XSDK_PRECISION} "DOUBLE" SMESH_HAVE_DOUBLE_PRECISION)
# string(COMPARE EQUAL ${XSDK_PRECISION} "SINGLE" SMESH_HAVE_SINGLE_PRECISION)
# string(COMPARE EQUAL ${XSDK_PRECISION} "QUAD" SMESH_HAVE_QUAD_PRECISION)

# # ##############################################################################
# # XSDK_INDEX_SIZE
# # ##############################################################################


# if(NOT SMESH_INDEX_BITSIZE)
#     set(SMESH_INDEX_BITSIZE 32 CACHE STRING "Choice of idx_t size between 32 or 64 bits" FORCE)
# endif()

# if(NOT SMESH_COUNT_BITSIZE)
#     set(SMESH_COUNT_BITSIZE 32 CACHE STRING "Choice of count_t size between 32 or 64 bits" FORCE)
# endif()


# if(USE_XSDK_DEFAULTS)
#     set(XSDK_INDEX_SIZE 32)
#     set(SMESH_INDEX_BITSIZE 32)
# else()
#     if(XSDK_INDEX_SIZE)
#         set(SMESH_INDEX_BITSIZE ${XSDK_INDEX_SIZE})
#     elseif(NOT SMESH_INDEX_BITSIZE)
#         set(XSDK_INDEX_SIZE 64)
#         set(SMESH_INDEX_BITSIZE 64)
#     endif()
# endif()