#ifndef SMESH_TYPES_HPP
#define SMESH_TYPES_HPP

#include "smesh_base.hpp"
#include <string_view>

namespace smesh {

using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using mask_t = char;
using block_idx_t = u16;

#if defined(__APPLE__)
using f16 = __fp16;
#else
using f16 = _Float16;
#endif

using f32 = float;
using f64 = double;

using real_t = f64;
using geom_t = f32;
using idx_t = i32;
using element_idx_t = i32;
using count_t = i32;

static const f16 F16_MAX = (f16)65504.0f;

template <typename T> struct TypeToString {
  static const std::string_view value() { return "raw"; }
};

template <> struct TypeToString<i16> {
  static const std::string_view value() { return "int16"; }
};

template <> struct TypeToString<i32> {
  static const std::string_view value() { return "int32"; }
};

template <> struct TypeToString<i64> {
  static const std::string_view value() { return "int64"; }
};

template <> struct TypeToString<f16> {
  static const std::string_view value() { return "float16"; }
};

template <> struct TypeToString<f32> {
  static const std::string_view value() { return "float32"; }
};

template <> struct TypeToString<f64> {
  static const std::string_view value() { return "float64"; }
};

enum PrimitiveType {
  SMESH_DEFAULT = 0,
  SMESH_FLOAT16 = 2,
  SMESH_FLOAT32 = 4,
  SMESH_FLOAT64 = 8,
  SMESH_INT8 = 10,
  SMESH_INT16 = 20,
  SMESH_INT32 = 40,
  SMESH_INT64 = 80,
  SMESH_UINT8 = 110,
  SMESH_UINT16 = 120,
  SMESH_UINT32 = 140,
  SMESH_UINT64 = 180,
  SMESH_TYPE_UNDEFINED = -1
};

template< typename T>
struct TypeToEnum {
  static enum PrimitiveType value() { return SMESH_TYPE_UNDEFINED; }
};

template <> struct TypeToEnum<f16> {
  static enum PrimitiveType value() { return SMESH_FLOAT16; }
};

template <> struct TypeToEnum<f32> {
  static enum PrimitiveType value() { return SMESH_FLOAT32; }
};

template <> struct TypeToEnum<f64> {
  static enum PrimitiveType value() { return SMESH_FLOAT64; }
};

template <> struct TypeToEnum<i16> {
  static enum PrimitiveType value() { return SMESH_INT16; }
};

template <> struct TypeToEnum<i32> {
  static enum PrimitiveType value() { return SMESH_INT32; }
};

template <> struct TypeToEnum<i64> {
  static enum PrimitiveType value() { return SMESH_INT64; }
};

template <> struct TypeToEnum<i8> {
  static enum PrimitiveType value() { return SMESH_INT8; }
};

template <> struct TypeToEnum<u8> {
  static enum PrimitiveType value() { return SMESH_UINT8; }
};

template <> struct TypeToEnum<u16> {
  static enum PrimitiveType value() { return SMESH_UINT16; }
};

template <> struct TypeToEnum<u32> {
  static enum PrimitiveType value() { return SMESH_UINT32; }
};

size_t num_bytes(enum PrimitiveType type) {
  switch (type) {
  case SMESH_FLOAT16:
    return sizeof(f16);
  case SMESH_FLOAT32:
    return sizeof(f32);
  case SMESH_FLOAT64:
    return sizeof(f64);
  case SMESH_INT16:
    return sizeof(i16);
  case SMESH_INT32:
    return sizeof(i32);
  case SMESH_INT64:
    return sizeof(i64);
  case SMESH_INT8:
    return sizeof(i8);
  case SMESH_UINT8:
    return sizeof(u8);
  case SMESH_UINT16:
    return sizeof(u16);
  case SMESH_UINT32:
    return sizeof(u32);
  case SMESH_UINT64:
    return sizeof(u64);
  default:
    SMESH_ERROR("Invalid primitive type: %d", type);
    return 0;
  }
}

inline std::string_view to_string(enum PrimitiveType type) {
  switch (type) {
  case SMESH_FLOAT16:
    return "float16";
  case SMESH_FLOAT32:
    return "float32";
  case SMESH_FLOAT64:
    return "float64";
  case SMESH_INT16:
    return "int16";
  case SMESH_INT32:
    return "int32";
  case SMESH_INT64:
    return "int64";
  case SMESH_UINT8:
    return "uint8";
  case SMESH_UINT16:
    return "uint16";
  case SMESH_UINT32:
    return "uint32";
  case SMESH_UINT64:
    return "uint64";
  default:
    SMESH_ERROR("Invalid primitive type: %d", type);
    return "undefined";
  }
}

inline PrimitiveType to_real_type(std::string_view type) {
  if (type == "float16") {
    return SMESH_FLOAT16;
  } else if (type == "float32") {
    return SMESH_FLOAT32;
  } else if (type == "float64") {
    return SMESH_FLOAT64;
  }
  return SMESH_TYPE_UNDEFINED;
}

inline PrimitiveType to_integer_type(std::string_view type) {
  if (type == "int16") {
    return SMESH_INT16;
  } else if (type == "int32") {
    return SMESH_INT32;
  } else if (type == "int64") {
    return SMESH_INT64;
  }
  return SMESH_TYPE_UNDEFINED;
}

template <typename T> T invalid_idx() { return static_cast<T>(-1); }

} // namespace smesh

#endif // SMESH_TYPES_HPP