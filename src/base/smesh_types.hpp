#ifndef SMESH_TYPES_HPP
#define SMESH_TYPES_HPP

#include "smesh_base.hpp"
#include <string_view>

namespace smesh {

using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using mask_t = char;


#if defined(__APPLE__)
using f16 = __fp16;
#else
using f16 = _Float16;
#endif

using f32 = float;
using f64 = double;
using real_t = f64;

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

enum RealType {
  SMESH_FLOAT16 = 2,
  SMESH_FLOAT32 = 4,
  SMESH_FLOAT64 = 8,
  SMESH_REAL_DEFAULT = 0,
  SMESH_REAL_UNDEFINED = -1
};

enum IntegerType {
  SMESH_INT16 = 2,
  SMESH_INT32 = 4,
  SMESH_INT64 = 8,
  SMESH_INT_DEFAULT = 0,
  SMESH_INT_UNDEFINED = -1
};

inline std::string_view to_string(enum RealType type) {
  switch (type) {
  case SMESH_FLOAT16:
    return "float16";
  case SMESH_FLOAT32:
    return "float32";
  case SMESH_FLOAT64:
    return "float64";
  default:
    return "undefined";
  }
}

inline std::string_view to_string(enum IntegerType type) {
  switch (type) {
  case SMESH_INT16:
    return "int16";
  case SMESH_INT32:
    return "int32";
  case SMESH_INT64:
    return "int64";
  default:
    return "undefined";
  }
}

inline RealType to_real_type(std::string_view type) {
  if (type == "float16") {
    return SMESH_FLOAT16;
  } else if (type == "float32") {
    return SMESH_FLOAT32;
  } else if (type == "float64") {
    return SMESH_FLOAT64;
  }
  return SMESH_REAL_UNDEFINED;
}

inline IntegerType to_integer_type(std::string_view type) {
  if (type == "int16") {
    return SMESH_INT16;
  } else if (type == "int32") {
    return SMESH_INT32;
  } else if (type == "int64") {
    return SMESH_INT64;
  }
  return SMESH_INT_UNDEFINED;
}

template <typename T>
T invalid_idx() {
  return static_cast<T>(-1);
}

} // namespace smesh

#endif // SMESH_TYPES_HPP