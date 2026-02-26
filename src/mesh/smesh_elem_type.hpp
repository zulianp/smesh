#ifndef SMESH_ELEM_TYPE_HPP
#define SMESH_ELEM_TYPE_HPP

#include "smesh_base.hpp"

#include <assert.h>
#include <string.h>

namespace smesh {

enum ElemType {
  NIL = 0,
  NODE1 = 1,
  EDGE2 = 2,
  EDGE3 = 11,
  EDGESHELL2 = 101,
  BEAM2 = 100002,
  TRI3 = 3,
  TRI6 = 6,
  TRI10 = 1010,
  TRISHELL3 = 103,
  TRISHELL6 = 106,
  TRISHELL10 = 110,
  QUAD4 = 40,
  QUADSHELL4 = 140,
  TET4 = 4,
  TET10 = 10,
  TET15 = 15,
  TET20 = 20,
  HEX8 = 8,
  WEDGE6 = 1006,
  MACRO = 200,
  MACRO_TRI3 = (MACRO + TRI3),
  MACRO_TRISHELL3 = (MACRO + TRISHELL3),
  MACRO_TET4 = (MACRO + TET4),
  SSTET4 = 4000,
  SSQUAD4 = 40000,
  SSQUADSHELL4 = 140000,
  SSHEX8 = 8000,
  PROTEUS_HEX8 = 80000,
  PROTEUS_HEX27 = 270000,
  PROTEUS_HEX64 = 640000,
  PROTEUS_HEX125 = 1250000,
  PROTEUS_HEX216 = 2160000,
  PROTEUS_HEX343 = 3430000,
  PROTEUS_HEX512 = 5120000,
  PROTEUS_HEX729 = 7290000,
  INVALID = -1
};

inline enum ElemType type_from_string(const char *str) {
  if (!strcmp(str, "NODE1"))
    return NODE1;
  if (!strcmp(str, "EDGE2"))
    return EDGE2;
  if (!strcmp(str, "EDGE3"))
    return EDGE3;
  if (!strcmp(str, "TRI3"))
    return TRI3;
  if (!strcmp(str, "TRI6"))
    return TRI6;
  if (!strcmp(str, "TRI10"))
    return TRI10;
  if (!strcmp(str, "TRISHELL3"))
    return TRISHELL3;
  if (!strcmp(str, "WEDGE6"))
    return WEDGE6;
  if (!strcmp(str, "QUAD4"))
    return QUAD4;
  if (!strcmp(str, "QUADSHELL4"))
    return QUADSHELL4;
  if (!strcmp(str, "SSQUAD4"))
    return SSQUAD4;
  if (!strcmp(str, "SSQUADSHELL4"))
    return SSQUADSHELL4;
  if (!strcmp(str, "TET4"))
    return TET4;
  if (!strcmp(str, "TET15"))
    return TET15;
  if (!strcmp(str, "TET10"))
    return TET10;
  if (!strcmp(str, "TET20"))
    return TET20;
  if (!strcmp(str, "MACRO_TRI3"))
    return MACRO_TRI3;
  if (!strcmp(str, "MACRO_TET4"))
    return MACRO_TET4;
  if (!strcmp(str, "HEX8"))
    return HEX8;
  if (!strcmp(str, "SSHEX8"))
    return SSHEX8;
  if (!strcmp(str, "PROTEUS_HEX8"))
    return PROTEUS_HEX8;
  if (!strcmp(str, "PROTEUS_HEX27"))
    return PROTEUS_HEX27;
  if (!strcmp(str, "PROTEUS_HEX64"))
    return PROTEUS_HEX64;
  if (!strcmp(str, "PROTEUS_HEX125"))
    return PROTEUS_HEX125;
  if (!strcmp(str, "PROTEUS_HEX216"))
    return PROTEUS_HEX216;
  if (!strcmp(str, "PROTEUS_HEX343"))
    return PROTEUS_HEX343;
  if (!strcmp(str, "PROTEUS_HEX512"))
    return PROTEUS_HEX512;
  if (!strcmp(str, "PROTEUS_HEX729"))
    return PROTEUS_HEX729;

  assert(0);
  return INVALID;
}

inline const char *type_to_string(enum ElemType type) {
  switch (type) {
  case NODE1:
    return "NODE1";
  case EDGE2:
    return "EDGE2";
  case EDGE3:
    return "EDGE3";
  case BEAM2:
    return "BEAM2";
  case TRI3:
    return "TRI3";
  case TRISHELL3:
    return "TRISHELL3";
  case WEDGE6:
    return "WEDGE6";
  case QUAD4:
    return "QUAD4";
  case QUADSHELL4:
    return "QUADSHELL4";
  case SSQUAD4:
    return "SSQUAD4";
  case SSQUADSHELL4:
    return "SSQUADSHELL4";
  case TET4:
    return "TET4";
  case TRI6:
    return "TRI6";
  case TRISHELL6:
    return "TRISHELL6";
  case TRI10:
    return "TRI10";
  case MACRO_TRI3:
    return "MACRO_TRI3";
  case MACRO_TRISHELL3:
    return "MACRO_TRISHELL3";
  case MACRO_TET4:
    return "MACRO_TET4";
  case HEX8:
    return "HEX8";
  case SSHEX8:
    return "SSHEX8";
  case TET10:
    return "TET10";
  case TET15:
    return "TET15";
  case TET20:
    return "TET20";
  case PROTEUS_HEX8:
    return "PROTEUS_HEX8";
  case PROTEUS_HEX27:
    return "PROTEUS_HEX27";
  case PROTEUS_HEX64:
    return "PROTEUS_HEX64";
  case PROTEUS_HEX125:
    return "PROTEUS_HEX125";
  case PROTEUS_HEX216:
    return "PROTEUS_HEX216";
  case PROTEUS_HEX343:
    return "PROTEUS_HEX343";
  case PROTEUS_HEX512:
    return "PROTEUS_HEX512";
  case PROTEUS_HEX729:
    return "PROTEUS_HEX729";
  default: {
    assert(0);
    return "INVALID";
  }
  }
}

inline enum ElemType side_type(const enum ElemType type) {
  switch (type) {
  case TRI3:
  case QUAD4:
    return EDGE2;
  case TRI6:
    return EDGE3;
  case TET4:
    return TRI3;
  case TET10:
    return TRI6;
  case TET20:
    return TRI10;
  case EDGE2:
    return NODE1;
  case TRISHELL3:
    return BEAM2;
  case QUADSHELL4:
    return BEAM2;
  case MACRO_TET4:
    return MACRO_TRI3;
  case HEX8:
    return QUAD4;
  case SSHEX8:
    return SSQUAD4;
  case PROTEUS_HEX216:
    return PROTEUS_HEX8;
  case PROTEUS_HEX343:
    return PROTEUS_HEX8;
  case PROTEUS_HEX512:
    return PROTEUS_HEX8;
  case PROTEUS_HEX729:
    return PROTEUS_HEX8;
  default: {
    assert(0);
    return INVALID;
  }
  }
}

inline enum ElemType shell_type(const enum ElemType type) {
  switch (type) {
  case TRI3:
    return TRISHELL3;
  case MACRO_TRI3:
    return MACRO_TRISHELL3;
  case TRI6:
    return TRISHELL6;
  case TRI10:
    return TRISHELL10;
  case TRISHELL3:
    return TRISHELL3;
  case TRISHELL6:
    return TRISHELL6;
  case TRISHELL10:
    return TRISHELL10;
  case EDGE2:
    return EDGESHELL2;
  case BEAM2:
    return BEAM2;
  case QUAD4:
    return QUADSHELL4;
  case QUADSHELL4:
    return QUADSHELL4;
  case SSQUAD4:
    return SSQUADSHELL4;
  case SSQUADSHELL4:
    return SSQUADSHELL4;
  case PROTEUS_HEX216:
    return PROTEUS_HEX8;
  case PROTEUS_HEX343:
    return PROTEUS_HEX8;
  case PROTEUS_HEX512:
    return PROTEUS_HEX8;
  case PROTEUS_HEX729:
    return PROTEUS_HEX8;
  default: {
    // assert(0);
    return INVALID;
  }
  }
}

inline enum ElemType elem_lower_order(const enum ElemType type) {
  switch (type) {
  case NIL:
    return NIL;
  case TRI6:
    return TRI3;
  case TET10:
    return TET4;
  case TET20:
    return TET10;
  case EDGE3:
    return EDGE2;
  default: {
    assert(0);
    return INVALID;
  }
  }
}

inline enum ElemType elem_higher_order(const enum ElemType type) {
  switch (type) {
  case NIL:
    return NIL;
  case TRI3:
    return TRI6;
  case TET4:
    return TET10;
  case TET10:
    return TET20;
  case EDGE2:
    return EDGE3;
  default: {
    assert(0);
    return INVALID;
  }
  }
}

inline int elem_num_nodes(const enum ElemType type) {
  switch (type) {
  case NIL:
    return 0;
  case NODE1:
    return 1;
  case EDGE2:
    return 2;
  case EDGE3:
    return 3;
  case TRI3:
    return 3;
  case TRISHELL3:
    return 3;
  case WEDGE6:
    return 6;
  case QUAD4:
    return 4;
  case QUADSHELL4:
    return 4;
  case TET4:
    return 4;
  case TRI6:
    return 6;
  case MACRO_TRI3:
    return 6;
  case MACRO_TET4:
    return 10;
  case HEX8:
    return 8;
  case TET10:
    return 10;
  case TET15:
    return 15;
  case TET20:
    return 20;
  case PROTEUS_HEX8:
    return 8;
  case PROTEUS_HEX27:
    return 27;
  case PROTEUS_HEX64:
    return 64;
  case PROTEUS_HEX125:
    return 125;
  case PROTEUS_HEX216:
    return 216;
  case PROTEUS_HEX343:
    return 343;
  case PROTEUS_HEX512:
    return 512;
  case PROTEUS_HEX729:
    return 729;
  default: {
    assert(0);
    return 0;
  }
  }
}

inline int elem_num_sides(const enum ElemType type) {
  switch (type) {
  case NIL:
    return 0;
  case EDGE2:
    return 2;
  case TRI3:
    return 3;
  case TRISHELL3:
    return 3;
  case MACRO_TRI3:
    return 3; // Really?
  case MACRO_TET4:
    return 4;
  case QUAD4:
    return 4;
  case QUADSHELL4:
    return 4;
  case TET4:
    return 4;
  case WEDGE6:
    return 5;
  case TRI6:
    return 3;
  case HEX8:
    return 6;
  case TET10:
    return 4;
  case TET15:
    return 4;
  case TET20:
    return 4;
  case PROTEUS_HEX8:
    return 8;
  case PROTEUS_HEX27:
    return 8;
  case PROTEUS_HEX64:
    return 8;
  case PROTEUS_HEX125:
    return 8;
  default: {
    assert(0);
    return 0;
  }
  }
}

inline int elem_manifold_dim(const enum ElemType type) {
  switch (type) {
  case NIL:
    return 0;
  case EDGE2:
    return 1;
  case TRI3:
    return 2;
  case QUAD4:
    return 2;
  case QUADSHELL4:
    return 2;
  case TET4:
    return 3;
  case WEDGE6:
    return 3;
  case TRI6:
    return 2;
  case MACRO_TRI3:
    return 2;
  case MACRO_TET4:
    return 3;
  case HEX8:
    return 3;
  case TET10:
    return 3;
  case TET15:
    return 3;
  case TET20:
    return 3;
  case PROTEUS_HEX8:
    return 3;
  case PROTEUS_HEX27:
    return 3;
  case PROTEUS_HEX64:
    return 3;
  case PROTEUS_HEX125:
    return 3;
  default: {
    assert(0);
    return INVALID;
  }
  }
}

inline enum ElemType macro_type_variant(const enum ElemType type) {
  switch (type) {
  case TET10:
    return MACRO_TET4;
  case TRI6:
    return MACRO_TRI3;

  default: {
    assert(0);
    return type;
  }
  }
}

inline enum ElemType macro_base_elem(const enum ElemType macro_type) {
  switch (macro_type) {
  case MACRO_TET4:
    return TET4;
  case MACRO_TRI3:
    return TRI3;
  case TET10:
    return TET4;
  case TRI6:
    return TRI3;
  case SSHEX8:
    return HEX8;
  case PROTEUS_HEX8:
    return PROTEUS_HEX8;
  case PROTEUS_HEX27:
    return PROTEUS_HEX8;
  case PROTEUS_HEX64:
    return PROTEUS_HEX8;
  case PROTEUS_HEX125:
    return PROTEUS_HEX8;
  default: {
    assert(0);
    return macro_type;
  }
  }
}

inline int is_second_order_lagrange(const enum ElemType type) {
  switch (type) {
  case TET10:
    return 1;
  case TRI6:
    return 1;
  default: {
    return 0;
  }
  }
}

inline enum ElemType proteus_hex_type(const int micro_elements_per_dim) {
  switch (micro_elements_per_dim) {
  case 1:
    return PROTEUS_HEX8;
  case 2:
    return PROTEUS_HEX27;
  case 3:
    return PROTEUS_HEX64;
  case 4:
    return PROTEUS_HEX125;
    case 5:
    return PROTEUS_HEX216;
  case 6:
    return PROTEUS_HEX343;
  case 7:
    return PROTEUS_HEX512;
  case 8:
    return PROTEUS_HEX729;
  default:
    assert(0);
    return INVALID;
  }
}

enum HEX8_Sides {
  HEX8_LEFT = 3,
  HEX8_RIGHT = 1,
  HEX8_BOTTOM = 4,
  HEX8_TOP = 5,
  HEX8_FRONT = 0,
  HEX8_BACK = 2
};

} // namespace smesh

#endif // SMESH_ELEM_TYPE_HPP