#ifndef SMESH_ELEM_TYPE_HPP
#define SMESH_ELEM_TYPE_HPP

#include "smesh_base.hpp"

#include <assert.h>
#include <string.h>

namespace smesh {

    enum ElemType {
        NIL             = 0,
        NODE1           = 1,
        EDGE2           = 2,
        EDGE3           = 11,
        EDGESHELL2      = 101,
        EDGESHELL3      = 111,
        BEAM2           = 100002,
        TRI3            = 3,
        TRI6            = 6,
        TRI10           = 1010,
        TRISHELL3       = 103,
        TRISHELL6       = 106,
        TRISHELL10      = 110,
        QUAD4           = 40,
        QUAD9           = 9,
        QUADSHELL4      = 140,
        QUADSHELL9      = 149,
        TET4            = 4,
        TET10           = 10,
        TET15           = 15,
        TET20           = 20,
        HEX8            = 8,
        HEX27           = 27,
        PYRAMID5        = 1005,
        WEDGE6          = 1006,
        MACRO           = 200,
        MACRO_TRI3      = (MACRO + TRI3),
        MACRO_TRISHELL3 = (MACRO + TRISHELL3),
        MACRO_TET4      = (MACRO + TET4),
        // Semistructured types only below this line
        SEMISTRUCTURED_ELEMENT = 100000,
        PROTEUS_HEX8           = 100008,
        PROTEUS_HEX27          = 270000,
        PROTEUS_HEX64          = 640000,
        PROTEUS_HEX125         = 1250000,
        PROTEUS_HEX216         = 2160000,
        PROTEUS_HEX343         = 3430000,
        PROTEUS_HEX512         = 5120000,
        PROTEUS_HEX729         = 7290000,
        PROTEUS_HEX4913        = 49130000,
        PROTEUS_QUAD4          = 400000,
        PROTEUS_QUAD9          = 900000,
        PROTEUS_QUAD16         = 1600000,
        PROTEUS_QUAD25         = 2500000,
        PROTEUS_QUAD36         = 3600000,
        PROTEUS_QUAD49         = 4900000,
        PROTEUS_QUAD64         = 6400000,
        PROTEUS_QUAD81         = 8100000,
        PROTEUS_QUAD289        = 28900000,
        PROTEUS_QUADSHELL4     = 400001,
        PROTEUS_QUADSHELL9     = 900001,
        PROTEUS_QUADSHELL16    = 1600001,
        PROTEUS_QUADSHELL25    = 2500001,
        PROTEUS_QUADSHELL36    = 3600001,
        PROTEUS_QUADSHELL49    = 4900001,
        PROTEUS_QUADSHELL64    = 6400001,
        PROTEUS_QUADSHELL81    = 8100001,
        PROTEUS_QUADSHELL289   = 28900001,
        // PROTEUS_TET{N}: N = (L+1)(L+2)(L+3)/6 lattice nodes. L=1/2 cannot use N*10000
        // (40000 < SEMISTRUCTURED_ELEMENT; 100000 collides with the sentinel).
        PROTEUS_TET4           = 100004,
        PROTEUS_TET10          = 100010,
        PROTEUS_TET20          = 200000,
        PROTEUS_TET35          = 350000,
        PROTEUS_TET56          = 560000,
        PROTEUS_TET84          = 840000,
        PROTEUS_TET120         = 1200000,
        PROTEUS_TET165         = 1650000,
        PROTEUS_TET969         = 9690000,
        // PROTEUS_WEDGE{N} / PROTEUS_PYRAMID{N}: N = lattice nodes. L=1 uses 100006 / 100005
        // (next to TET4/HEX8). L>=2 uses dedicated bands (N*10000 collides, e.g. WEDGE40=QUAD4).
        PROTEUS_WEDGE6         = 100006,
        PROTEUS_WEDGE18        = 8000018,
        PROTEUS_WEDGE40        = 8000040,
        PROTEUS_WEDGE75        = 8000075,
        PROTEUS_WEDGE126       = 8000126,
        PROTEUS_WEDGE196       = 8000196,
        PROTEUS_WEDGE288       = 8000288,
        PROTEUS_WEDGE405       = 8000405,
        PROTEUS_WEDGE2601      = 8002601,
        PROTEUS_PYRAMID5       = 100005,
        PROTEUS_PYRAMID14      = 9000014,
        PROTEUS_PYRAMID30      = 9000030,
        PROTEUS_PYRAMID55      = 9000055,
        PROTEUS_PYRAMID91      = 9000091,
        PROTEUS_PYRAMID140     = 9000140,
        PROTEUS_PYRAMID204     = 9000204,
        PROTEUS_PYRAMID285     = 9000285,
        PROTEUS_PYRAMID1785    = 9001785,
        INVALID                = -1
    };

    /// True if \c type is in the semistructured range (\c type >= SEMISTRUCTURED_ELEMENT), including PROTEUS_HEX8.
    inline bool is_semistructured_type(const enum ElemType type) { return type >= SEMISTRUCTURED_ELEMENT; }

    inline enum ElemType type_from_string(const char *str) {
        if (!strcmp(str, "NODE1")) return NODE1;
        if (!strcmp(str, "EDGE2")) return EDGE2;
        if (!strcmp(str, "EDGE3")) return EDGE3;
        if (!strcmp(str, "EDGESHELL2")) return EDGESHELL2;
        if (!strcmp(str, "EDGESHELL3")) return EDGESHELL3;
        if (!strcmp(str, "TRI3")) return TRI3;
        if (!strcmp(str, "TRI6")) return TRI6;
        if (!strcmp(str, "TRI10")) return TRI10;
        if (!strcmp(str, "TRISHELL3")) return TRISHELL3;
        if (!strcmp(str, "PYRAMID5")) return PYRAMID5;
        if (!strcmp(str, "WEDGE6")) return WEDGE6;
        if (!strcmp(str, "QUAD4")) return QUAD4;
        if (!strcmp(str, "QUAD9")) return QUAD9;
        if (!strcmp(str, "QUADSHELL4")) return QUADSHELL4;
        if (!strcmp(str, "QUADSHELL9")) return QUADSHELL9;
        if (!strcmp(str, "TET4")) return TET4;
        if (!strcmp(str, "TET15")) return TET15;
        if (!strcmp(str, "TET10")) return TET10;
        if (!strcmp(str, "TET20")) return TET20;
        if (!strcmp(str, "MACRO_TRI3")) return MACRO_TRI3;
        if (!strcmp(str, "MACRO_TET4")) return MACRO_TET4;
        if (!strcmp(str, "HEX8")) return HEX8;
        if (!strcmp(str, "HEX27")) return HEX27;
        if (!strcmp(str, "PROTEUS_HEX8")) return PROTEUS_HEX8;
        if (!strcmp(str, "PROTEUS_HEX27")) return PROTEUS_HEX27;
        if (!strcmp(str, "PROTEUS_HEX64")) return PROTEUS_HEX64;
        if (!strcmp(str, "PROTEUS_HEX125")) return PROTEUS_HEX125;
        if (!strcmp(str, "PROTEUS_HEX216")) return PROTEUS_HEX216;
        if (!strcmp(str, "PROTEUS_HEX343")) return PROTEUS_HEX343;
        if (!strcmp(str, "PROTEUS_HEX512")) return PROTEUS_HEX512;
        if (!strcmp(str, "PROTEUS_HEX729")) return PROTEUS_HEX729;
        if (!strcmp(str, "PROTEUS_HEX4913")) return PROTEUS_HEX4913;
        if (!strcmp(str, "PROTEUS_QUAD4")) return PROTEUS_QUAD4;
        if (!strcmp(str, "PROTEUS_QUAD9")) return PROTEUS_QUAD9;
        if (!strcmp(str, "PROTEUS_QUAD16")) return PROTEUS_QUAD16;
        if (!strcmp(str, "PROTEUS_QUAD25")) return PROTEUS_QUAD25;
        if (!strcmp(str, "PROTEUS_QUAD36")) return PROTEUS_QUAD36;
        if (!strcmp(str, "PROTEUS_QUAD49")) return PROTEUS_QUAD49;
        if (!strcmp(str, "PROTEUS_QUAD64")) return PROTEUS_QUAD64;
        if (!strcmp(str, "PROTEUS_QUAD81")) return PROTEUS_QUAD81;
        if (!strcmp(str, "PROTEUS_QUAD289")) return PROTEUS_QUAD289;
        if (!strcmp(str, "PROTEUS_QUADSHELL4")) return PROTEUS_QUADSHELL4;
        if (!strcmp(str, "PROTEUS_QUADSHELL9")) return PROTEUS_QUADSHELL9;
        if (!strcmp(str, "PROTEUS_QUADSHELL16")) return PROTEUS_QUADSHELL16;
        if (!strcmp(str, "PROTEUS_QUADSHELL25")) return PROTEUS_QUADSHELL25;
        if (!strcmp(str, "PROTEUS_QUADSHELL36")) return PROTEUS_QUADSHELL36;
        if (!strcmp(str, "PROTEUS_QUADSHELL49")) return PROTEUS_QUADSHELL49;
        if (!strcmp(str, "PROTEUS_QUADSHELL64")) return PROTEUS_QUADSHELL64;
        if (!strcmp(str, "PROTEUS_QUADSHELL81")) return PROTEUS_QUADSHELL81;
        if (!strcmp(str, "PROTEUS_QUADSHELL289")) return PROTEUS_QUADSHELL289;
        if (!strcmp(str, "PROTEUS_TET4")) return PROTEUS_TET4;
        if (!strcmp(str, "PROTEUS_TET10")) return PROTEUS_TET10;
        if (!strcmp(str, "PROTEUS_TET20")) return PROTEUS_TET20;
        if (!strcmp(str, "PROTEUS_TET35")) return PROTEUS_TET35;
        if (!strcmp(str, "PROTEUS_TET56")) return PROTEUS_TET56;
        if (!strcmp(str, "PROTEUS_TET84")) return PROTEUS_TET84;
        if (!strcmp(str, "PROTEUS_TET120")) return PROTEUS_TET120;
        if (!strcmp(str, "PROTEUS_TET165")) return PROTEUS_TET165;
        if (!strcmp(str, "PROTEUS_TET969")) return PROTEUS_TET969;
        if (!strcmp(str, "PROTEUS_WEDGE6")) return PROTEUS_WEDGE6;
        if (!strcmp(str, "PROTEUS_WEDGE18")) return PROTEUS_WEDGE18;
        if (!strcmp(str, "PROTEUS_WEDGE40")) return PROTEUS_WEDGE40;
        if (!strcmp(str, "PROTEUS_WEDGE75")) return PROTEUS_WEDGE75;
        if (!strcmp(str, "PROTEUS_WEDGE126")) return PROTEUS_WEDGE126;
        if (!strcmp(str, "PROTEUS_WEDGE196")) return PROTEUS_WEDGE196;
        if (!strcmp(str, "PROTEUS_WEDGE288")) return PROTEUS_WEDGE288;
        if (!strcmp(str, "PROTEUS_WEDGE405")) return PROTEUS_WEDGE405;
        if (!strcmp(str, "PROTEUS_WEDGE2601")) return PROTEUS_WEDGE2601;
        if (!strcmp(str, "PROTEUS_PYRAMID5")) return PROTEUS_PYRAMID5;
        if (!strcmp(str, "PROTEUS_PYRAMID14")) return PROTEUS_PYRAMID14;
        if (!strcmp(str, "PROTEUS_PYRAMID30")) return PROTEUS_PYRAMID30;
        if (!strcmp(str, "PROTEUS_PYRAMID55")) return PROTEUS_PYRAMID55;
        if (!strcmp(str, "PROTEUS_PYRAMID91")) return PROTEUS_PYRAMID91;
        if (!strcmp(str, "PROTEUS_PYRAMID140")) return PROTEUS_PYRAMID140;
        if (!strcmp(str, "PROTEUS_PYRAMID204")) return PROTEUS_PYRAMID204;
        if (!strcmp(str, "PROTEUS_PYRAMID285")) return PROTEUS_PYRAMID285;
        if (!strcmp(str, "PROTEUS_PYRAMID1785")) return PROTEUS_PYRAMID1785;
        SMESH_ERROR("No element type found for string: %s\n", str);
        return INVALID;
    }

    inline const char *type_to_string(enum ElemType type) {
        switch (type) {
            case NODE1:
                return "NODE1";
            case EDGE2:
                return "EDGE2";
            case EDGESHELL2:
                return "EDGESHELL2";
            case EDGESHELL3:
                return "EDGESHELL3";
            case EDGE3:
                return "EDGE3";
            case BEAM2:
                return "BEAM2";
            case TRI3:
                return "TRI3";
            case TRISHELL3:
                return "TRISHELL3";
            case PYRAMID5:
                return "PYRAMID5";
            case WEDGE6:
                return "WEDGE6";
            case QUAD4:
                return "QUAD4";
            case QUAD9:
                return "QUAD9";
            case QUADSHELL4:
                return "QUADSHELL4";
            case QUADSHELL9:
                return "QUADSHELL9";
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
            case HEX27:
                return "HEX27";
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
            case PROTEUS_HEX4913:
                return "PROTEUS_HEX4913";
            case PROTEUS_QUAD4:
                return "PROTEUS_QUAD4";
            case PROTEUS_QUAD9:
                return "PROTEUS_QUAD9";
            case PROTEUS_QUAD16:
                return "PROTEUS_QUAD16";
            case PROTEUS_QUAD25:
                return "PROTEUS_QUAD25";
            case PROTEUS_QUAD36:
                return "PROTEUS_QUAD36";
            case PROTEUS_QUAD49:
                return "PROTEUS_QUAD49";
            case PROTEUS_QUAD64:
                return "PROTEUS_QUAD64";
            case PROTEUS_QUAD81:
                return "PROTEUS_QUAD81";
            case PROTEUS_QUAD289:
                return "PROTEUS_QUAD289";
            case PROTEUS_QUADSHELL4:
                return "PROTEUS_QUADSHELL4";
            case PROTEUS_QUADSHELL9:
                return "PROTEUS_QUADSHELL9";
            case PROTEUS_QUADSHELL16:
                return "PROTEUS_QUADSHELL16";
            case PROTEUS_QUADSHELL25:
                return "PROTEUS_QUADSHELL25";
            case PROTEUS_QUADSHELL36:
                return "PROTEUS_QUADSHELL36";
            case PROTEUS_QUADSHELL49:
                return "PROTEUS_QUADSHELL49";
            case PROTEUS_QUADSHELL64:
                return "PROTEUS_QUADSHELL64";
            case PROTEUS_QUADSHELL81:
                return "PROTEUS_QUADSHELL81";
            case PROTEUS_QUADSHELL289:
                return "PROTEUS_QUADSHELL289";
            case PROTEUS_TET4:
                return "PROTEUS_TET4";
            case PROTEUS_TET10:
                return "PROTEUS_TET10";
            case PROTEUS_TET20:
                return "PROTEUS_TET20";
            case PROTEUS_TET35:
                return "PROTEUS_TET35";
            case PROTEUS_TET56:
                return "PROTEUS_TET56";
            case PROTEUS_TET84:
                return "PROTEUS_TET84";
            case PROTEUS_TET120:
                return "PROTEUS_TET120";
            case PROTEUS_TET165:
                return "PROTEUS_TET165";
            case PROTEUS_TET969:
                return "PROTEUS_TET969";
            case PROTEUS_WEDGE6:
                return "PROTEUS_WEDGE6";
            case PROTEUS_WEDGE18:
                return "PROTEUS_WEDGE18";
            case PROTEUS_WEDGE40:
                return "PROTEUS_WEDGE40";
            case PROTEUS_WEDGE75:
                return "PROTEUS_WEDGE75";
            case PROTEUS_WEDGE126:
                return "PROTEUS_WEDGE126";
            case PROTEUS_WEDGE196:
                return "PROTEUS_WEDGE196";
            case PROTEUS_WEDGE288:
                return "PROTEUS_WEDGE288";
            case PROTEUS_WEDGE405:
                return "PROTEUS_WEDGE405";
            case PROTEUS_WEDGE2601:
                return "PROTEUS_WEDGE2601";
            case PROTEUS_PYRAMID5:
                return "PROTEUS_PYRAMID5";
            case PROTEUS_PYRAMID14:
                return "PROTEUS_PYRAMID14";
            case PROTEUS_PYRAMID30:
                return "PROTEUS_PYRAMID30";
            case PROTEUS_PYRAMID55:
                return "PROTEUS_PYRAMID55";
            case PROTEUS_PYRAMID91:
                return "PROTEUS_PYRAMID91";
            case PROTEUS_PYRAMID140:
                return "PROTEUS_PYRAMID140";
            case PROTEUS_PYRAMID204:
                return "PROTEUS_PYRAMID204";
            case PROTEUS_PYRAMID285:
                return "PROTEUS_PYRAMID285";
            case PROTEUS_PYRAMID1785:
                return "PROTEUS_PYRAMID1785";
            default: {
                SMESH_ERROR("No element type found for type: %d\n", type);
                return "INVALID";
            }
        }
    }

    inline enum ElemType side_type(const enum ElemType type) {
        switch (type) {
            case TRI3:
            case QUAD4:
            case PROTEUS_QUAD4:
            case PROTEUS_QUAD9:
            case PROTEUS_QUAD16:
            case PROTEUS_QUAD25:
            case PROTEUS_QUAD36:
            case PROTEUS_QUAD49:
            case PROTEUS_QUAD64:
            case PROTEUS_QUAD81:
            case PROTEUS_QUAD289:
                return EDGE2;
            case QUAD9:
                return EDGE3;
            case QUADSHELL9:
                return EDGE3;
            case TRI6:
                return EDGE3;
            case TET4:
                return TRI3;
            case TET10:
                return TRI6;
            case TET20:
                return TRI10;
            case EDGE2:
            case EDGE3:
            case EDGESHELL2:
            case EDGESHELL3:
                return NODE1;
            case TRISHELL3:
                return BEAM2;
            case QUADSHELL4:
                return BEAM2;
            case MACRO_TET4:
                return MACRO_TRI3;
            case HEX8:
                return QUAD4;
            case HEX27:
                return QUAD9;
            case PROTEUS_HEX8:
                return PROTEUS_QUAD4;
            case PROTEUS_HEX27:
                return PROTEUS_QUAD9;
            case PROTEUS_HEX64:
                return PROTEUS_QUAD16;
            case PROTEUS_HEX125:
                return PROTEUS_QUAD25;
            case PROTEUS_HEX216:
                return PROTEUS_QUAD36;
            case PROTEUS_HEX343:
                return PROTEUS_QUAD49;
            case PROTEUS_HEX512:
                return PROTEUS_QUAD64;
            case PROTEUS_HEX729:
                return PROTEUS_QUAD81;
            case PROTEUS_HEX4913:
                return PROTEUS_QUAD289;
            case PROTEUS_TET4:
            case PROTEUS_TET10:
            case PROTEUS_TET20:
            case PROTEUS_TET35:
            case PROTEUS_TET56:
            case PROTEUS_TET84:
            case PROTEUS_TET120:
            case PROTEUS_TET165:
            case PROTEUS_TET969:
                return TRI3;
            default: {
                SMESH_ERROR("No side type found for type: %s\n", type_to_string(type));
                return INVALID;
            }
        }
    }

    inline bool elem_sides_homogeneous(const enum ElemType type) {
        switch (type) {
            case WEDGE6:
            case PYRAMID5:
            case PROTEUS_WEDGE6:
            case PROTEUS_WEDGE18:
            case PROTEUS_WEDGE40:
            case PROTEUS_WEDGE75:
            case PROTEUS_WEDGE126:
            case PROTEUS_WEDGE196:
            case PROTEUS_WEDGE288:
            case PROTEUS_WEDGE405:
            case PROTEUS_WEDGE2601:
            case PROTEUS_PYRAMID5:
            case PROTEUS_PYRAMID14:
            case PROTEUS_PYRAMID30:
            case PROTEUS_PYRAMID55:
            case PROTEUS_PYRAMID91:
            case PROTEUS_PYRAMID140:
            case PROTEUS_PYRAMID204:
            case PROTEUS_PYRAMID285:
            case PROTEUS_PYRAMID1785:
                return false;
            default:
                return true;
        }
    }

    /// Per-side face type. WEDGE6/PYRAMID5 mix TRI3 and QUAD4; homogeneous types
    /// ignore \c s and match \c side_type(type).
    inline enum ElemType side_type(const enum ElemType type, const int s) {
        switch (type) {
            case WEDGE6:
            case PROTEUS_WEDGE6:
            case PROTEUS_WEDGE18:
            case PROTEUS_WEDGE40:
            case PROTEUS_WEDGE75:
            case PROTEUS_WEDGE126:
            case PROTEUS_WEDGE196:
            case PROTEUS_WEDGE288:
            case PROTEUS_WEDGE405:
            case PROTEUS_WEDGE2601:
                return (s < 3) ? QUAD4 : TRI3;
            case PYRAMID5:
            case PROTEUS_PYRAMID5:
            case PROTEUS_PYRAMID14:
            case PROTEUS_PYRAMID30:
            case PROTEUS_PYRAMID55:
            case PROTEUS_PYRAMID91:
            case PROTEUS_PYRAMID140:
            case PROTEUS_PYRAMID204:
            case PROTEUS_PYRAMID285:
            case PROTEUS_PYRAMID1785:
                return (s < 4) ? TRI3 : QUAD4;
            default:
                (void)s;
                return side_type(type);
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
            case EDGE3:
                return EDGESHELL3;
            case EDGESHELL2:
                return EDGESHELL2;
            case EDGESHELL3:
                return EDGESHELL3;
            case BEAM2:
                return BEAM2;
            case QUAD4:
                return QUADSHELL4;
            case QUADSHELL4:
                return QUADSHELL4;
            case QUAD9:
                return QUADSHELL9;
            case QUADSHELL9:
                return QUADSHELL9;
            case PROTEUS_HEX8:
                return PROTEUS_QUADSHELL4;
            case PROTEUS_HEX27:
                return PROTEUS_QUADSHELL9;
            case PROTEUS_HEX64:
                return PROTEUS_QUADSHELL16;
            case PROTEUS_HEX125:
                return PROTEUS_QUADSHELL25;
            case PROTEUS_HEX216:
                return PROTEUS_QUADSHELL36;
            case PROTEUS_HEX343:
                return PROTEUS_QUADSHELL49;
            case PROTEUS_HEX512:
                return PROTEUS_QUADSHELL64;
            case PROTEUS_HEX729:
                return PROTEUS_QUADSHELL81;
            case PROTEUS_HEX4913:
                return PROTEUS_QUADSHELL289;
            case PROTEUS_QUAD4:
                return PROTEUS_QUADSHELL4;
            case PROTEUS_QUAD9:
                return PROTEUS_QUADSHELL9;
            case PROTEUS_QUAD16:
                return PROTEUS_QUADSHELL16;
            case PROTEUS_QUAD25:
                return PROTEUS_QUADSHELL25;
            case PROTEUS_QUAD36:
                return PROTEUS_QUADSHELL36;
            case PROTEUS_QUAD49:
                return PROTEUS_QUADSHELL49;
            case PROTEUS_QUAD64:
                return PROTEUS_QUADSHELL64;
            case PROTEUS_QUAD81:
                return PROTEUS_QUADSHELL81;
            case PROTEUS_QUAD289:
                return PROTEUS_QUADSHELL289;
            case PROTEUS_QUADSHELL4:
                return PROTEUS_QUADSHELL4;
            case PROTEUS_QUADSHELL9:
                return PROTEUS_QUADSHELL9;
            case PROTEUS_QUADSHELL16:
                return PROTEUS_QUADSHELL16;
            case PROTEUS_QUADSHELL25:
                return PROTEUS_QUADSHELL25;
            case PROTEUS_QUADSHELL36:
                return PROTEUS_QUADSHELL36;
            case PROTEUS_QUADSHELL289:
                return PROTEUS_QUADSHELL289;
            case PROTEUS_QUADSHELL64:
                return PROTEUS_QUADSHELL64;
            case PROTEUS_QUADSHELL81:
                return PROTEUS_QUADSHELL81;
            default: {
                SMESH_ERROR("No shell type found for type: %s\n", type_to_string(type));
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
            case QUAD9:
                return QUAD4;
            case TET10:
                return TET4;
            case HEX27:
                return HEX8;
            case TET20:
                return TET10;
            case EDGE3:
                return EDGE2;
            case EDGESHELL3:
                return EDGESHELL2;
            default: {
                SMESH_ERROR("No lower order type found for type: %s\n", type_to_string(type));
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
            case QUAD4:
                return QUAD9;
            case TET4:
                return TET10;
            case HEX8:
                return HEX27;
            case TET10:
                return TET20;
            case EDGE2:
                return EDGE3;
            case EDGESHELL2:
                return EDGESHELL3;
            default: {
                SMESH_ERROR("No higher order type found for type: %s\n", type_to_string(type));
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
            case EDGESHELL2:
                return 2;
            case EDGE3:
            case EDGESHELL3:
                return 3;
            case TRI3:
                return 3;
            case TRISHELL3:
                return 3;
            case TRISHELL6:
                return 6;
            case PYRAMID5:
                return 5;
            case WEDGE6:
                return 6;
            case QUAD4:
                return 4;
            case QUAD9:
                return 9;
            case QUADSHELL4:
                return 4;
            case QUADSHELL9:
                return 9;
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
            case HEX27:
                return 27;
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
            case PROTEUS_QUAD4:
                return 4;
            case PROTEUS_QUAD9:
                return 9;
            case PROTEUS_QUAD16:
                return 16;
            case PROTEUS_QUAD25:
                return 25;
            case PROTEUS_QUAD36:
                return 36;
            case PROTEUS_QUAD49:
                return 49;
            case PROTEUS_QUAD64:
                return 64;
            case PROTEUS_QUAD81:
                return 81;
            case PROTEUS_QUADSHELL4:
                return 4;
            case PROTEUS_QUADSHELL9:
                return 9;
            case PROTEUS_QUADSHELL16:
                return 16;
            case PROTEUS_QUADSHELL25:
                return 25;
            case PROTEUS_QUADSHELL36:
                return 36;
            case PROTEUS_QUADSHELL49:
                return 49;
            case PROTEUS_QUADSHELL64:
                return 64;
            case PROTEUS_QUADSHELL81:
                return 81;
            case PROTEUS_HEX4913:
                return 4913;
            case PROTEUS_QUAD289:
                return 289;
            case PROTEUS_TET4:
                return 4;
            case PROTEUS_TET10:
                return 10;
            case PROTEUS_TET20:
                return 20;
            case PROTEUS_TET35:
                return 35;
            case PROTEUS_TET56:
                return 56;
            case PROTEUS_TET84:
                return 84;
            case PROTEUS_TET120:
                return 120;
            case PROTEUS_TET165:
                return 165;
            case PROTEUS_TET969:
                return 969;
            case PROTEUS_WEDGE6:
                return 6;
            case PROTEUS_WEDGE18:
                return 18;
            case PROTEUS_WEDGE40:
                return 40;
            case PROTEUS_WEDGE75:
                return 75;
            case PROTEUS_WEDGE126:
                return 126;
            case PROTEUS_WEDGE196:
                return 196;
            case PROTEUS_WEDGE288:
                return 288;
            case PROTEUS_WEDGE405:
                return 405;
            case PROTEUS_WEDGE2601:
                return 2601;
            case PROTEUS_PYRAMID5:
                return 5;
            case PROTEUS_PYRAMID14:
                return 14;
            case PROTEUS_PYRAMID30:
                return 30;
            case PROTEUS_PYRAMID55:
                return 55;
            case PROTEUS_PYRAMID91:
                return 91;
            case PROTEUS_PYRAMID140:
                return 140;
            case PROTEUS_PYRAMID204:
                return 204;
            case PROTEUS_PYRAMID285:
                return 285;
            case PROTEUS_PYRAMID1785:
                return 1785;
            default: {
                SMESH_ERROR("No number of nodes found for type: %s\n", type_to_string(type));
                return 0;
            }
        }
    }

    inline int elem_num_sides(const enum ElemType type) {
        switch (type) {
            case NIL:
                return 0;
            case EDGE2:
            case EDGESHELL2:
            case EDGE3:
            case EDGESHELL3:
                return 2;
            case TRI3:
                return 3;
            case TRISHELL3:
                return 3;
            case MACRO_TRI3:
                return 3;  // Really?
            case MACRO_TET4:
                return 4;
            case QUAD4:
                return 4;
            case QUAD9:
                return 4;
            case QUADSHELL4:
                return 4;
            case QUADSHELL9:
                return 4;
            case TET4:
                return 4;
            case PYRAMID5:
                return 5;
            case WEDGE6:
                return 5;
            case TRI6:
                return 3;
            case HEX8:
                return 6;
            case HEX27:
                return 6;
            case TET10:
                return 4;
            case TET15:
                return 4;
            case TET20:
                return 4;
            case PROTEUS_HEX8:
                return 6;
            case PROTEUS_HEX27:
                return 6;
            case PROTEUS_HEX64:
                return 6;
            case PROTEUS_HEX125:
                return 6;
            case PROTEUS_HEX4913:
                return 6;
            case PROTEUS_QUAD289:
                return 4;
            case PROTEUS_TET4:
            case PROTEUS_TET10:
            case PROTEUS_TET20:
            case PROTEUS_TET35:
            case PROTEUS_TET56:
            case PROTEUS_TET84:
            case PROTEUS_TET120:
            case PROTEUS_TET165:
            case PROTEUS_TET969:
                return 4;
            case PROTEUS_WEDGE6:
            case PROTEUS_WEDGE18:
            case PROTEUS_WEDGE40:
            case PROTEUS_WEDGE75:
            case PROTEUS_WEDGE126:
            case PROTEUS_WEDGE196:
            case PROTEUS_WEDGE288:
            case PROTEUS_WEDGE405:
            case PROTEUS_WEDGE2601:
            case PROTEUS_PYRAMID5:
            case PROTEUS_PYRAMID14:
            case PROTEUS_PYRAMID30:
            case PROTEUS_PYRAMID55:
            case PROTEUS_PYRAMID91:
            case PROTEUS_PYRAMID140:
            case PROTEUS_PYRAMID204:
            case PROTEUS_PYRAMID285:
            case PROTEUS_PYRAMID1785:
                return 5;
            default: {
                SMESH_ERROR("No number of sides found for type: %s\n", type_to_string(type));
                return 0;
            }
        }
    }

    inline int elem_manifold_dim(const enum ElemType type) {
        switch (type) {
            case NIL:
                return 0;
            case EDGE2:
            case EDGESHELL2:
            case EDGE3:
            case EDGESHELL3:
                return 1;
            case TRI3:
                return 2;
            case QUAD4:
                return 2;
            case QUAD9:
                return 2;
            case QUADSHELL4:
                return 2;
            case QUADSHELL9:
                return 2;
            case TET4:
                return 3;
            case PYRAMID5:
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
            case HEX27:
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
            case PROTEUS_HEX4913:
                return 3;
            case PROTEUS_QUAD289:
                return 2;
            case PROTEUS_TET4:
            case PROTEUS_TET10:
            case PROTEUS_TET20:
            case PROTEUS_TET35:
            case PROTEUS_TET56:
            case PROTEUS_TET84:
            case PROTEUS_TET120:
            case PROTEUS_TET165:
            case PROTEUS_TET969:
                return 3;
            case PROTEUS_WEDGE6:
            case PROTEUS_WEDGE18:
            case PROTEUS_WEDGE40:
            case PROTEUS_WEDGE75:
            case PROTEUS_WEDGE126:
            case PROTEUS_WEDGE196:
            case PROTEUS_WEDGE288:
            case PROTEUS_WEDGE405:
            case PROTEUS_WEDGE2601:
            case PROTEUS_PYRAMID5:
            case PROTEUS_PYRAMID14:
            case PROTEUS_PYRAMID30:
            case PROTEUS_PYRAMID55:
            case PROTEUS_PYRAMID91:
            case PROTEUS_PYRAMID140:
            case PROTEUS_PYRAMID204:
            case PROTEUS_PYRAMID285:
            case PROTEUS_PYRAMID1785:
                return 3;
            default: {
                SMESH_ERROR("No manifold dimension found for type: %s\n", type_to_string(type));
                return INVALID;
            }
        }
    }

    /// True if \c overlap is the node count of some side of \c type.
    inline bool elem_overlap_is_full_side(const enum ElemType type, const int overlap) {
        switch (type) {
            case WEDGE6:
            case PYRAMID5:
            case PROTEUS_WEDGE6:
            case PROTEUS_WEDGE18:
            case PROTEUS_WEDGE40:
            case PROTEUS_WEDGE75:
            case PROTEUS_WEDGE126:
            case PROTEUS_WEDGE196:
            case PROTEUS_WEDGE288:
            case PROTEUS_WEDGE405:
            case PROTEUS_WEDGE2601:
            case PROTEUS_PYRAMID5:
            case PROTEUS_PYRAMID14:
            case PROTEUS_PYRAMID30:
            case PROTEUS_PYRAMID55:
            case PROTEUS_PYRAMID91:
            case PROTEUS_PYRAMID140:
            case PROTEUS_PYRAMID204:
            case PROTEUS_PYRAMID285:
            case PROTEUS_PYRAMID1785:
                return overlap == 3 || overlap == 4;
            case TET10:
                return overlap == 3;
            default:
                return overlap == elem_num_nodes(side_type(type));
        }
    }

    inline enum ElemType macro_type_variant(const enum ElemType type) {
        switch (type) {
            case TET10:
                return MACRO_TET4;
            case TRI6:
                return MACRO_TRI3;

            default: {
                SMESH_ERROR("No macro type variant found for type: %s\n", type_to_string(type));
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
            case PROTEUS_HEX8:
                return PROTEUS_HEX8;
            case PROTEUS_HEX27:
                return PROTEUS_HEX8;
            case PROTEUS_HEX64:
                return PROTEUS_HEX8;
            case PROTEUS_HEX125:
                return PROTEUS_HEX8;
            case PROTEUS_HEX216:
                return PROTEUS_HEX8;
            case PROTEUS_HEX343:
                return PROTEUS_HEX8;
            case PROTEUS_HEX512:
                return PROTEUS_HEX8;
            case PROTEUS_HEX729:
                return PROTEUS_HEX8;
            case PROTEUS_HEX4913:
                return PROTEUS_HEX8;
            case PROTEUS_TET4:
                return PROTEUS_TET4;
            case PROTEUS_TET10:
            case PROTEUS_TET20:
            case PROTEUS_TET35:
            case PROTEUS_TET56:
            case PROTEUS_TET84:
            case PROTEUS_TET120:
            case PROTEUS_TET165:
            case PROTEUS_TET969:
                return PROTEUS_TET4;
            default: {
                SMESH_ERROR("No macro base elem found for type: %s\n", type_to_string(macro_type));
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
            case QUAD9:
                return 1;
            case EDGE3:
                return 1;
            case EDGESHELL3:
                return 1;
            case HEX27:
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
            case 16:
                return PROTEUS_HEX4913;
            default:
                SMESH_ERROR("proteus_hex_type:Invalid element setup for proteus hex: %d", micro_elements_per_dim);
                return INVALID;
        }
    }

    /// Macro family for SS conversion. HEX8 / PROTEUS_HEX* map to HEX8.
    /// Do not use \c macro_base_elem on HEX8 (that path errors).
    /// Mixed-family unique numbering (shared HEX–TET edges) cannot be two
    /// independent per-family concats; that is a follow-on generator, not A5.
    inline enum ElemType ss_source_family(const enum ElemType type) {
        switch (type) {
            case HEX8:
            case PROTEUS_HEX8:
            case PROTEUS_HEX27:
            case PROTEUS_HEX64:
            case PROTEUS_HEX125:
            case PROTEUS_HEX216:
            case PROTEUS_HEX343:
            case PROTEUS_HEX512:
            case PROTEUS_HEX729:
            case PROTEUS_HEX4913:
                return HEX8;
            case TET4:
            case MACRO_TET4:
            case PROTEUS_TET4:
            case PROTEUS_TET10:
            case PROTEUS_TET20:
            case PROTEUS_TET35:
            case PROTEUS_TET56:
            case PROTEUS_TET84:
            case PROTEUS_TET120:
            case PROTEUS_TET165:
            case PROTEUS_TET969:
                return TET4;
            case QUAD4:
            case QUADSHELL4:
            case PROTEUS_QUAD4:
            case PROTEUS_QUAD9:
            case PROTEUS_QUAD16:
            case PROTEUS_QUAD25:
            case PROTEUS_QUAD36:
            case PROTEUS_QUAD49:
            case PROTEUS_QUAD64:
            case PROTEUS_QUAD81:
            case PROTEUS_QUAD289:
            case PROTEUS_QUADSHELL4:
            case PROTEUS_QUADSHELL9:
            case PROTEUS_QUADSHELL16:
            case PROTEUS_QUADSHELL25:
            case PROTEUS_QUADSHELL36:
            case PROTEUS_QUADSHELL49:
            case PROTEUS_QUADSHELL64:
            case PROTEUS_QUADSHELL81:
            case PROTEUS_QUADSHELL289:
                return QUAD4;
            case WEDGE6:
            case PROTEUS_WEDGE6:
            case PROTEUS_WEDGE18:
            case PROTEUS_WEDGE40:
            case PROTEUS_WEDGE75:
            case PROTEUS_WEDGE126:
            case PROTEUS_WEDGE196:
            case PROTEUS_WEDGE288:
            case PROTEUS_WEDGE405:
            case PROTEUS_WEDGE2601:
                return WEDGE6;
            case PYRAMID5:
            case PROTEUS_PYRAMID5:
            case PROTEUS_PYRAMID14:
            case PROTEUS_PYRAMID30:
            case PROTEUS_PYRAMID55:
            case PROTEUS_PYRAMID91:
            case PROTEUS_PYRAMID140:
            case PROTEUS_PYRAMID204:
            case PROTEUS_PYRAMID285:
            case PROTEUS_PYRAMID1785:
                return PYRAMID5;
            default:
                return type;
        }
    }

    inline bool is_hex_ss_family(const enum ElemType type) { return ss_source_family(type) == HEX8; }

    inline bool is_tet_ss_family(const enum ElemType type) { return ss_source_family(type) == TET4; }

    inline bool is_quad_ss_family(const enum ElemType type) { return ss_source_family(type) == QUAD4; }

    inline bool is_wedge_ss_family(const enum ElemType type) { return ss_source_family(type) == WEDGE6; }

    inline bool is_pyramid_ss_family(const enum ElemType type) { return ss_source_family(type) == PYRAMID5; }

    inline enum ElemType proteus_tet_type(const int micro_elements_per_dim) {
        switch (micro_elements_per_dim) {
            case 1:
                return PROTEUS_TET4;
            case 2:
                return PROTEUS_TET10;
            case 3:
                return PROTEUS_TET20;
            case 4:
                return PROTEUS_TET35;
            case 5:
                return PROTEUS_TET56;
            case 6:
                return PROTEUS_TET84;
            case 7:
                return PROTEUS_TET120;
            case 8:
                return PROTEUS_TET165;
            case 16:
                return PROTEUS_TET969;
            default:
                SMESH_ERROR("proteus_tet_type: Invalid element setup for proteus tet: %d", micro_elements_per_dim);
                return INVALID;
        }
    }

    inline enum ElemType proteus_quad_type(const int micro_elements_per_dim) {
        switch (micro_elements_per_dim) {
            case 1:
                return PROTEUS_QUAD4;
            case 2:
                return PROTEUS_QUAD9;
            case 3:
                return PROTEUS_QUAD16;
            case 4:
                return PROTEUS_QUAD25;
            case 5:
                return PROTEUS_QUAD36;
            case 6:
                return PROTEUS_QUAD49;
            case 7:
                return PROTEUS_QUAD64;
            case 8:
                return PROTEUS_QUAD81;
            case 16:
                return PROTEUS_QUAD289;
            default:
                SMESH_ERROR("proteus_quad_type: Invalid element setup for proteus quad: %d",
                            micro_elements_per_dim);
                return INVALID;
        }
    }

    inline enum ElemType proteus_wedge_type(const int micro_elements_per_dim) {
        switch (micro_elements_per_dim) {
            case 1:
                return PROTEUS_WEDGE6;
            case 2:
                return PROTEUS_WEDGE18;
            case 3:
                return PROTEUS_WEDGE40;
            case 4:
                return PROTEUS_WEDGE75;
            case 5:
                return PROTEUS_WEDGE126;
            case 6:
                return PROTEUS_WEDGE196;
            case 7:
                return PROTEUS_WEDGE288;
            case 8:
                return PROTEUS_WEDGE405;
            case 16:
                return PROTEUS_WEDGE2601;
            default:
                SMESH_ERROR("proteus_wedge_type: Invalid element setup for proteus wedge: %d",
                            micro_elements_per_dim);
                return INVALID;
        }
    }

    inline enum ElemType proteus_pyramid_type(const int micro_elements_per_dim) {
        switch (micro_elements_per_dim) {
            case 1:
                return PROTEUS_PYRAMID5;
            case 2:
                return PROTEUS_PYRAMID14;
            case 3:
                return PROTEUS_PYRAMID30;
            case 4:
                return PROTEUS_PYRAMID55;
            case 5:
                return PROTEUS_PYRAMID91;
            case 6:
                return PROTEUS_PYRAMID140;
            case 7:
                return PROTEUS_PYRAMID204;
            case 8:
                return PROTEUS_PYRAMID285;
            case 16:
                return PROTEUS_PYRAMID1785;
            default:
                SMESH_ERROR("proteus_pyramid_type: Invalid element setup for proteus pyramid: %d",
                            micro_elements_per_dim);
                return INVALID;
        }
    }

    inline enum ElemType semistructured_type(const enum ElemType type, const int micro_elements_per_dim) {
        switch (type) {
            case HEX8: {
                return proteus_hex_type(micro_elements_per_dim);
            }
            case PROTEUS_HEX8: {
                return proteus_hex_type(micro_elements_per_dim);
            }
            case TET4:
            case PROTEUS_TET4: {
                return proteus_tet_type(micro_elements_per_dim);
            }
            case QUAD4:
            case PROTEUS_QUAD4: {
                return proteus_quad_type(micro_elements_per_dim);
            }
            case QUADSHELL4:
            case PROTEUS_QUADSHELL4: {
                return shell_type(proteus_quad_type(micro_elements_per_dim));
            }
            case WEDGE6:
            case PROTEUS_WEDGE6: {
                return proteus_wedge_type(micro_elements_per_dim);
            }
            case PYRAMID5:
            case PROTEUS_PYRAMID5: {
                return proteus_pyramid_type(micro_elements_per_dim);
            }
            default: {
                SMESH_ERROR("semistructured_type: Invalid element type %d", type);
                return INVALID;
            }
        }
    }

    inline int proteus_hex_micro_elements_per_dim(const enum ElemType type) {
        switch (type) {
            case PROTEUS_HEX8:
                return 1;
            case PROTEUS_HEX27:
                return 2;
            case PROTEUS_HEX64:
                return 3;
            case PROTEUS_HEX125:
                return 4;
            case PROTEUS_HEX216:
                return 5;
            case PROTEUS_HEX343:
                return 6;
            case PROTEUS_HEX512:
                return 7;
            case PROTEUS_HEX729:
                return 8;
            case PROTEUS_HEX4913:
                return 16;
            default: {
                SMESH_ERROR("proteus_hex_micro_elements_per_dim: Invalid element type %d", type);
                return INVALID;
            }
        }
    }

    inline int proteus_tet_micro_elements_per_dim(const enum ElemType type) {
        switch (type) {
            case PROTEUS_TET4:
                return 1;
            case PROTEUS_TET10:
                return 2;
            case PROTEUS_TET20:
                return 3;
            case PROTEUS_TET35:
                return 4;
            case PROTEUS_TET56:
                return 5;
            case PROTEUS_TET84:
                return 6;
            case PROTEUS_TET120:
                return 7;
            case PROTEUS_TET165:
                return 8;
            case PROTEUS_TET969:
                return 16;
            default: {
                SMESH_ERROR("proteus_tet_micro_elements_per_dim: Invalid element type %d", type);
                return INVALID;
            }
        }
    }

    inline int proteus_quad_micro_elements_per_dim(const enum ElemType type) {
        switch (type) {
            case PROTEUS_QUAD4:
                return 1;
            case PROTEUS_QUAD9:
                return 2;
            case PROTEUS_QUAD16:
                return 3;
            case PROTEUS_QUAD25:
                return 4;
            case PROTEUS_QUAD36:
                return 5;
            case PROTEUS_QUAD49:
                return 6;
            case PROTEUS_QUAD64:
                return 7;
            case PROTEUS_QUAD81:
                return 8;
            case PROTEUS_QUAD289:
                return 16;
            case PROTEUS_QUADSHELL4:
                return 1;
            case PROTEUS_QUADSHELL9:
                return 2;
            case PROTEUS_QUADSHELL16:
                return 3;
            case PROTEUS_QUADSHELL25:
                return 4;
            case PROTEUS_QUADSHELL36:
                return 5;
            case PROTEUS_QUADSHELL49:
                return 6;
            case PROTEUS_QUADSHELL64:
                return 7;
            case PROTEUS_QUADSHELL81:
                return 8;
            case PROTEUS_QUADSHELL289:
                return 16;
            default: {
                SMESH_ERROR("proteus_quad_micro_elements_per_dim: Invalid element type %s", type_to_string(type));
                return INVALID;
            }
        }
    }

    inline int proteus_wedge_micro_elements_per_dim(const enum ElemType type) {
        switch (type) {
            case PROTEUS_WEDGE6:
                return 1;
            case PROTEUS_WEDGE18:
                return 2;
            case PROTEUS_WEDGE40:
                return 3;
            case PROTEUS_WEDGE75:
                return 4;
            case PROTEUS_WEDGE126:
                return 5;
            case PROTEUS_WEDGE196:
                return 6;
            case PROTEUS_WEDGE288:
                return 7;
            case PROTEUS_WEDGE405:
                return 8;
            case PROTEUS_WEDGE2601:
                return 16;
            default: {
                SMESH_ERROR("proteus_wedge_micro_elements_per_dim: Invalid element type %s",
                            type_to_string(type));
                return INVALID;
            }
        }
    }

    inline int proteus_pyramid_micro_elements_per_dim(const enum ElemType type) {
        switch (type) {
            case PROTEUS_PYRAMID5:
                return 1;
            case PROTEUS_PYRAMID14:
                return 2;
            case PROTEUS_PYRAMID30:
                return 3;
            case PROTEUS_PYRAMID55:
                return 4;
            case PROTEUS_PYRAMID91:
                return 5;
            case PROTEUS_PYRAMID140:
                return 6;
            case PROTEUS_PYRAMID204:
                return 7;
            case PROTEUS_PYRAMID285:
                return 8;
            case PROTEUS_PYRAMID1785:
                return 16;
            default: {
                SMESH_ERROR("proteus_pyramid_micro_elements_per_dim: Invalid element type %s",
                            type_to_string(type));
                return INVALID;
            }
        }
    }

    /// Semi-structured refinement level encoded in Proteus types (e.g. PROTEUS_HEX27 -> 2).
    inline int semistructured_level(const enum ElemType type) {
        if (is_tet_ss_family(type) && is_semistructured_type(type)) {
            return proteus_tet_micro_elements_per_dim(type);
        }
        if (is_quad_ss_family(type) && is_semistructured_type(type)) {
            return proteus_quad_micro_elements_per_dim(type);
        }
        if (is_wedge_ss_family(type) && is_semistructured_type(type)) {
            return proteus_wedge_micro_elements_per_dim(type);
        }
        if (is_pyramid_ss_family(type) && is_semistructured_type(type)) {
            return proteus_pyramid_micro_elements_per_dim(type);
        }
        return proteus_hex_micro_elements_per_dim(type);
    }

    enum HEX8_Sides { HEX8_LEFT = 3, HEX8_RIGHT = 1, HEX8_BOTTOM = 4, HEX8_TOP = 5, HEX8_FRONT = 0, HEX8_BACK = 2 };

}  // namespace smesh

#endif  // SMESH_ELEM_TYPE_HPP
