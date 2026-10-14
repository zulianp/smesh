#ifndef SMESH_JACOBIANS_HPP
#define SMESH_JACOBIANS_HPP

#include "smesh_elem_type.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_types.hpp"

#include <stddef.h>

namespace smesh {

    /// Packed adjugate SoA length: 4 for planar TRI3/QUAD4, 9 for 3D families.
    inline int jacobian_adjugate_components(const enum ElemType element_type) {
        switch (element_type) {
            case TRI3:
            case QUAD4:
                return 4;
            default:
                return 9;
        }
    }

    template <typename AdjugateType, typename DeterminantType>
    int tri3_adjugate_fill(const ptrdiff_t                                          nelements,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                           const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                           const ptrdiff_t                                          stride,
                           AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
                           DeterminantType *const SMESH_RESTRICT                    determinant);

    template <typename AdjugateType, typename DeterminantType>
    int quad4_adjugate_fill(const ptrdiff_t                                          nelements,
                            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                            const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                            const geom_t                                             qx,
                            const geom_t                                             qy,
                            const ptrdiff_t                                          stride,
                            AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
                            DeterminantType *const SMESH_RESTRICT                    determinant);

    template <typename AdjugateType, typename DeterminantType>
    int tet4_adjugate_fill(const ptrdiff_t                                          nelements,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                           const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                           const ptrdiff_t                                          stride,
                           AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
                           DeterminantType *const SMESH_RESTRICT                    determinant);

    template <typename AdjugateType, typename DeterminantType>
    int hex8_adjugate_fill(const ptrdiff_t                                          nelements,
                           const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                           const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                           const geom_t                                             qx,
                           const geom_t                                             qy,
                           const geom_t                                             qz,
                           const ptrdiff_t                                          stride,
                           AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
                           DeterminantType *const SMESH_RESTRICT                    determinant);

    /// Per-macro-element adjugate / det J at reference point (1/2,1/2,1/2) from SSHEX8 macro corners.
    template <typename AdjugateType, typename DeterminantType>
    int sshex8_macro_adjugate_fill(const int                                                level,
                                   const ptrdiff_t                                          nelements,
                                   const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                                   const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                                   const ptrdiff_t                                          stride,
                                   AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
                                   DeterminantType *const SMESH_RESTRICT                    determinant);

    template <typename AdjugateType, typename DeterminantType>
    int adjugate_fill(const enum ElemType                                      element_type,
                      const ptrdiff_t                                          nelements,
                      const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                      const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                      const ptrdiff_t                                          stride,
                      AdjugateType *const SMESH_RESTRICT *const SMESH_RESTRICT adjugate,
                      DeterminantType *const SMESH_RESTRICT                    determinant) {
        if (is_semistructured_type(element_type) && is_quad_ss_family(element_type)) {
            // The macro quadrilateral through its four corner nodes, as the hex
            // family below does with its eight.
            const int L = proteus_quad_micro_elements_per_dim(element_type);
            const idx_t *const SMESH_RESTRICT macro_elements[4] = {elements[ssquad4_lidx(L, 0, 0)],
                                                                   elements[ssquad4_lidx(L, L, 0)],
                                                                   elements[ssquad4_lidx(L, L, L)],
                                                                   elements[ssquad4_lidx(L, 0, L)]};
            return quad4_adjugate_fill(nelements, macro_elements, points, 0.5, 0.5, stride, adjugate, determinant);
        }

        if (is_semistructured_type(element_type)) {
            return sshex8_macro_adjugate_fill(
                    proteus_hex_micro_elements_per_dim(element_type), nelements, elements, points, stride, adjugate, determinant);
        }

        switch (element_type) {
            case TRI3:
                return tri3_adjugate_fill(nelements, elements, points, stride, adjugate, determinant);
            case QUAD4:
                return quad4_adjugate_fill(nelements, elements, points, 0.5, 0.5, stride, adjugate, determinant);
            case TET4:
            case TET10:
                return tet4_adjugate_fill(nelements, elements, points, stride, adjugate, determinant);
            case HEX8:
                return hex8_adjugate_fill(nelements, elements, points, 0.5, 0.5, 0.5, stride, adjugate, determinant);
            default:
                SMESH_ERROR("adjugate_fill: Unsupported element type: %s\n", type_to_string(element_type));
                return SMESH_FAILURE;
        }
    }
}  // namespace smesh

#endif  // SMESH_JACOBIANS_HPP
