#ifndef SMESH_EXTRACT_SHAP_FEATURES_HPP
#define SMESH_EXTRACT_SHAP_FEATURES_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {
    // TODO: add declarations from the impl.hpp file
    template <typename idx_t, typename geom_t, typename count_t>
    int extract_sharp_edges(
        const enum ElemType element_type, const ptrdiff_t nelements,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
        const ptrdiff_t nnodes,
        geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
        const count_t *const SMESH_RESTRICT rowptr,
        const idx_t *const SMESH_RESTRICT colidx, const geom_t angle_threshold,
        ptrdiff_t *out_n_sharp_edges, idx_t **out_e0, idx_t **out_e1);

    template <typename idx_t, typename count_t>
    int extract_sharp_corners(const ptrdiff_t nnodes, const ptrdiff_t n_sharp_edges,
                              idx_t *const SMESH_RESTRICT e0,
                              idx_t *const SMESH_RESTRICT e1,
                              ptrdiff_t *out_ncorners, idx_t **out_corners,
                              int edge_clean_up);

    template <typename idx_t, typename count_t, typename element_idx_t>
    int extract_disconnected_faces(
        const enum ElemType element_type, const ptrdiff_t nelements,
        const ptrdiff_t nnodes,
        idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
        const ptrdiff_t n_sharp_edges, const idx_t *const SMESH_RESTRICT e0,
        const idx_t *const SMESH_RESTRICT e1,
        ptrdiff_t *out_n_disconnected_elements,
        element_idx_t **out_disconnected_elements);
} // namespace smesh

#endif // SMESH_EXTRACT_SHAP_FEATURES_HPP