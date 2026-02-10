#ifndef SMESH_PROLONGATION_HPP
#define SMESH_PROLONGATION_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {
template <typename idx_t, typename T>
int hierarchical_prolongation(
    const enum ElemType from_element, const enum ElemType to_element,
    const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const int vec_size, const T *const SMESH_RESTRICT from,
    T *const SMESH_RESTRICT to);

} // namespace smesh

#endif // SMESH_PROLONGATION_HPP
