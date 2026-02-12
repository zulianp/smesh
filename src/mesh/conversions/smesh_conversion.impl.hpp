#ifndef SMESH_CONVERSION_IMPL_HPP
#define SMESH_CONVERSION_IMPL_HPP

#include "smesh_conversion.hpp"

namespace smesh {

template <typename idx_t>
void mesh_hex8_to_6x_tet4(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements) {
  // #pragma omp parallel for
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    idx_t ii[8];
    for (ptrdiff_t d = 0; d < 8; d++) {
      ii[d] = hex8_elements[d][e];
    }

    idx_t hex8[6][4] = {
        {ii[0], ii[1], ii[3], ii[7]}, {ii[0], ii[1], ii[7], ii[5]},
        {ii[0], ii[4], ii[5], ii[7]}, {ii[1], ii[2], ii[3], ii[6]},
        {ii[1], ii[3], ii[7], ii[6]}, {ii[1], ii[5], ii[6], ii[7]}};

    for (int sub_e = 0; sub_e < 6; sub_e++) {
      for (int node = 0; node < 4; node++) {
        tet4_elements[node][e * 6 + sub_e] = hex8[sub_e][node];
      }
    }
  }
}

template <typename idx_t>
void mesh_tet15_to_4x_hex8(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet15_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements) {
  // #pragma omp parallel for
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    idx_t ii[15];
    for (ptrdiff_t d = 0; d < 15; d++) {
      ii[d] = tet15_elements[d][e];
    }

    idx_t hex8[4][8] = {
        // HEX8(0)
        {ii[0], ii[4], ii[13], ii[6], ii[7], ii[10], ii[14], ii[12]},
        // HEX8(1)
        {ii[4], ii[1], ii[5], ii[13], ii[10], ii[8], ii[11], ii[14]},
        // HEX8(2)
        {ii[13], ii[5], ii[2], ii[6], ii[14], ii[11], ii[9], ii[12]},
        // HEX8(3)
        {ii[7], ii[10], ii[14], ii[12], ii[3], ii[8], ii[11], ii[9]}};

    for (int sub_e = 0; sub_e < 4; sub_e++) {
      for (int node = 0; node < 8; node++) {
        hex8_elements[node][e * 4 + sub_e] = hex8[sub_e][node];
      }
    }
  }
}
template <typename idx_t>
void mesh_wedge6_to_3x_tet4(
    const ptrdiff_t n_elements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT wedge6_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT tet4_elements) {
  for (ptrdiff_t e = 0; e < n_elements; e++) {
    const idx_t i0 = wedge6_elements[0][2];
    const idx_t i1 = wedge6_elements[1][2];
    const idx_t i2 = wedge6_elements[2][2];
    const idx_t i3 = wedge6_elements[3][2];
    const idx_t i4 = wedge6_elements[4][2];
    const idx_t i5 = wedge6_elements[5][2];

    idx_t *node0 = &tet4_elements[0][e * 3];
    idx_t *node1 = &tet4_elements[1][e * 3];
    idx_t *node2 = &tet4_elements[2][e * 3];
    idx_t *node3 = &tet4_elements[3][e * 3];

    // tet 0
    node0[0] = i0;
    node1[0] = i1;
    node2[0] = i2;
    node3[0] = i3;

    // tet 1
    node0[1] = i1;
    node1[1] = i4;
    node2[1] = i5;
    node3[1] = i3;

    // tet 2
    node0[2] = i2;
    node1[2] = i1;
    node2[2] = i5;
    node3[2] = i3;
  }
}

} // namespace smesh

#endif // SMESH_CONVERSION_IMPL_HPP
