#ifndef SMESH_SSQUAD4_PROLONGATION_MATRIX_IMPL_HPP
#define SMESH_SSQUAD4_PROLONGATION_MATRIX_IMPL_HPP

#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_prolongation_matrix.hpp"

namespace smesh {
    template <typename IDXType, typename RowPtrType, typename ColIdxType, typename T>
    int ssquad4_prolongation_matrix_symb(const int                                                 level,
                                         const ptrdiff_t                                           nelements,
                                         const ptrdiff_t                                           nnodes,
                                         const IDXType *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         RowPtrType *const SMESH_RESTRICT                          rowptr) {
        rowptr[0] = 0;

        const int corners[4] = {ssquad4_lidx(level, 0, 0),
                                ssquad4_lidx(level, level, 0),
                                ssquad4_lidx(level, level, level),
                                ssquad4_lidx(level, 0, level)};

#pragma omp parallel for
        for (int d = 0; d < 4; d++) {
            const int c = corners[d];
            for (ptrdiff_t e = 0; e < nelements; e++) {
                const IDXType node = elements[c][e];
                rowptr[node + 1]   = 1;
            }
        }

#pragma omp parallel for
        for (int vi = 1; vi < level; vi++) {
            const int edges[4] = {
                    ssquad4_lidx(level, 0, vi),
                    ssquad4_lidx(level, level, vi),
                    ssquad4_lidx(level, vi, 0),
                    ssquad4_lidx(level, vi, level),
            };

            for (int d = 0; d < 4; d++) {
                const int edge = edges[d];
                for (ptrdiff_t e = 0; e < nelements; e++) {
                    const IDXType node = elements[edge][e];
                    rowptr[node + 1]   = 2;
                }
            }
        }

#pragma omp parallel for collapse(2)
        for (int yi = 1; yi < level; yi++) {
            for (int xi = 1; xi < level; xi++) {
                const int ii = ssquad4_lidx(level, xi, yi);
                for (ptrdiff_t e = 0; e < nelements; e++) {
                    const IDXType node = elements[ii][e];
                    rowptr[node + 1]   = 4;
                }
            }
        }

        for (ptrdiff_t i = 0; i < nnodes; i++) {
            SMESH_ASSERT(rowptr[i + 1] != 0);
            rowptr[i + 1] += rowptr[i];
        }

        return SMESH_SUCCESS;
    }

    template <typename IDXType, typename RowPtrType, typename ColIdxType, typename T>
    int ssquad4_prolongation_matrix_fill(const int                                                 level,
                                         const ptrdiff_t                                           nelements,
                                         const ptrdiff_t                                           nnodes,
                                         const IDXType *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                                         const RowPtrType *const SMESH_RESTRICT                    rowptr,
                                         ColIdxType *const SMESH_RESTRICT                          colidx,
                                         T *const SMESH_RESTRICT                                   values) {
        (void)nnodes;

        const int corners[4] = {ssquad4_lidx(level, 0, 0),
                                ssquad4_lidx(level, level, 0),
                                ssquad4_lidx(level, level, level),
                                ssquad4_lidx(level, 0, level)};

#pragma omp parallel for
        for (int d = 0; d < 4; d++) {
            const int c = corners[d];
            for (ptrdiff_t e = 0; e < nelements; e++) {
                const IDXType    node = elements[c][e];
                const RowPtrType k    = rowptr[node];
                colidx[k]             = node;
                values[k]             = T(1);
            }
        }

        const T h = T(1) / level;

        const int ev0[4] = {corners[0], corners[1], corners[0], corners[3]};
        const int ev1[4] = {corners[3], corners[2], corners[1], corners[2]};

        const auto sort2 = [](IDXType &ca, IDXType &cb, T &wa, T &wb) {
            if (cb < ca) {
                const IDXType temp_c = ca;
                ca                   = cb;
                cb                   = temp_c;

                const T temp_w = wa;
                wa             = wb;
                wb             = temp_w;
            }
        };

#pragma omp parallel for
        for (int vi = 1; vi < level; vi++) {
            const int edges[4] = {
                    ssquad4_lidx(level, 0, vi),
                    ssquad4_lidx(level, level, vi),
                    ssquad4_lidx(level, vi, 0),
                    ssquad4_lidx(level, vi, level),
            };

            const T t  = vi * h;
            const T w0 = T(1) - t;
            const T w1 = t;

            for (int d = 0; d < 4; d++) {
                const int eidx = edges[d];
                const int v0   = ev0[d];
                const int v1   = ev1[d];

                for (ptrdiff_t e = 0; e < nelements; e++) {
                    const IDXType    node = elements[eidx][e];
                    const RowPtrType k    = rowptr[node];
                    IDXType          node0 = elements[v0][e];
                    IDXType          node1 = elements[v1][e];
                    T                w0i   = w0;
                    T                w1i   = w1;

                    if (node1 < node0) {
                        const IDXType temp_node = node0;
                        node0                   = node1;
                        node1                   = temp_node;

                        const T temp_weight = w0i;
                        w0i                 = w1i;
                        w1i                 = temp_weight;
                    }

                    colidx[k]     = node0;
                    values[k]     = w0i;
                    colidx[k + 1] = node1;
                    values[k + 1] = w1i;
                }
            }
        }

#pragma omp parallel for collapse(2)
        for (int yi = 1; yi < level; yi++) {
            for (int xi = 1; xi < level; xi++) {
                const int ii = ssquad4_lidx(level, xi, yi);

                const T x = xi * h;
                const T y = yi * h;
                const T w[4] = {
                        (T(1) - x) * (T(1) - y),
                        x * (T(1) - y),
                        x * y,
                        (T(1) - x) * y,
                };

                for (ptrdiff_t e = 0; e < nelements; e++) {
                    const IDXType    node = elements[ii][e];
                    const RowPtrType k    = rowptr[node];

                    IDXType c0 = elements[corners[0]][e];
                    IDXType c1 = elements[corners[1]][e];
                    IDXType c2 = elements[corners[2]][e];
                    IDXType c3 = elements[corners[3]][e];
                    T       w0 = w[0];
                    T       w1 = w[1];
                    T       w2 = w[2];
                    T       w3 = w[3];

                    sort2(c0, c1, w0, w1);
                    sort2(c2, c3, w2, w3);
                    sort2(c0, c2, w0, w2);
                    sort2(c1, c3, w1, w3);
                    sort2(c1, c2, w1, w2);

                    colidx[k]     = c0;
                    values[k]     = w0;
                    colidx[k + 1] = c1;
                    values[k + 1] = w1;
                    colidx[k + 2] = c2;
                    values[k + 2] = w2;
                    colidx[k + 3] = c3;
                    values[k + 3] = w3;
                }
            }
        }

        return SMESH_SUCCESS;
    }
}  // namespace smesh

#endif  // SMESH_SSQUAD4_PROLONGATION_MATRIX_IMPL_HPP
