#ifndef SMESH_SSTET4_TRANSFER_IMPL_HPP
#define SMESH_SSTET4_TRANSFER_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_sstet4.hpp"

namespace smesh {
namespace sstet4_transfer {

    static SMESH_INLINE bool is_power_of_two(const int x) { return x > 0 && (x & (x - 1)) == 0; }

    static SMESH_INLINE void fill_lattice_coords(const int L, int *const x, int *const y, int *const z) {
        for (int zi = 0; zi <= L; ++zi) {
            for (int yi = 0; yi <= L - zi; ++yi) {
                for (int xi = 0; xi <= L - zi - yi; ++xi) {
                    const int lidx = sstet4_lidx(L, xi, yi, zi);
                    x[lidx]        = xi;
                    y[lidx]        = yi;
                    z[lidx]        = zi;
                }
            }
        }
    }

    static SMESH_INLINE void midpoint(const int *const a, const int *const b, int *const m) {
        m[0] = (a[0] + b[0]) / 2;
        m[1] = (a[1] + b[1]) / 2;
        m[2] = (a[2] + b[2]) / 2;
    }

    template <typename Func>
    static SMESH_INLINE void for_each_red_microtet_rec(const int global_L,
                                                       const int edge_L,
                                                       const int *const p0,
                                                       const int *const p1,
                                                       const int *const p2,
                                                       const int *const p3,
                                                       Func &func) {
        if (edge_L == 1) {
            int ev[4] = {sstet4_lidx(global_L, p0[0], p0[1], p0[2]),
                         sstet4_lidx(global_L, p1[0], p1[1], p1[2]),
                         sstet4_lidx(global_L, p2[0], p2[1], p2[2]),
                         sstet4_lidx(global_L, p3[0], p3[1], p3[2])};
            func(ev);
            return;
        }

        const int half_L = edge_L / 2;
        int       m01[3], m12[3], m02[3], m03[3], m13[3], m23[3];
        midpoint(p0, p1, m01);
        midpoint(p1, p2, m12);
        midpoint(p0, p2, m02);
        midpoint(p0, p3, m03);
        midpoint(p1, p3, m13);
        midpoint(p2, p3, m23);

        const int *const pts[10] = {p0, p1, p2, p3, m01, m12, m02, m03, m13, m23};
        static const int pattern[8][4] = {{0, 4, 6, 7},
                                          {4, 1, 5, 8},
                                          {6, 5, 2, 9},
                                          {7, 8, 9, 3},
                                          {4, 5, 6, 8},
                                          {7, 4, 6, 8},
                                          {6, 5, 9, 8},
                                          {7, 6, 9, 8}};

        for (int e = 0; e < 8; ++e) {
            for_each_red_microtet_rec(global_L,
                                      half_L,
                                      pts[pattern[e][0]],
                                      pts[pattern[e][1]],
                                      pts[pattern[e][2]],
                                      pts[pattern[e][3]],
                                      func);
        }
    }

    template <typename Func>
    static SMESH_INLINE void for_each_red_microtet(const int L, Func &func) {
        const int p0[3] = {0, 0, 0};
        const int p1[3] = {L, 0, 0};
        const int p2[3] = {0, L, 0};
        const int p3[3] = {0, 0, L};
        for_each_red_microtet_rec(L, L, p0, p1, p2, p3, func);
    }

    template <typename Func>
    static SMESH_INLINE void for_each_microtet(const int L, Func &&func) {
        int ev[4];

        if (is_power_of_two(L)) {
            for_each_red_microtet(L, func);
            return;
        }

        const int n = L + 1;

        if (L == 1) {
            ev[0] = sstet4_lidx(1, 0, 0, 0);
            ev[1] = sstet4_lidx(1, 1, 0, 0);
            ev[2] = sstet4_lidx(1, 0, 1, 0);
            ev[3] = sstet4_lidx(1, 0, 0, 1);
            func(ev);
            return;
        }

        {
            int p = 0;
            for (int i = 0; i < n - 1; i++) {
                const int layer_items = (n - i + 1) * (n - i) / 2;
                for (int j = 0; j < n - i - 1; j++) {
                    for (int k = 0; k < n - i - j - 1; k++) {
                        ev[0] = p;
                        ev[1] = p + 1;
                        ev[2] = p + n - i - j;
                        ev[3] = p + layer_items - j;
                        func(ev);
                        p++;
                    }
                    p++;
                }
                p++;
            }
        }

        {
            int p = 0;
            for (int i = 0; i < n - 1; i++) {
                const int layer_items = (n - i) * (n - i - 1) / 2;
                for (int j = 0; j < n - i - 1; j++) {
                    p++;
                    for (int k = 1; k < n - i - j - 1; k++) {
                        ev[0] = p;
                        ev[1] = p + layer_items + n - i - j - 1;
                        ev[2] = p + layer_items + n - i - j;
                        ev[3] = p + layer_items + n - i - j - 1 + n - i - j - 1;
                        func(ev);
                        p++;
                    }
                    p++;
                }
                p++;
            }
        }

        {
            int p = 0;
            for (int i = 0; i < n - 1; i++) {
                const int layer_items = (n - i) * (n - i - 1) / 2;
                for (int j = 0; j < n - i - 1; j++) {
                    p++;
                    for (int k = 1; k < n - i - j - 1; k++) {
                        ev[0] = p;
                        ev[1] = p + n - i - j;
                        ev[2] = p + layer_items + n - i - j - 1 + n - i - j - 1;
                        ev[3] = p + layer_items + n - i - j;
                        func(ev);
                        p++;
                    }
                    p++;
                }
                p++;
            }
        }

        {
            int p = 0;
            for (int i = 0; i < n - 1; i++) {
                const int layer_items = (n - i) * (n - i - 1) / 2;
                for (int j = 0; j < n - i - 1; j++) {
                    p++;
                    for (int k = 1; k < n - i - j - 1; k++) {
                        ev[0] = p;
                        ev[1] = p + n - i - j - 1;
                        ev[2] = p + layer_items + n - i - j - 1;
                        ev[3] = p + layer_items + n - i - j - 1 + n - i - j - 1;
                        func(ev);
                        p++;
                    }
                    p++;
                }
                p++;
            }
        }

        {
            int p = 0;
            for (int i = 1; i < n - 1; i++) {
                p += n - i + 1;
                const int layer_items = (n - i) * (n - i - 1) / 2;
                for (int j = 0; j < n - i - 1; j++) {
                    p++;
                    for (int k = 1; k < n - i - j - 1; k++) {
                        ev[0] = p;
                        ev[1] = p + layer_items + n - i;
                        ev[2] = p + layer_items + n - i - j + n - i;
                        ev[3] = p + layer_items + n - i - j + n - i - 1;
                        func(ev);
                        p++;
                    }
                    p++;
                }
                p++;
            }
        }

        {
            int p = 0;
            for (int i = 0; i < n - 1; i++) {
                const int layer_items = (n - i) * (n - i - 1) / 2;
                for (int j = 0; j < n - i - 1; j++) {
                    p++;
                    for (int k = 1; k < n - i - j - 1; k++) {
                        ev[0] = p;
                        ev[1] = p + n - i - j - 1;
                        ev[2] = p + layer_items + n - i - j - 1 + n - i - j - 1;
                        ev[3] = p + n - i - j;
                        func(ev);
                        p++;
                    }
                    p++;
                }
                p++;
            }
        }
    }

    template <typename T>
    static SMESH_INLINE T det3(const T ax,
                               const T ay,
                               const T az,
                               const T bx,
                               const T by,
                               const T bz,
                               const T cx,
                               const T cy,
                               const T cz) {
        return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
    }

    template <typename T>
    static SMESH_INLINE bool barycentric(const T px,
                                         const T py,
                                         const T pz,
                                         const int *const x,
                                         const int *const y,
                                         const int *const z,
                                         const int *const ev,
                                         T *const w) {
        const T x0 = static_cast<T>(x[ev[0]]);
        const T y0 = static_cast<T>(y[ev[0]]);
        const T z0 = static_cast<T>(z[ev[0]]);
        const T ax = static_cast<T>(x[ev[1]]) - x0;
        const T ay = static_cast<T>(y[ev[1]]) - y0;
        const T az = static_cast<T>(z[ev[1]]) - z0;
        const T bx = static_cast<T>(x[ev[2]]) - x0;
        const T by = static_cast<T>(y[ev[2]]) - y0;
        const T bz = static_cast<T>(z[ev[2]]) - z0;
        const T cx = static_cast<T>(x[ev[3]]) - x0;
        const T cy = static_cast<T>(y[ev[3]]) - y0;
        const T cz = static_cast<T>(z[ev[3]]) - z0;
        const T qx = px - x0;
        const T qy = py - y0;
        const T qz = pz - z0;

        const T den = det3(ax, ay, az, bx, by, bz, cx, cy, cz);
        if (den == T(0)) {
            return false;
        }

        w[1] = det3(qx, qy, qz, bx, by, bz, cx, cy, cz) / den;
        w[2] = det3(ax, ay, az, qx, qy, qz, cx, cy, cz) / den;
        w[3] = det3(ax, ay, az, bx, by, bz, qx, qy, qz) / den;
        w[0] = T(1) - w[1] - w[2] - w[3];

        const T tol = T(1e-10);
        return w[0] >= -tol && w[1] >= -tol && w[2] >= -tol && w[3] >= -tol && w[0] <= T(1) + tol &&
               w[1] <= T(1) + tol && w[2] <= T(1) + tol && w[3] <= T(1) + tol;
    }

    template <typename T>
    static SMESH_INLINE int find_weights(const int      L,
                                         const T        px,
                                         const T        py,
                                         const T        pz,
                                         const int *const x,
                                         const int *const y,
                                         const int *const z,
                                         int *const      coarse_lidx,
                                         T *const        weights) {
        int found = 0;
        for_each_microtet(L, [&](const int *const ev) {
            if (!found && barycentric(px, py, pz, x, y, z, ev, weights)) {
                coarse_lidx[0] = ev[0];
                coarse_lidx[1] = ev[1];
                coarse_lidx[2] = ev[2];
                coarse_lidx[3] = ev[3];
                found          = 1;
            }
        });
        return found ? SMESH_SUCCESS : SMESH_FAILURE;
    }

}  // namespace sstet4_transfer
}  // namespace smesh

#endif  // SMESH_SSTET4_TRANSFER_IMPL_HPP
