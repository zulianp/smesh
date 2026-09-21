#ifndef SMESH_HEX8_NOZZLE_GRID_HPP
#define SMESH_HEX8_NOZZLE_GRID_HPP

#include "smesh_base.hpp"

#include <cmath>
#include <vector>

namespace smesh {

template <typename Scalar>
inline bool nozzle_arguments_valid(const std::vector<Scalar> &x_breaks, const std::vector<Scalar> &bore_radius,
                                   const std::vector<ptrdiff_t> &n_axial, const ptrdiff_t expansion,
                                   const Scalar expanded_radius, const ptrdiff_t n_core, const ptrdiff_t n_bore,
                                   const ptrdiff_t n_outer, const Scalar core_fraction,
                                   const Scalar radial_grading, const Scalar axial_grading,
                                   const char *const who) {
    const ptrdiff_t n_segments = (ptrdiff_t)n_axial.size();
    if (n_segments < 1 || (ptrdiff_t)x_breaks.size() != n_segments + 1 ||
        (ptrdiff_t)bore_radius.size() != n_segments + 1) {
        SMESH_ERROR("%s: need S >= 1 segments, S+1 breaks and S+1 radii", who);
        return false;
    }
    for (ptrdiff_t s = 0; s < n_segments; ++s) {
        if (n_axial[(size_t)s] < 1 || !(x_breaks[(size_t)s + 1] > x_breaks[(size_t)s])) {
            SMESH_ERROR("%s: segment %td needs n_axial >= 1 and increasing x", who, s);
            return false;
        }
    }
    for (ptrdiff_t k = 0; k <= n_segments; ++k) {
        if (!(bore_radius[(size_t)k] > 0)) {
            SMESH_ERROR("%s: bore radius must be positive at break %td", who, k);
            return false;
        }
    }
    if (n_core < 2 || n_core % 2 != 0 || n_bore < 1 || !(core_fraction > 0) || !(core_fraction < 1)) {
        SMESH_ERROR("%s: need an even n_core >= 2, n_bore >= 1 and "
                    "0 < core_fraction < 1",
                    who);
        return false;
    }
    // Zero is ungraded and is the default; negative would invert the stretch and fold the
    // layers it is meant to cluster. The upper bound is not a taste: tanh(8) is 1 - 2.3e-7, so
    // beyond it every interior layer lands on the wall to within single precision and the mesh
    // degenerates silently instead of failing.
    if (!(radial_grading >= 0) || !(axial_grading >= 0) || radial_grading > 8 || axial_grading > 8) {
        SMESH_ERROR("%s: grading must lie in [0, 8]; got radial %g, axial %g", who,
                    (double)radial_grading, (double)axial_grading);
        return false;
    }
    const bool has_expansion = expansion >= 0 && expansion < n_segments;
    if (has_expansion) {
        if (n_outer < 1) {
            SMESH_ERROR("%s: an expansion needs n_outer >= 1", who);
            return false;
        }
        for (ptrdiff_t k = expansion; k <= n_segments; ++k) {
            if (!(bore_radius[(size_t)k] < expanded_radius)) {
                SMESH_ERROR("%s: downstream of the expansion the bore "
                            "(break %td) must be narrower than the expanded radius",
                            who, k);
                return false;
            }
        }
    }
    return true;
}

// A one-sided tanh stretch of a normalised coordinate, clustering cells toward t = 1. On the
// ring layers that puts resolution at the bore wall; on the axial planes it puts it at the
// downstream end of each segment, which for the FDA nozzle is the throat and the expansion --
// the two places the gradients are.
//
// beta == 0 RETURNS t, as the first statement and before any arithmetic. That is not a
// convenience: it is what makes an ungraded mesh bit-for-bit the mesh this generator produced
// before grading existed. A form that merely tends to the identity as beta falls --
// t * (1 + beta * f(t)), say -- would move every node by a rounding error, and every number
// ever recorded on this geometry with it.
//
// The map fixes both ends exactly (tanh(0) = 0, and the ratio is 1 at t = 1) and is monotone
// for beta > 0, so it cannot fold a cell or reorder the grid. Vinokur's one-sided form
// (J. Comput. Phys. 50, 1983); the two-sided variant is not needed because each axial segment
// already ends at a break.
inline double nozzle_stretch(const double t, const double beta) {
    if (beta == 0) return t;
    return std::tanh(beta * t) / std::tanh(beta);
}

struct NozzleSection {
    ptrdiff_t n{0}, n_bore{0}, n_outer{0}, n_perim{0}, n_loops{0}, n_core_nodes{0}, n2d{0}, n2d_inner{0};
    double    c{0};
    double    beta_r{0};  // radial grading toward the bore wall; 0 is ungraded, exactly

    NozzleSection(const ptrdiff_t n_core, const ptrdiff_t nb, const ptrdiff_t no, const bool has_expansion,
                  const double core_fraction, const double radial_grading = 0)
        : n(n_core),
          n_bore(nb),
          n_outer(no),
          n_perim(4 * n_core),
          n_loops(nb + (has_expansion ? no : 0)),
          n_core_nodes((n_core + 1) * (n_core + 1)),
          n2d((n_core + 1) * (n_core + 1) + (nb + (has_expansion ? no : 0)) * 4 * n_core),
          n2d_inner((n_core + 1) * (n_core + 1) + nb * 4 * n_core),
          c(core_fraction),
          beta_r(radial_grading) {}

    ptrdiff_t core_node(const ptrdiff_t i, const ptrdiff_t j) const { return i + j * (n + 1); }

    void perim_ij(const ptrdiff_t m, ptrdiff_t &i, ptrdiff_t &j) const {
        if (m < n) {
            i = n;
            j = m;
        } else if (m < 2 * n) {
            i = 2 * n - m;
            j = n;
        } else if (m < 3 * n) {
            i = 0;
            j = 3 * n - m;
        } else {
            i = m - 3 * n;
            j = 0;
        }
    }

    ptrdiff_t loop_node(const ptrdiff_t l, const ptrdiff_t m) const {
        const ptrdiff_t mm = ((m % n_perim) + n_perim) % n_perim;
        if (l == 0) {
            ptrdiff_t i, j;
            perim_ij(mm, i, j);
            return core_node(i, j);
        }
        return n_core_nodes + (l - 1) * n_perim + mm;
    }

    void core(const double i, const double j, double &by, double &bz, double &oy, double &oz) const {
        by = c * (-1.0 + 2.0 * i / (double)n);
        bz = c * (-1.0 + 2.0 * j / (double)n);
        oy = 0;
        oz = 0;
    }

    void ring(const double l, const double m, const bool outer, double &by, double &bz, double &oy, double &oz) const {
        const double dn = (double)n;
        double       mm = std::fmod(m, (double)n_perim);
        if (mm < 0)
            mm += (double)n_perim;
        double si, sj;
        if (mm < dn) {
            si = dn;
            sj = mm;
        } else if (mm < 2 * dn) {
            si = 2 * dn - mm;
            sj = dn;
        } else if (mm < 3 * dn) {
            si = 0;
            sj = 3 * dn - mm;
        } else {
            si = mm - 3 * dn;
            sj = 0;
        }
        const double sq_y = -1.0 + 2.0 * si / dn;
        const double sq_z = -1.0 + 2.0 * sj / dn;
        const double th   = -0.25 * M_PI + 0.5 * M_PI * m / dn;
        const double cy = std::cos(th), cz = std::sin(th);
        // Graded in s -- the normalised position ACROSS the layers -- and not in l. The stretch
        // fixes s = 0 and s = 1, so the core boundary stays on the core square and the outermost
        // layer stays on the bore; only the interior layers move. That is why core() needs no
        // grading of its own, and why a node on the core boundary is still placed identically
        // whichever element places it.
        if (!outer) {
            const double s = nozzle_stretch(l / (double)n_bore, beta_r);
            by             = (1.0 - s) * c * sq_y + s * cy;
            bz             = (1.0 - s) * c * sq_z + s * cz;
            oy             = 0;
            oz             = 0;
        } else {
            const double s = nozzle_stretch((l - (double)n_bore) / (double)n_outer, beta_r);
            by             = cy;
            bz             = cz;
            oy             = s * cy;
            oz             = s * cz;
        }
    }

    void reference(std::vector<double> &by, std::vector<double> &bz, std::vector<double> &oy,
                   std::vector<double> &oz) const {
        by.assign((size_t)n2d, 0);
        bz.assign((size_t)n2d, 0);
        oy.assign((size_t)n2d, 0);
        oz.assign((size_t)n2d, 0);
        for (ptrdiff_t j = 0; j <= n; ++j)
            for (ptrdiff_t i = 0; i <= n; ++i) {
                const ptrdiff_t v = core_node(i, j);
                core((double)i, (double)j, by[(size_t)v], bz[(size_t)v], oy[(size_t)v], oz[(size_t)v]);
            }
        for (ptrdiff_t l = 1; l <= n_loops; ++l)
            for (ptrdiff_t m = 0; m < n_perim; ++m) {
                const ptrdiff_t v = loop_node(l, m);
                ring((double)l, (double)m, l > n_bore, by[(size_t)v], bz[(size_t)v], oy[(size_t)v], oz[(size_t)v]);
            }
    }
};

template <typename Scalar>
struct NozzlePlanes {
    ptrdiff_t              n_cells_x{0}, n_planes{0};
    std::vector<double>    plane_x, plane_rb;
    std::vector<ptrdiff_t> cell_segment;
    // The segment table, kept rather than discarded after the constructor, so that at() can
    // evaluate the same map at a FRACTIONAL plane coordinate.
    std::vector<double>    seg_x0, seg_x1, seg_rb0, seg_rb1;
    std::vector<ptrdiff_t> seg_first_cell, seg_cells;
    double                 beta_a{0};  // axial grading; 0 is ungraded, exactly

    NozzlePlanes(const std::vector<Scalar> &x_breaks, const std::vector<Scalar> &bore_radius,
                 const std::vector<ptrdiff_t> &n_axial, const double axial_grading = 0)
        : beta_a(axial_grading) {
        const ptrdiff_t n_segments = (ptrdiff_t)n_axial.size();
        for (auto na : n_axial)
            n_cells_x += na;
        n_planes = n_cells_x + 1;
        plane_x.resize((size_t)n_planes);
        plane_rb.resize((size_t)n_planes);
        cell_segment.resize((size_t)n_cells_x);
        ptrdiff_t k = 0;
        for (ptrdiff_t s = 0; s < n_segments; ++s) {
            const ptrdiff_t na = n_axial[(size_t)s];
            seg_first_cell.push_back(k);
            seg_cells.push_back(na);
            seg_x0.push_back((double)x_breaks[(size_t)s]);
            seg_x1.push_back((double)x_breaks[(size_t)s + 1]);
            seg_rb0.push_back((double)bore_radius[(size_t)s]);
            seg_rb1.push_back((double)bore_radius[(size_t)s + 1]);
            for (ptrdiff_t a = 0; a < na; ++a, ++k) {
                const double t          = nozzle_stretch((double)a / (double)na, beta_a);
                plane_x[(size_t)k]      = (1 - t) * x_breaks[(size_t)s] + t * x_breaks[(size_t)s + 1];
                plane_rb[(size_t)k]     = (1 - t) * bore_radius[(size_t)s] + t * bore_radius[(size_t)s + 1];
                cell_segment[(size_t)k] = s;
            }
        }
        plane_x[(size_t)k]  = x_breaks[(size_t)n_segments];
        plane_rb[(size_t)k] = bore_radius[(size_t)n_segments];
    }

    // The axial map at a CONTINUOUS plane coordinate.
    //
    // The integer samples are the same expression, with the same operand order, that filled
    // plane_x and plane_rb above -- so the generator, which places nodes at integers, and the
    // warp, which places them at the fractions a lattice sits at, agree at every node they
    // share. That is the property the nozzle's header comment claims for the cross-section map,
    // extended to the axial one.
    //
    // Without it the planes could be graded while the micro nodes between them stayed uniformly
    // spaced: grading that is piecewise-linear across macro cells, with a kink at every cell
    // boundary, and the finer the lattice the more of the mesh sits on chords rather than on
    // the nozzle.
    void at(const double k, double &x, double &rb) const {
        ptrdiff_t cell = (ptrdiff_t)std::floor(k);
        if (cell < 0) cell = 0;
        if (cell > n_cells_x - 1) cell = n_cells_x - 1;
        const ptrdiff_t s = cell_segment[(size_t)cell];
        const double    a = k - (double)seg_first_cell[(size_t)s];
        const double    t = nozzle_stretch(a / (double)seg_cells[(size_t)s], beta_a);
        x  = (1 - t) * seg_x0[(size_t)s] + t * seg_x1[(size_t)s];
        rb = (1 - t) * seg_rb0[(size_t)s] + t * seg_rb1[(size_t)s];
    }
};

struct NozzleQuad {
    ptrdiff_t v[4];
    bool      outer;
};

inline void nozzle_build_quads(const NozzleSection &sec, std::vector<NozzleQuad> &quads) {
    const ptrdiff_t n       = sec.n;
    const ptrdiff_t n_perim = sec.n_perim;
    quads.clear();
    quads.reserve((size_t)(n * n + sec.n_loops * n_perim));
    for (ptrdiff_t j = 0; j < n; ++j)
        for (ptrdiff_t i = 0; i < n; ++i)
            quads.push_back({{sec.core_node(i, j), sec.core_node(i + 1, j), sec.core_node(i + 1, j + 1),
                              sec.core_node(i, j + 1)},
                             false});
    for (ptrdiff_t l = 0; l < sec.n_loops; ++l)
        for (ptrdiff_t m = 0; m < n_perim; ++m)
            quads.push_back({{sec.loop_node(l, m), sec.loop_node(l + 1, m), sec.loop_node(l + 1, m + 1),
                              sec.loop_node(l, m + 1)},
                             l >= sec.n_bore});
}

inline ptrdiff_t nozzle_first_downstream_cell(const std::vector<ptrdiff_t> &n_axial, const ptrdiff_t expansion) {
    const ptrdiff_t n_segments = (ptrdiff_t)n_axial.size();
    if (!(expansion >= 0 && expansion < n_segments)) {
        return 0;
    }
    ptrdiff_t k = 0;
    for (ptrdiff_t s = 0; s < expansion; ++s) {
        k += n_axial[(size_t)s];
    }
    return k;
}

}  // namespace smesh

#endif
