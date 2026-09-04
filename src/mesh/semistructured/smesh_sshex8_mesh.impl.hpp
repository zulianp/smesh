#ifndef SMESH_SSHEX8_MESH_IMPL_HPP
#define SMESH_SSHEX8_MESH_IMPL_HPP

#include "smesh_sshex8.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_ssquad4.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>

namespace smesh {

// Deterministic node ownership for the coordinate fills below.
//
// A lattice node on a macro-element face, edge or corner belongs to several
// macro-elements, and every one of them can compute its coordinate. The value is the
// same mathematically, but the arithmetic is not: each element interpolates from its own
// eight corners with its own weights, so the results differ in the last bits. The fills
// wrote every node from every element containing it, which means the value that survived
// was whichever store landed last, and that follows the thread schedule.
//
// The visible effect is a mesh whose coordinates change between runs. Measured on a
// 27-macro-element mesh, the x coordinate checksum varied between runs on 8 threads
// (16562.000025525689, 16562.000030174851) and was stable on one (16562.000023022294);
// y and z were unaffected only because their values happened to be exactly representable,
// so the order could not matter. Everything downstream inherits it -- fields seeded from
// the coordinates differ before any solver runs, and iteration counts vary with them.
//
// Ownership removes the race rather than ordering it: the lowest-numbered element
// containing a node writes it and the others skip it, so each node is written exactly
// once. The minimum is order independent, so building the table is itself deterministic,
// and a single-threaded run produces exactly what it produced before, since element 0 was
// already the first to write. It is also less work: a shared node is now interpolated
// once instead of once per incident element.
//
// SMESH_DETERMINISTIC_COORDS=0 restores the previous behaviour exactly -- no ownership
// table, every element writing every node it contains -- for comparing against results
// produced before this change.
inline bool sshex8_deterministic_coords() {
  static const bool on = [] {
    const char *e = std::getenv("SMESH_DETERMINISTIC_COORDS");
    return !e || std::atoi(e) != 0;
  }();
  return on;
}

template <typename idx_t>
static void sshex8_build_node_owner(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    std::vector<ptrdiff_t> &owner) {
  const int nxe = (level + 1) * (level + 1) * (level + 1);

  ptrdiff_t max_node = -1;
  for (int l = 0; l < nxe; l++) {
    for (ptrdiff_t e = 0; e < nelements; e++) {
      max_node = std::max(max_node, (ptrdiff_t)elements[l][e]);
    }
  }

  owner.assign((size_t)(max_node + 1), nelements);
  for (ptrdiff_t e = 0; e < nelements; e++) {
    for (int l = 0; l < nxe; l++) {
      const ptrdiff_t g = (ptrdiff_t)elements[l][e];
      if (e < owner[(size_t)g]) owner[(size_t)g] = e;
    }
  }
}

template <typename T>
SMESH_INLINE void
hex8_eval_f(const T x, const T y, const T z, T *const SMESH_RESTRICT f0,
            T *const SMESH_RESTRICT f1, T *const SMESH_RESTRICT f2,
            T *const SMESH_RESTRICT f3, T *const SMESH_RESTRICT f4,
            T *const SMESH_RESTRICT f5, T *const SMESH_RESTRICT f6,
            T *const SMESH_RESTRICT f7) {
  const T xm = (T(1) - x);
  const T ym = (T(1) - y);
  const T zm = (T(1) - z);

  *f0 = xm * ym * zm; // (0, 0, 0)
  *f1 = x * ym * zm;  // (1, 0, 0)
  *f2 = x * y * zm;   // (1, 1, 0)
  *f3 = xm * y * zm;  // (0, 1, 0)
  *f4 = xm * ym * z;  // (0, 0, 1)
  *f5 = x * ym * z;   // (1, 0, 1)
  *f6 = x * y * z;    // (1, 1, 1)
  *f7 = xm * y * z;   // (0, 1, 1)
}

template <typename T>
SMESH_INLINE void
proteus_hex8_eval_f(const T x, const T y, const T z, T *const SMESH_RESTRICT f0,
                    T *const SMESH_RESTRICT f1, T *const SMESH_RESTRICT f2,
                    T *const SMESH_RESTRICT f3, T *const SMESH_RESTRICT f4,
                    T *const SMESH_RESTRICT f5, T *const SMESH_RESTRICT f6,
                    T *const SMESH_RESTRICT f7) {
  const T xm = (T(1) - x);
  const T ym = (T(1) - y);
  const T zm = (T(1) - z);

  *f0 = xm * ym * zm; // (0, 0, 0)
  *f1 = x * ym * zm;  // (1, 0, 0)
  *f2 = xm * y * zm;  // (0, 1, 0)
  *f3 = x * y * zm;   // (1, 1, 0)
  *f4 = xm * ym * z;  // (0, 0, 1)
  *f5 = x * ym * z;   // (1, 0, 1)
  *f6 = xm * y * z;   // (0, 1, 1)
  *f7 = x * y * z;    // (1, 1, 1)
}

template <typename idx_t>
int sshex8_to_standard_hex8_mesh(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    idx_t *SMESH_RESTRICT *const SMESH_RESTRICT hex8_elements) {
  const int txe = sshex8_txe(level);

  int lnode[8];
  for (int zi = 0; zi < level; zi++) {
    for (int yi = 0; yi < level; yi++) {
      for (int xi = 0; xi < level; xi++) {
        lnode[0] = sshex8_lidx(level, xi, yi, zi);
        lnode[1] = sshex8_lidx(level, xi + 1, yi, zi);
        lnode[2] = sshex8_lidx(level, xi + 1, yi + 1, zi);
        lnode[3] = sshex8_lidx(level, xi, yi + 1, zi);

        lnode[4] = sshex8_lidx(level, xi, yi, zi + 1);
        lnode[5] = sshex8_lidx(level, xi + 1, yi, zi + 1);
        lnode[6] = sshex8_lidx(level, xi + 1, yi + 1, zi + 1);
        lnode[7] = sshex8_lidx(level, xi, yi + 1, zi + 1);

        int le = zi * level * level + yi * level + xi;
        SMESH_ASSERT(le < txe);

        for (int l = 0; l < 8; l++) {
          for (ptrdiff_t e = 0; e < nelements; e++) {
            idx_t node = elements[lnode[l]][e];
            hex8_elements[l][e * txe + le] = node;
          }
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t>
int sshex8_fill_points(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
    geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points) {
  const int proteus_to_std_hex8_corners[8] = {
      // Bottom
      sshex8_lidx(level, 0, 0, 0), sshex8_lidx(level, level, 0, 0),
      sshex8_lidx(level, level, level, 0), sshex8_lidx(level, 0, level, 0),

      // Top
      sshex8_lidx(level, 0, 0, level), sshex8_lidx(level, level, 0, level),
      sshex8_lidx(level, level, level, level),
      sshex8_lidx(level, 0, level, level)};

  // Nodes
  const geom_t h = geom_t(1) / geom_t(level);

  const bool             deterministic = sshex8_deterministic_coords();
  std::vector<ptrdiff_t> owner;
  if (deterministic) {
    sshex8_build_node_owner<idx_t>(level, nelements, elements, owner);
  }

#pragma omp parallel for collapse(4)
  for (int zi = 0; zi < level + 1; zi++) {
    for (int yi = 0; yi < level + 1; yi++) {
      for (int xi = 0; xi < level + 1; xi++) {
        for (int d = 0; d < 3; d++) {
          geom_t f[8];
          hex8_eval_f<geom_t>(geom_t(xi) * h, geom_t(yi) * h, geom_t(zi) * h, &f[0],
                      &f[1], &f[2], &f[3], &f[4], &f[5], &f[6], &f[7]);

          int lidx = sshex8_lidx(level, xi, yi, zi);

          for (ptrdiff_t e = 0; e < nelements; e++) {
            const ptrdiff_t g = (ptrdiff_t)elements[lidx][e];
            if (deterministic && owner[(size_t)g] != e) {
              // Another element owns this node and will write it.
              continue;
            }

            geom_t acc = 0;

            for (int lnode = 0; lnode < 8; lnode++) {
              const int corner_idx = proteus_to_std_hex8_corners[lnode];
              geom_t p = macro_mesh_points[d][elements[corner_idx][e]];
              acc += p * f[lnode];
            }

            points[d][g] = acc;
          }
        }
      }
    }
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename geom_t, typename ref_t>
int sshex8_fill_points_1D_map(
    const int level, const ptrdiff_t nelements,
    const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT macro_mesh_points,
    const ref_t *const SMESH_RESTRICT ref_points,
    geom_t *SMESH_RESTRICT *const SMESH_RESTRICT points) {
  const int proteus_to_std_hex8_corners[8] = {
      // Bottom
      sshex8_lidx(level, 0, 0, 0), sshex8_lidx(level, level, 0, 0),
      sshex8_lidx(level, level, level, 0), sshex8_lidx(level, 0, level, 0),

      // Top
      sshex8_lidx(level, 0, 0, level), sshex8_lidx(level, level, 0, level),
      sshex8_lidx(level, level, level, level),
      sshex8_lidx(level, 0, level, level)};

  const bool             deterministic = sshex8_deterministic_coords();
  std::vector<ptrdiff_t> owner;
  if (deterministic) {
    sshex8_build_node_owner<idx_t>(level, nelements, elements, owner);
  }

#pragma omp parallel for collapse(4)
  for (int zi = 0; zi < level + 1; zi++) {
    for (int yi = 0; yi < level + 1; yi++) {
      for (int xi = 0; xi < level + 1; xi++) {
        for (int d = 0; d < 3; d++) {
          geom_t f[8];
          hex8_eval_f<geom_t>(ref_points[xi], ref_points[yi], ref_points[zi], &f[0],
                      &f[1], &f[2], &f[3], &f[4], &f[5], &f[6], &f[7]);
          int lidx = sshex8_lidx(level, xi, yi, zi);

          for (ptrdiff_t e = 0; e < nelements; e++) {
            const ptrdiff_t g = (ptrdiff_t)elements[lidx][e];
            if (deterministic && owner[(size_t)g] != e) {
              // Another element owns this node and will write it.
              continue;
            }

            geom_t acc = 0;

            for (int lnode = 0; lnode < 8; lnode++) {
              const int corner_idx = proteus_to_std_hex8_corners[lnode];
              geom_t p = macro_mesh_points[d][elements[corner_idx][e]];
              acc += p * f[lnode];
            }

            points[d][g] = acc;
          }
        }
      }
    }
  }

  return SMESH_SUCCESS;
}
} // namespace smesh

#endif // SMESH_SSHEX8_MESH_IMPL_HPP
