#ifndef SMESH_BUILD_IMPL_HPP
#define SMESH_BUILD_IMPL_HPP

#include "smesh_base.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
void mesh_fill_hex8_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax, idx_t **SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {

  const ptrdiff_t ldz = (ny + 1) * (nx + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;
  const double hz = (zmax - zmin) * 1. / nz;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);
  SMESH_ASSERT(hz > 0);

  for (ptrdiff_t zi = 0; zi < nz; zi++) {
    for (ptrdiff_t yi = 0; yi < ny; yi++) {
      for (ptrdiff_t xi = 0; xi < nx; xi++) {
        const ptrdiff_t e = zi * ny * nx + yi * nx + xi;

        const idx_t i0 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i1 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i2 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;
        const idx_t i3 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;

        const idx_t i4 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i5 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i6 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;
        const idx_t i7 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;

        elements[0][e] = i0;
        elements[1][e] = i1;
        elements[2][e] = i2;
        elements[3][e] = i3;

        elements[4][e] = i4;
        elements[5][e] = i5;
        elements[6][e] = i6;
        elements[7][e] = i7;
      }
    }
  }

  for (ptrdiff_t zi = 0; zi <= nz; zi++) {
    for (ptrdiff_t yi = 0; yi <= ny; yi++) {
      for (ptrdiff_t xi = 0; xi <= nx; xi++) {
        ptrdiff_t node = xi * ldx + yi * ldy + zi * ldz;
        points[0][node] = (geom_t)(xmin + xi * hx);
        points[1][node] = (geom_t)(ymin + yi * hy);
        points[2][node] = (geom_t)(zmin + zi * hz);
      }
    }
  }
}

void mesh_fill_tri3_square(
    const int nx, const int ny, const geom_t xmin, const geom_t ymin,
    const geom_t xmax, const geom_t ymax,
    idx_t **SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  const ptrdiff_t nelements = 2 * nx * ny;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);

  for (ptrdiff_t yi = 0; yi < ny; yi++) {
    for (ptrdiff_t xi = 0; xi < nx; xi++) {
      const ptrdiff_t e = 2 * (yi * nx + xi);

      const idx_t i0 = (xi + 0) * ldx + (yi + 0) * ldy;
      const idx_t i1 = (xi + 1) * ldx + (yi + 0) * ldy;
      const idx_t i2 = (xi + 1) * ldx + (yi + 1) * ldy;
      const idx_t i3 = (xi + 0) * ldx + (yi + 1) * ldy;

      elements[0][e] = i0;
      elements[1][e] = i1;
      elements[2][e] = i3;

      elements[0][e + 1] = i1;
      elements[1][e + 1] = i2;
      elements[2][e + 1] = i3;
    }
  }

  for (ptrdiff_t yi = 0; yi <= ny; yi++) {
    for (ptrdiff_t xi = 0; xi <= nx; xi++) {
      ptrdiff_t node = xi * ldx + yi * ldy;
      points[0][node] = (geom_t)(xmin + xi * hx);
      points[1][node] = (geom_t)(ymin + yi * hy);
    }
  }
}

void mesh_fill_quad4_square(
    const int nx, const int ny, const geom_t xmin, const geom_t ymin,
    const geom_t xmax, const geom_t ymax,
    idx_t **SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  const ptrdiff_t nelements = nx * ny;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);

  for (ptrdiff_t yi = 0; yi < ny; yi++) {
    for (ptrdiff_t xi = 0; xi < nx; xi++) {
      const ptrdiff_t e = yi * nx + xi;

      const idx_t i0 = (xi + 0) * ldx + (yi + 0) * ldy;
      const idx_t i1 = (xi + 1) * ldx + (yi + 0) * ldy;
      const idx_t i2 = (xi + 1) * ldx + (yi + 1) * ldy;
      const idx_t i3 = (xi + 0) * ldx + (yi + 1) * ldy;

      elements[0][e] = i0;
      elements[1][e] = i1;
      elements[2][e] = i2;
      elements[3][e] = i3;
    }
  }

  for (ptrdiff_t yi = 0; yi <= ny; yi++) {
    for (ptrdiff_t xi = 0; xi <= nx; xi++) {
      ptrdiff_t node = xi * ldx + yi * ldy;
      points[0][node] = (geom_t)(xmin + xi * hx);
      points[1][node] = (geom_t)(ymin + yi * hy);
    }
  }
}

} // namespace smesh

#endif // SMESH_BUILD_IMPL_HPP