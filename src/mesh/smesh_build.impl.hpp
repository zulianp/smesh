#ifndef SMESH_BUILD_IMPL_HPP
#define SMESH_BUILD_IMPL_HPP

#include "smesh_base.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
void mesh_fill_hex8_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
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

  // #pragma omp parallel for collapse(3)
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

  // #pragma omp parallel for collapse(3)
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

template <typename idx_t, typename geom_t>
void mesh_fill_tri3_square(
    const int nx, const int ny, const geom_t xmin, const geom_t ymin,
    const geom_t xmax, const geom_t ymax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  //   const ptrdiff_t nelements = 2 * nx * ny;
  //   const ptrdiff_t nnodes = (nx + 1) * (ny + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);

  // #pragma omp parallel for collapse(2)
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

  // #pragma omp parallel for collapse(2)
  for (ptrdiff_t yi = 0; yi <= ny; yi++) {
    for (ptrdiff_t xi = 0; xi <= nx; xi++) {
      ptrdiff_t node = xi * ldx + yi * ldy;
      points[0][node] = (geom_t)(xmin + xi * hx);
      points[1][node] = (geom_t)(ymin + yi * hy);
    }
  }
}

template <typename idx_t, typename geom_t>
void mesh_fill_quad4_square(
    const int nx, const int ny, const geom_t xmin, const geom_t ymin,
    const geom_t xmax, const geom_t ymax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  //   const ptrdiff_t nelements = nx * ny;
  //   const ptrdiff_t nnodes = (nx + 1) * (ny + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);

  // #pragma omp parallel for collapse(2)
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

  // #pragma omp parallel for collapse(2)
  for (ptrdiff_t yi = 0; yi <= ny; yi++) {
    for (ptrdiff_t xi = 0; xi <= nx; xi++) {
      ptrdiff_t node = xi * ldx + yi * ldy;
      points[0][node] = (geom_t)(xmin + xi * hx);
      points[1][node] = (geom_t)(ymin + yi * hy);
    }
  }
}

template <typename idx_t, typename geom_t>
void mesh_fill_tet4_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  // const ptrdiff_t nelements = ptrdiff_t(nx) * ptrdiff_t(ny) * ptrdiff_t(nz);
  // const ptrdiff_t nnodes_total = nnodes_vertices + nelements;
  const ptrdiff_t nnodes_vertices =
      ptrdiff_t(nx + 1) * ptrdiff_t(ny + 1) * ptrdiff_t(nz + 1);
  const ptrdiff_t ldz = (ny + 1) * (nx + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;
  const double hz = (zmax - zmin) * 1. / nz;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);
  SMESH_ASSERT(hz > 0);

  static const int face_nodes[6][4] = {{0, 1, 2, 3}, {4, 7, 6, 5},
                                       {0, 4, 5, 1}, {3, 2, 6, 7},
                                       {0, 3, 7, 4}, {1, 5, 6, 2}};

  for (ptrdiff_t zi = 0; zi < nz; zi++) {
    for (ptrdiff_t yi = 0; yi < ny; yi++) {
      for (ptrdiff_t xi = 0; xi < nx; xi++) {
        const idx_t i0 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i1 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i2 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;
        const idx_t i3 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;

        const idx_t i4 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i5 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i6 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;
        const idx_t i7 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;

        const idx_t cube_nodes[8] = {i0, i1, i2, i3, i4, i5, i6, i7};
        const ptrdiff_t elem_index = zi * ny * nx + yi * nx + xi;
        const idx_t center_idx = nnodes_vertices + elem_index;
        const ptrdiff_t base = elem_index * 12;

        for (int face = 0; face < 6; face++) {
          const int *fn = face_nodes[face];

          const ptrdiff_t t0 = base + face * 2 + 0;
          elements[0][t0] = cube_nodes[fn[0]];
          elements[1][t0] = cube_nodes[fn[1]];
          elements[2][t0] = cube_nodes[fn[2]];
          elements[3][t0] = center_idx;

          const ptrdiff_t t1 = base + face * 2 + 1;
          elements[0][t1] = cube_nodes[fn[0]];
          elements[1][t1] = cube_nodes[fn[2]];
          elements[2][t1] = cube_nodes[fn[3]];
          elements[3][t1] = center_idx;
        }
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

  for (ptrdiff_t zi = 0; zi < nz; zi++) {
    for (ptrdiff_t yi = 0; yi < ny; yi++) {
      for (ptrdiff_t xi = 0; xi < nx; xi++) {
        const ptrdiff_t elem_index = zi * ny * nx + yi * nx + xi;
        const idx_t center_idx = nnodes_vertices + elem_index;

        points[0][center_idx] = (geom_t)(xmin + (xi + 0.5) * hx);
        points[1][center_idx] = (geom_t)(ymin + (yi + 0.5) * hy);
        points[2][center_idx] = (geom_t)(zmin + (zi + 0.5) * hz);
      }
    }
  }
}

template <typename idx_t, typename geom_t>
void mesh_fill_hex8_bidomain_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax, const int dim_split, const idx_t split_index,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT left_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT right_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  const ptrdiff_t nelements = ptrdiff_t(nx) * ptrdiff_t(ny) * ptrdiff_t(nz);

  SMESH_ASSERT(dim_split >= 0);
  SMESH_ASSERT(dim_split <= 2);
  const ptrdiff_t split = ptrdiff_t(split_index);
  SMESH_ASSERT(split >= 0);
  const ptrdiff_t split_dim = (dim_split == 0   ? ptrdiff_t(nx)
                               : dim_split == 1 ? ptrdiff_t(ny)
                                                : ptrdiff_t(nz));
  SMESH_ASSERT(split <= split_dim);
  SMESH_UNUSED(split_dim);
  SMESH_UNUSED(nelements);

  const ptrdiff_t ldz = (ny + 1) * (nx + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const geom_t hx = (xmax - xmin) * 1. / nx;
  const geom_t hy = (ymax - ymin) * 1. / ny;
  const geom_t hz = (zmax - zmin) * 1. / nz;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);
  SMESH_ASSERT(hz > 0);

  ptrdiff_t left_elements_count = 0;
  ptrdiff_t right_elements_count = 0;
  for (ptrdiff_t zi = 0; zi < nz; zi++) {
    for (ptrdiff_t yi = 0; yi < ny; yi++) {
      for (ptrdiff_t xi = 0; xi < nx; xi++) {
        const idx_t i0 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i1 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i2 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;
        const idx_t i3 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;

        const idx_t i4 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i5 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i6 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;
        const idx_t i7 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;

        ptrdiff_t idx3[3] = {xi, yi, zi};

        if (idx3[dim_split] < split_index) {
          left_elements[0][left_elements_count] = i0;
          left_elements[1][left_elements_count] = i1;
          left_elements[2][left_elements_count] = i2;
          left_elements[3][left_elements_count] = i3;

          left_elements[4][left_elements_count] = i4;
          left_elements[5][left_elements_count] = i5;
          left_elements[6][left_elements_count] = i6;
          left_elements[7][left_elements_count] = i7;

          left_elements_count++;
        } else {
          right_elements[0][right_elements_count] = i0;
          right_elements[1][right_elements_count] = i1;
          right_elements[2][right_elements_count] = i2;
          right_elements[3][right_elements_count] = i3;

          right_elements[4][right_elements_count] = i4;
          right_elements[5][right_elements_count] = i5;
          right_elements[6][right_elements_count] = i6;
          right_elements[7][right_elements_count] = i7;

          right_elements_count++;
        }
      }
    }
  }

  const ptrdiff_t expected_left = [&]() -> ptrdiff_t {
    switch (dim_split) {
    case 0:
      return split * ptrdiff_t(ny) * ptrdiff_t(nz);
    case 1:
      return ptrdiff_t(nx) * split * ptrdiff_t(nz);
    default:
      return ptrdiff_t(nx) * ptrdiff_t(ny) * split;
    }
  }();
  const ptrdiff_t expected_right = nelements - expected_left;
  SMESH_UNUSED(expected_left);
  SMESH_UNUSED(expected_right);

  SMESH_ASSERT(left_elements_count == expected_left);
  SMESH_ASSERT(right_elements_count == expected_right);

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

template <typename idx_t, typename geom_t>
void mesh_fill_hex8_checkerboard_cube(
    const int nx, const int ny, const int nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT white_elements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT black_elements,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {

  const ptrdiff_t nelements = ptrdiff_t(nx) * ptrdiff_t(ny) * ptrdiff_t(nz);

  if (nx % 2 != 0 || ny % 2 != 0 || nz % 2 != 0) {
    SMESH_ERROR("nx, ny, and nz must be even");
  }

  const ptrdiff_t ldz = (ny + 1) * (nx + 1);
  const ptrdiff_t ldy = nx + 1;
  const ptrdiff_t ldx = 1;

  const geom_t hx = (xmax - xmin) * geom_t(1) / geom_t(nx);
  const geom_t hy = (ymax - ymin) * geom_t(1) / geom_t(ny);
  const geom_t hz = (zmax - zmin) * geom_t(1) / geom_t(nz);

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);
  SMESH_ASSERT(hz > 0);

  ptrdiff_t white_elements_count = 0;
  ptrdiff_t black_elements_count = 0;

  for (ptrdiff_t zi = 0; zi < nz; zi++) {
    for (ptrdiff_t yi = 0; yi < ny; yi++) {
      for (ptrdiff_t xi = 0; xi < nx; xi++) {

        const idx_t i0 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i1 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 0) * ldz;
        const idx_t i2 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;
        const idx_t i3 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 0) * ldz;

        const idx_t i4 = (xi + 0) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i5 = (xi + 1) * ldx + (yi + 0) * ldy + (zi + 1) * ldz;
        const idx_t i6 = (xi + 1) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;
        const idx_t i7 = (xi + 0) * ldx + (yi + 1) * ldy + (zi + 1) * ldz;

        if ((xi + yi + zi) % 2 == 0) {
          white_elements[0][white_elements_count] = i0;
          white_elements[1][white_elements_count] = i1;
          white_elements[2][white_elements_count] = i2;
          white_elements[3][white_elements_count] = i3;

          white_elements[4][white_elements_count] = i4;
          white_elements[5][white_elements_count] = i5;
          white_elements[6][white_elements_count] = i6;
          white_elements[7][white_elements_count] = i7;

          white_elements_count++;
        } else {
          black_elements[0][black_elements_count] = i0;
          black_elements[1][black_elements_count] = i1;
          black_elements[2][black_elements_count] = i2;
          black_elements[3][black_elements_count] = i3;

          black_elements[4][black_elements_count] = i4;
          black_elements[5][black_elements_count] = i5;
          black_elements[6][black_elements_count] = i6;
          black_elements[7][black_elements_count] = i7;

          black_elements_count++;
        }
      }
    }
  }

  SMESH_ASSERT(white_elements_count == nelements / 2);
  SMESH_ASSERT(black_elements_count == nelements / 2);
  SMESH_UNUSED(nelements);

  for (ptrdiff_t zi = 0; zi <= nz; zi++) {
    for (ptrdiff_t yi = 0; yi <= ny; yi++) {
      for (ptrdiff_t xi = 0; xi <= nx; xi++) {
        ptrdiff_t node = xi * ldx + yi * ldy + zi * ldz;
        points[0][node] = xmin + geom_t(xi) * hx;
        points[1][node] = ymin + geom_t(yi) * hy;
        points[2][node] = zmin + geom_t(zi) * hz;
      }
    }
  }
}

template <typename idx_t, typename geom_t>
void mesh_fill_quad4_ring(const geom_t inner_radius, const geom_t outer_radius,
                          const ptrdiff_t nlayers, const ptrdiff_t nelements,
                          idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
                          geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points) {
  // ptrdiff_t nnodes = nelements * 2;
  const geom_t dangle = 2 * geom_t(M_PI) / geom_t(nelements);
  const geom_t dh = (outer_radius - inner_radius) / nlayers;

  for (ptrdiff_t l = 0; l <= nlayers; l++) {
    for (ptrdiff_t i = 0; i < nelements; i++) {
      ptrdiff_t idx = l * nelements + i;
      points[0][idx] = cos(dangle * i) * (inner_radius + dh * l);
      points[1][idx] = sin(dangle * i) * (inner_radius + dh * l);
    }
  }

  for (ptrdiff_t l = 0; l < nlayers; l++) {
    for (ptrdiff_t i = 0; i < nelements; i++) {
      ptrdiff_t idx = l * nelements + i;
      elements[0][idx] = l * nelements + (i + 1) % nelements;
      elements[1][idx] = l * nelements + i;
      elements[2][idx] = elements[1][idx] + nelements;
      elements[3][idx] = elements[0][idx] + nelements;
    }
  }
}

} // namespace smesh

#endif // SMESH_BUILD_IMPL_HPP