#ifndef SMESH_GENCUBE_IMPL_HPP
#define SMESH_GENCUBE_IMPL_HPP

#include "smesh_gencube.hpp"

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_glob.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"

namespace smesh {

template <typename idx_t, typename geom_t>
int mesh_hex8_cube_to_folder(const Path &folder, const ptrdiff_t nx,
                             const ptrdiff_t ny, const ptrdiff_t nz,
                             const geom_t xmin, const geom_t ymin,
                             const geom_t zmin, const geom_t xmax,
                             const geom_t ymax, const geom_t zmax,
                             const ptrdiff_t z_chunk_size) {
  SMESH_TRACE_SCOPE("mesh_hex8_cube_to_folder");

  const ptrdiff_t ldz = (ny + 1) * (nx + 1);
  const ptrdiff_t ldy = (nx + 1);
  const ptrdiff_t ldx = 1;

  SMESH_ASSERT(nx > 0 && ny > 0 && nz > 0);
  const double hx = (xmax - xmin) * 1. / nx;
  const double hy = (ymax - ymin) * 1. / ny;
  const double hz = (zmax - zmin) * 1. / nz;

  SMESH_ASSERT(hx > 0);
  SMESH_ASSERT(hy > 0);
  SMESH_ASSERT(hz > 0);

  printf("writing to %s\n", folder.c_str());

  create_directory(folder);

  mesh_write_yaml_basic(folder, HEX8, nx * ny * nz, 3,
                        (nx + 1) * (ny + 1) * (nz + 1));

  const ptrdiff_t ix[8] = {0, 1, 1, 0, 0, 1, 1, 0};
  const ptrdiff_t iy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
  const ptrdiff_t iz[8] = {0, 0, 0, 0, 1, 1, 1, 1};

  idx_t *elements_chunk =
      (idx_t *)malloc(z_chunk_size * ny * nx * sizeof(idx_t));

  for (int v = 0; v < 8; v++) {
    Path path = folder / Path("i" + std::to_string(v) + "." +
                              std::string(TypeToString<idx_t>::value()));
    FILE *f = fopen(path.c_str(), "wb");
    if (!f) {
      free(elements_chunk);
      SMESH_ERROR("mesh_hex8_cube_to_folder: Failed to open file %s\n",
                  path.c_str());
      return SMESH_FAILURE;
    }

    for (ptrdiff_t zi_chunk = 0; zi_chunk < nz; zi_chunk += z_chunk_size) {
      const ptrdiff_t zi_end = std::min(zi_chunk + z_chunk_size, nz);

#pragma omp parallel for collapse(2)
      for (ptrdiff_t zi = zi_chunk; zi < zi_end; zi++) {
        for (ptrdiff_t yi = 0; yi < ny; yi++) {
          for (ptrdiff_t xi = 0; xi < nx; xi++) {
            const ptrdiff_t e = (zi - zi_chunk) * (ny * nx) + yi * (nx) + xi;

            const idx_t iv =
                (xi + ix[v]) * ldx + (yi + iy[v]) * ldy + (zi + iz[v]) * ldz;
            elements_chunk[e] = iv;
          }
        }
      }
      const ptrdiff_t chunk_len = zi_end - zi_chunk;
      fwrite(elements_chunk, sizeof(idx_t), chunk_len * ny * nx, f);
    }

    fclose(f);
  }

  free(elements_chunk);

  geom_t *points_chunk =
      (geom_t *)malloc(z_chunk_size * (ny + 1) * (nx + 1) * sizeof(geom_t));

  auto write_points_chunk = [&](const std::string &coord, const geom_t shift,
                                const geom_t hx, const geom_t hy,
                                const geom_t hz) -> int {
    Path path =
        folder / Path(coord + "." + std::string(TypeToString<geom_t>::value()));
    FILE *f = fopen(path.c_str(), "wb");
    if (!f) {
      SMESH_ERROR("mesh_hex8_cube_to_folder: Failed to open file %s\n",
                  path.c_str());
      return SMESH_FAILURE;
    }
    for (ptrdiff_t zi_chunk = 0; zi_chunk < nz + 1; zi_chunk += z_chunk_size) {
      const ptrdiff_t zi_end = std::min(zi_chunk + z_chunk_size, nz + 1);

#pragma omp parallel for collapse(2)
      for (ptrdiff_t zi = zi_chunk; zi < zi_end; zi++) {
        for (ptrdiff_t yi = 0; yi <= ny; yi++) {
          for (ptrdiff_t xi = 0; xi <= nx; xi++) {
            const ptrdiff_t node = (zi - zi_chunk) * ldz + yi * ldy + xi * ldx;
            points_chunk[node] = (geom_t)(shift + xi * hx + yi * hy + zi * hz);
          }
        }
      }

      const ptrdiff_t chunk_len = zi_end - zi_chunk;
      fwrite(points_chunk, sizeof(geom_t), chunk_len * (ny + 1) * (nx + 1), f);
    }

    fclose(f);
    return SMESH_SUCCESS;
  };

  int ret = write_points_chunk("x", xmin, hx, 0, 0) != SMESH_SUCCESS ||
            write_points_chunk("y", ymin, 0, hy, 0) != SMESH_SUCCESS ||
            write_points_chunk("z", zmin, 0, 0, hz) != SMESH_SUCCESS;

  free(points_chunk);

  return ret ? SMESH_FAILURE : SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_GENCUBE_HPP
