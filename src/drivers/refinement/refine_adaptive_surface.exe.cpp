#include "smesh_alloc.hpp"
#include "smesh_context.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cmath>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <type_traits>

#ifdef SMESH_HAS_SSDF
#include "ssdf.hpp"
#endif

using namespace smesh;

#ifdef SMESH_HAS_SSDF

static int parse_level(const char *text, int *level) {
  char *end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (end == text || *end != '\0' || value < 0 ||
      value > static_cast<long>(std::numeric_limits<int>::max())) {
    return SMESH_FAILURE;
  }
  *level = static_cast<int>(value);
  return SMESH_SUCCESS;
}

static int parse_distance(const char *text, geom_t *distance) {
  char *end = nullptr;
  const double value = std::strtod(text, &end);
  if (end == text || *end != '\0' || !(value >= 0.0)) {
    return SMESH_FAILURE;
  }
  *distance = static_cast<geom_t>(value);
  return SMESH_SUCCESS;
}

static int query_surface_distance(const ptrdiff_t n_query, const geom_t *qx,
                                  const geom_t *qy, const geom_t *qz,
                                  const ptrdiff_t n_tri, const idx_t *s0,
                                  const idx_t *s1, const idx_t *s2,
                                  const ptrdiff_t n_surf_nodes, const geom_t *sx,
                                  const geom_t *sy, const geom_t *sz,
                                  geom_t *dist) {
  if constexpr (std::is_same<idx_t, int>::value) {
    return ssdf::edf_select<geom_t, geom_t, idx_t>(
        n_query, qx, qy, qz, n_tri, s0, s1, s2, n_surf_nodes, sx, sy, sz, dist);
  } else {
    return ssdf::edf_celllist<geom_t, geom_t, idx_t>(
        n_query, qx, qy, qz, n_tri, s0, s1, s2, n_surf_nodes, sx, sy, sz, dist);
  }
}

static int mark_near_surface(const Mesh &mesh, const Mesh &surface,
                             const geom_t distance, element_idx_t **marked,
                             ptrdiff_t *n_marked) {
  *marked = nullptr;
  *n_marked = 0;

  if (surface.spatial_dimension() != 3 || surface.n_nodes() <= 0) {
    fprintf(stderr, "refine_adaptive_surface: surface must be a 3D triangle mesh\n");
    return SMESH_FAILURE;
  }

  ptrdiff_t n_tri = 0;
  for (size_t b = 0; b < surface.n_blocks(); ++b) {
    const enum ElemType type = surface.element_type(static_cast<block_idx_t>(b));
    if (type != TRI3 && type != TRISHELL3) {
      fprintf(stderr,
              "refine_adaptive_surface: surface block %zu is %s, expected TRI3 or TRISHELL3\n",
              b, type_to_string(type));
      return SMESH_FAILURE;
    }
    n_tri += surface.n_elements(static_cast<block_idx_t>(b));
  }
  if (n_tri <= 0) {
    fprintf(stderr, "refine_adaptive_surface: surface has no triangles\n");
    return SMESH_FAILURE;
  }

  const ptrdiff_t n_nodes = mesh.n_nodes();
  ptrdiff_t n_elem = 0;
  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    n_elem += mesh.n_elements(static_cast<block_idx_t>(b));
  }
  if (n_nodes <= 0 || n_elem <= 0) {
    return SMESH_SUCCESS;
  }
  if (n_elem > static_cast<ptrdiff_t>(std::numeric_limits<element_idx_t>::max())) {
    fprintf(stderr, "refine_adaptive_surface: element count exceeds element_idx_t\n");
    return SMESH_FAILURE;
  }

  const ptrdiff_t n_query = n_nodes + n_elem;
  geom_t *qx = static_cast<geom_t *>(SMESH_ALLOC(static_cast<size_t>(n_query) * sizeof(geom_t)));
  geom_t *qy = static_cast<geom_t *>(SMESH_ALLOC(static_cast<size_t>(n_query) * sizeof(geom_t)));
  geom_t *qz = static_cast<geom_t *>(SMESH_ALLOC(static_cast<size_t>(n_query) * sizeof(geom_t)));
  geom_t *dist = static_cast<geom_t *>(SMESH_ALLOC(static_cast<size_t>(n_query) * sizeof(geom_t)));
  idx_t *s0 = static_cast<idx_t *>(SMESH_ALLOC(static_cast<size_t>(n_tri) * sizeof(idx_t)));
  idx_t *s1 = static_cast<idx_t *>(SMESH_ALLOC(static_cast<size_t>(n_tri) * sizeof(idx_t)));
  idx_t *s2 = static_cast<idx_t *>(SMESH_ALLOC(static_cast<size_t>(n_tri) * sizeof(idx_t)));
  if (!qx || !qy || !qz || !dist || !s0 || !s1 || !s2) {
    fprintf(stderr, "refine_adaptive_surface: allocation failed\n");
    SMESH_FREE(qx);
    SMESH_FREE(qy);
    SMESH_FREE(qz);
    SMESH_FREE(dist);
    SMESH_FREE(s0);
    SMESH_FREE(s1);
    SMESH_FREE(s2);
    return SMESH_FAILURE;
  }

  ptrdiff_t tri_cursor = 0;
  for (size_t b = 0; b < surface.n_blocks(); ++b) {
    const auto elems = surface.elements(static_cast<block_idx_t>(b));
    const idx_t *const b0 = elems->data()[0];
    const idx_t *const b1 = elems->data()[1];
    const idx_t *const b2 = elems->data()[2];
    const ptrdiff_t ne = surface.n_elements(static_cast<block_idx_t>(b));
    for (ptrdiff_t e = 0; e < ne; ++e) {
      s0[tri_cursor] = b0[e];
      s1[tri_cursor] = b1[e];
      s2[tri_cursor] = b2[e];
      ++tri_cursor;
    }
  }

  const auto pts = mesh.points();
  const geom_t *const px = pts->data()[0];
  const geom_t *const py = pts->data()[1];
  const geom_t *const pz = pts->data()[2];
  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    qx[i] = px[i];
    qy[i] = py[i];
    qz[i] = pz[i];
  }

  ptrdiff_t elem_base = 0;
  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    const auto vol = mesh.elements(static_cast<block_idx_t>(b));
    const idx_t *const e0 = vol->data()[0];
    const idx_t *const e1 = vol->data()[1];
    const idx_t *const e2 = vol->data()[2];
    const idx_t *const e3 = vol->data()[3];
    const ptrdiff_t ne = mesh.n_elements(static_cast<block_idx_t>(b));
    for (ptrdiff_t e = 0; e < ne; ++e) {
      const idx_t n0 = e0[e];
      const idx_t n1 = e1[e];
      const idx_t n2 = e2[e];
      const idx_t n3 = e3[e];
      const ptrdiff_t q = n_nodes + elem_base + e;
      qx[q] = (px[n0] + px[n1] + px[n2] + px[n3]) * geom_t(0.25);
      qy[q] = (py[n0] + py[n1] + py[n2] + py[n3]) * geom_t(0.25);
      qz[q] = (pz[n0] + pz[n1] + pz[n2] + pz[n3]) * geom_t(0.25);
    }
    elem_base += ne;
  }

  const geom_t huge =
      std::sqrt(std::numeric_limits<geom_t>::max()) * geom_t(0.25);
  for (ptrdiff_t i = 0; i < n_query; ++i) {
    dist[i] = huge;
  }

  const auto surf_pts = surface.points();
  const int query_rc = query_surface_distance(
      n_query, qx, qy, qz, n_tri, s0, s1, s2, surface.n_nodes(),
      surf_pts->data()[0], surf_pts->data()[1], surf_pts->data()[2], dist);

  SMESH_FREE(qx);
  SMESH_FREE(qy);
  SMESH_FREE(qz);
  SMESH_FREE(s0);
  SMESH_FREE(s1);
  SMESH_FREE(s2);

  if (query_rc != 0) {
    fprintf(stderr, "refine_adaptive_surface: surface distance query failed\n");
    SMESH_FREE(dist);
    return SMESH_FAILURE;
  }

  ptrdiff_t count = 0;
  elem_base = 0;
  for (size_t b = 0; b < mesh.n_blocks(); ++b) {
    const auto vol = mesh.elements(static_cast<block_idx_t>(b));
    const idx_t *const e0 = vol->data()[0];
    const idx_t *const e1 = vol->data()[1];
    const idx_t *const e2 = vol->data()[2];
    const idx_t *const e3 = vol->data()[3];
    const ptrdiff_t ne = mesh.n_elements(static_cast<block_idx_t>(b));
    for (ptrdiff_t e = 0; e < ne; ++e) {
      geom_t d = dist[n_nodes + elem_base + e];
      const geom_t d0 = dist[e0[e]];
      const geom_t d1 = dist[e1[e]];
      const geom_t d2 = dist[e2[e]];
      const geom_t d3 = dist[e3[e]];
      d = d < d0 ? d : d0;
      d = d < d1 ? d : d1;
      d = d < d2 ? d : d2;
      d = d < d3 ? d : d3;
      count += (d <= distance);
    }
    elem_base += ne;
  }

  element_idx_t *marks = nullptr;
  if (count > 0) {
    marks = static_cast<element_idx_t *>(
        SMESH_ALLOC(static_cast<size_t>(count) * sizeof(element_idx_t)));
    if (!marks) {
      fprintf(stderr, "refine_adaptive_surface: allocation failed\n");
      SMESH_FREE(dist);
      return SMESH_FAILURE;
    }
    ptrdiff_t cursor = 0;
    elem_base = 0;
    for (size_t b = 0; b < mesh.n_blocks(); ++b) {
      const auto vol = mesh.elements(static_cast<block_idx_t>(b));
      const idx_t *const e0 = vol->data()[0];
      const idx_t *const e1 = vol->data()[1];
      const idx_t *const e2 = vol->data()[2];
      const idx_t *const e3 = vol->data()[3];
      const ptrdiff_t ne = mesh.n_elements(static_cast<block_idx_t>(b));
      for (ptrdiff_t e = 0; e < ne; ++e) {
        geom_t d = dist[n_nodes + elem_base + e];
        const geom_t d0 = dist[e0[e]];
        const geom_t d1 = dist[e1[e]];
        const geom_t d2 = dist[e2[e]];
        const geom_t d3 = dist[e3[e]];
        d = d < d0 ? d : d0;
        d = d < d1 ? d : d1;
        d = d < d2 ? d : d2;
        d = d < d3 ? d : d3;
        if (d <= distance) {
          marks[cursor++] = static_cast<element_idx_t>(elem_base + e);
        }
      }
      elem_base += ne;
    }
  }

  SMESH_FREE(dist);
  *marked = marks;
  *n_marked = count;
  return SMESH_SUCCESS;
}

#endif

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("refine_adaptive_surface.exe");
  auto ctx = initialize(argc, argv);
  auto comm = ctx->communicator();

  if (argc != 6) {
    if (!comm->rank()) {
      fprintf(stderr,
              "Usage: %s <mesh_folder> <level> <surface_folder> <distance> <output_folder>\n"
              "Marks TET4 elements whose nodes or centroid lie within <distance> of the\n"
              "triangle surface and red-refines them to <level>. Each such element\n"
              "raises a face neighbor to at least one level less, and that repeats.\n",
              argv[0]);
    }
    return SMESH_FAILURE;
  }

#ifndef SMESH_HAS_SSDF
  if (!comm->rank()) {
    fprintf(stderr,
            "refine_adaptive_surface: built without ssdf. Configure with -DSMESH_ENABLE_SSDF=ON.\n");
  }
  return SMESH_FAILURE;
#else
  int level = 0;
  geom_t distance = 0;
  if (parse_level(argv[2], &level) != SMESH_SUCCESS ||
      parse_distance(argv[4], &distance) != SMESH_SUCCESS) {
    if (!comm->rank()) {
      fprintf(stderr,
              "refine_adaptive_surface: level must be an integer >= 0 and distance a number >= 0\n");
    }
    return SMESH_FAILURE;
  }

  auto mesh = Mesh::create_from_file(comm, Path(argv[1]));
  if (!mesh) {
    return SMESH_FAILURE;
  }
  if (mesh->n_blocks() < 1 || mesh->spatial_dimension() != 3) {
    if (!comm->rank()) {
      fprintf(stderr, "refine_adaptive_surface: mesh must be TET4 blocks in 3D\n");
    }
    return SMESH_FAILURE;
  }
  for (size_t b = 0; b < mesh->n_blocks(); ++b) {
    if (mesh->element_type(static_cast<block_idx_t>(b)) != TET4) {
      if (!comm->rank()) {
        fprintf(stderr, "refine_adaptive_surface: mesh must be TET4 blocks in 3D\n");
      }
      return SMESH_FAILURE;
    }
  }

  auto surface = Mesh::create_from_file(comm, Path(argv[3]));
  if (!surface) {
    return SMESH_FAILURE;
  }

  element_idx_t *marked = nullptr;
  ptrdiff_t n_marked = 0;
  if (mark_near_surface(*mesh, *surface, distance, &marked, &n_marked) !=
      SMESH_SUCCESS) {
    return SMESH_FAILURE;
  }

  if (!comm->rank()) {
    fprintf(stdout,
            "refine_adaptive_surface: marked %td / %td elements (distance <= %g)\n",
            n_marked, mesh->n_elements(), static_cast<double>(distance));
  }

  auto refined = refine_adaptive(mesh, marked, n_marked, level);
  SMESH_FREE(marked);
  if (!refined) {
    return SMESH_FAILURE;
  }
  return refined->write(Path(argv[5]));
#endif
}
