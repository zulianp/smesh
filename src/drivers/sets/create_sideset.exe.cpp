#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mask.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sidesets.impl.hpp"
#include "smesh_tracer.hpp"
#include <cmath>
#include <cstdlib>
#include <stdio.h>

using namespace smesh;

static inline int side_num_corner_nodes(const enum ElemType side_type) {
  // We only use corner nodes for geometry/angles.
  switch (side_type) {
  case NODE1:
    return 1;
  case EDGE2:
  case EDGE3:
  case BEAM2:
    return 2;
  case TRI3:
  case TRI6:
  case TRI10:
    return 3;
  case QUAD4:
  case QUADSHELL4:
    return 4;
  default:
    return -1;
  }
}

static inline bool compute_side_normal(const int dim,
                                      const enum ElemType element_type,
                                      const idx_t *const SMESH_RESTRICT
                                          *const SMESH_RESTRICT elements,
                                      const geom_t *const SMESH_RESTRICT
                                          *const SMESH_RESTRICT points,
                                      const element_idx_t parent,
                                      const i16 lfi, real_t *out_unit_normal) {
  LocalSideTable lst;
  lst.fill(element_type);

  const enum ElemType st = side_type(element_type);
  const int ncorner = side_num_corner_nodes(st);
  if (ncorner < 0) {
    return false;
  }

  if (dim == 2) {
    if (ncorner != 2) {
      return false;
    }

    const idx_t n0 = elements[lst(lfi, 0)][parent];
    const idx_t n1 = elements[lst(lfi, 1)][parent];

    const real_t x0 = points[0][n0];
    const real_t y0 = points[1][n0];
    const real_t x1 = points[0][n1];
    const real_t y1 = points[1][n1];

    // Perpendicular to edge (p1 - p0), normalized.
    const real_t nx = -(y1 - y0);
    const real_t ny = (x1 - x0);
    const real_t len = std::sqrt(nx * nx + ny * ny);
    if (len <= 0) {
      return false;
    }

    out_unit_normal[0] = nx / len;
    out_unit_normal[1] = ny / len;
    out_unit_normal[2] = 0;
    return true;
  }

  if (dim == 3) {
    if (ncorner < 3) {
      return false;
    }

    const idx_t n0 = elements[lst(lfi, 0)][parent];
    const idx_t n1 = elements[lst(lfi, 1)][parent];
    const idx_t n2 = elements[lst(lfi, 2)][parent];

    const real_t p0[3] = {points[0][n0], points[1][n0], points[2][n0]};
    const real_t p1[3] = {points[0][n1], points[1][n1], points[2][n1]};
    const real_t p2[3] = {points[0][n2], points[1][n2], points[2][n2]};

    const real_t u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    const real_t v[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};

    const real_t nx = (u[1] * v[2]) - (u[2] * v[1]);
    const real_t ny = (u[2] * v[0]) - (u[0] * v[2]);
    const real_t nz = (u[0] * v[1]) - (u[1] * v[0]);
    const real_t len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len <= 0) {
      return false;
    }

    out_unit_normal[0] = nx / len;
    out_unit_normal[1] = ny / len;
    out_unit_normal[2] = nz / len;
    return true;
  }

  return false;
}

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("create_sideset.exe");
  auto ctx = smesh::initialize_serial(argc, argv);

  if (argc != 7) {
    fprintf(stderr,
            "Usage: %s <mesh_folder> <x> <y> <z> <angle_threshold> "
            "<output_folder>\n",
            argv[0]);
    return SMESH_FAILURE;
  }

  int ret = SMESH_SUCCESS;
  {
    const geom_t roi[3] = {(geom_t)std::atof(argv[2]), (geom_t)std::atof(argv[3]),
                           (geom_t)std::atof(argv[4])};
    // NOTE: We follow the reference: compare fabs(dot(n0,n1)) > threshold,
    // so the threshold is a cosine similarity bound in [0, 1].
    const real_t angle_threshold = (real_t)std::atof(argv[5]);
    const Path output_folder(argv[6]);

    auto mesh = Mesh::create_from_file(ctx->communicator(), Path(argv[1]));
    auto sideset = skin_sideset(mesh);

    const ptrdiff_t n_surf = sideset->parent()->size();
    const int dim = mesh->spatial_dimension();
    auto points = mesh->points()->data();
    auto elements = mesh->elements()->data();
    const enum ElemType element_type = mesh->element_type();

    // TODO: Use node_to_element_graph instead of creating a new one
    // Build node-to-element CSR once (reference behavior).

    auto n2e = mesh->node_to_element_graph();
  
    // -------------------------------------------------------------------------
    // Find approximately closest surface side (seed).
    // -------------------------------------------------------------------------
    const auto surf_parent = sideset->parent()->data();
    const auto surf_lfi = sideset->lfi()->data();

    LocalSideTable lst;
    lst.fill(element_type);

    const enum ElemType st = side_type(element_type);
    const int ncorner = side_num_corner_nodes(st);
    if (ncorner < 0) {
      SMESH_ERROR("Unsupported element/side type for create_sideset\n");
      return SMESH_FAILURE;
    }

    element_idx_t seed_side = invalid_idx<element_idx_t>();
    real_t best_sq_dist = (real_t)1e300;

    for (ptrdiff_t si = 0; si < n_surf; ++si) {
      const element_idx_t e = surf_parent[si];
      const i16 s = surf_lfi[si];

      real_t bary[3] = {0, 0, 0};
      for (int ln = 0; ln < ncorner; ++ln) {
        const idx_t node = elements[lst(s, ln)][e];
        bary[0] += points[0][node];
        bary[1] += points[1][node];
        if (dim == 3) {
          bary[2] += points[2][node];
        }
      }

      const real_t inv = (real_t)1 / (real_t)ncorner;
      bary[0] *= inv;
      bary[1] *= inv;
      bary[2] *= inv;

      const real_t dx = bary[0] - (real_t)roi[0];
      const real_t dy = bary[1] - (real_t)roi[1];
      const real_t dz = (dim == 3) ? (bary[2] - (real_t)roi[2]) : (real_t)0;
      const real_t sq_dist = dx * dx + dy * dy + dz * dz;

      if (sq_dist < best_sq_dist) {
        best_sq_dist = sq_dist;
        seed_side = (element_idx_t)si;
      }
    }

    if (seed_side == invalid_idx<element_idx_t>()) {
      SMESH_ERROR("Unable to find seed side\n");
      return SMESH_FAILURE;
    }

    auto selected =
        create_host_buffer<mask_t>(mask_count(sideset->parent()->size()));

    const int err = sideset_select_propagate(
        sideset->parent()->size(), sideset->parent()->data(),
        sideset->lfi()->data(), n2e->rowptr()->data(), n2e->colidx()->data(), element_type,
        mesh->n_elements(), elements, seed_side, selected->data(),
        [&](const ptrdiff_t prev_side, const ptrdiff_t next_side) {
          const element_idx_t pe = surf_parent[prev_side];
          const element_idx_t ne = surf_parent[next_side];
          const i16 ps = surf_lfi[prev_side];
          const i16 ns = surf_lfi[next_side];

          real_t n0[3] = {0, 0, 0};
          real_t n1[3] = {0, 0, 0};
          if (!compute_side_normal(dim, element_type, elements, points, pe, ps,
                                   n0)) {
            return false;
          }
          if (!compute_side_normal(dim, element_type, elements, points, ne, ns,
                                   n1)) {
            return false;
          }

          const real_t cos_angle =
              std::fabs(n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]);
          return cos_angle > angle_threshold;
        });


    if (err != SMESH_SUCCESS) {
      SMESH_ERROR("Failed to select sideset!\n");
      return SMESH_FAILURE;
    }

    ptrdiff_t n_selected = 0;
    for (ptrdiff_t i = 0; i < n_surf; ++i) {
      n_selected += (mask_get(i, selected->data()) != 0);
    }

    auto out_parent = create_host_buffer<element_idx_t>(n_selected);
    auto out_lfi = create_host_buffer<i16>(n_selected);

    ptrdiff_t ins = 0;
    for (ptrdiff_t i = 0; i < n_surf; ++i) {
      if (!mask_get(i, selected->data())) {
        continue;
      }
      out_parent->data()[ins] = surf_parent[i];
      out_lfi->data()[ins] = surf_lfi[i];
      ++ins;
    }

    auto out_ss =
        Sideset::create(mesh->comm(), out_parent, out_lfi, sideset->block_id());
    if (out_ss->write(output_folder) != SMESH_SUCCESS) {
      SMESH_ERROR("Failed to write output sideset\n");
      return SMESH_FAILURE;
    }
  }

  return ret;
}

// Reference implementation:

// #include <math.h>
// // (avoid hash maps; keep data in CSR/arrays per original approach)
// #include <stdint.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <sys/stat.h>

// #include "array_dtof.h"
// #include "matrixio_array.h"
// #include "matrixio_crs.h"
// #include "utils.h"

// #include "crs_graph.h"
// #include "read_mesh.h"
// #include "sfem_base.h"
// #include "sfem_mesh_write.h"

// #include "sfem_defs.h"

// #include "argsort.h"

// #include "sfem_API.hpp"

// static SFEM_INLINE void normal(real_t u[3], real_t v[3], real_t *n) {
//     const real_t u_len = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
//     const real_t v_len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

//     for (int d = 0; d < 3; ++d) {
//         u[d] /= u_len;
//         v[d] /= v_len;
//     }

//     n[0] = (u[1] * v[2]) - (u[2] * v[1]);
//     n[1] = (u[2] * v[0]) - (u[0] * v[2]);
//     n[2] = (u[0] * v[1]) - (u[1] * v[0]);

//     const real_t n_len = sqrt((n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]));

//     for (int d = 0; d < 3; ++d) {
//         n[d] /= n_len;
//     }
// }

// static SFEM_INLINE void normal2(real_t p0[2], real_t p1[2], real_t *n) {
//     n[0] = -p1[1] + p0[1];
//     n[1] = p1[0] - p0[0];

//     const real_t len_n = sqrt(n[0] * n[0] + n[1] * n[1]);
//     n[0] /= len_n;
//     n[1] /= len_n;
// }

// int main(int argc, char *argv[]) {
//     sfem::Context context(argc, argv);
//     {
//         auto comm = context.communicator();

//         int rank, size;
//         MPI_Comm_rank(comm->get(), &rank);
//         MPI_Comm_size(comm->get(), &size);

//         if (argc != 7) {
//             if (!rank) {
//                 fprintf(stderr, "usage: %s <folder> <x> <y> <z>
//                 <angle_threshold> <output_folder>\n", argv[0]);
//             }

//             return EXIT_FAILURE;
//         }

//         const char  *path_mesh       = argv[1];
//         const geom_t roi[3]          = {(geom_t)atof(argv[2]),
//         (geom_t)atof(argv[3]), (geom_t)atof(argv[4])}; const geom_t
//         angle_threshold = atof(argv[5]); std::string  output_folder   =
//         argv[6];

//         int SFEM_DEBUG = 0;
//         SFEM_READ_ENV(SFEM_DEBUG, atoi);

//         if (!rank) {
//             fprintf(stderr,
//                     "%s %s %g %g %g %g %s\n",
//                     argv[0],
//                     path_mesh,
//                     (double)roi[0],
//                     (double)roi[1],
//                     (double)roi[2],
//                     (double)angle_threshold,
//                     output_folder.c_str());
//         }

//         double tick = MPI_Wtime();

//         ///////////////////////////////////////////////////////////////////////////////
//         // Read data
//         ///////////////////////////////////////////////////////////////////////////////

//         auto mesh = sfem::Mesh::create_from_file(comm, path_mesh);

//         ///////////////////////////////////////////////////////////////////////////////
//         // Extract buffers and values
//         ///////////////////////////////////////////////////////////////////////////////

//         const ptrdiff_t n_elements = mesh->n_elements();
//         const ptrdiff_t n_nodes    = mesh->n_nodes();
//         auto            elements   = mesh->elements()->data();
//         auto            points     = mesh->points()->data();

//         ///////////////////////////////////////////////////////////////////////////////
//         // Create skin sideset
//         ///////////////////////////////////////////////////////////////////////////////

//         // Reduce cost of computation by exploiting low-order representation
//         enum ElemType element_type_hack = mesh->element_type();
//         switch (element_type_hack) {
//             case TRI6: {
//                 element_type_hack = TRI3;
//                 break;
//             }
//             case TET10: {
//                 element_type_hack = TET4;
//                 break;
//             }
//             case EDGE3: {
//                 element_type_hack = EDGE2;
//                 break;
//             }
//             default:
//                 break;
//         }

//         const int nsxe = elem_num_sides(element_type_hack);

//         enum ElemType st   = side_type(element_type_hack);
//         const int     nnxs = elem_num_nodes(st);
//         const int     dim  = mesh->spatial_dimension();
//         const int     ns   = elem_num_sides(element_type_hack);

//         ptrdiff_t      n_surf_elements = 0;
//         element_idx_t *parent_buff     = 0;
//         int16_t       *side_idx_buff   = 0;
//         element_idx_t *table_buff      = 0;

//         create_element_adj_table(n_elements, n_nodes, element_type_hack,
//         elements, &table_buff);

//         if (extract_sideset_from_adj_table(
//                     element_type_hack, n_elements, table_buff,
//                     &n_surf_elements, &parent_buff, &side_idx_buff) !=
//                     SFEM_SUCCESS) {
//             SFEM_ERROR("Failed to extract
//             extract_sideset_from_adj_table!\n");
//         }

//         auto parent              =
//         sfem::manage_host_buffer<element_idx_t>(n_surf_elements,
//         parent_buff); auto side_idx            =
//         sfem::manage_host_buffer<int16_t>(n_surf_elements, side_idx_buff);
//         auto table               =
//         sfem::manage_host_buffer<element_idx_t>(n_elements * nsxe,
//         table_buff); auto element_mapping_ptr =
//         sfem::create_host_buffer<element_idx_t>(n_elements + 1);

//         auto local_side_table = sfem::create_host_buffer<int>(nsxe * nnxs);
//         fill_local_side_table(element_type_hack, local_side_table->data());

//         ///////////////////////////////////////////////////////////////////////////////
//         // Find approximately closest elemenent
//         ///////////////////////////////////////////////////////////////////////////////

//         element_idx_t *surf_parents = parent->data();
//         int16_t       *surf_idx     = side_idx->data();
//         auto           lst          = local_side_table->data();

//         element_idx_t closest_element = SFEM_ELEMENT_IDX_INVALID;
//         int16_t       closest_side    = -1;
//         real_t        closest_sq_dist = 1000000;

//         auto emap_ptr = element_mapping_ptr->data();

// #pragma omp parallel for
//         for (ptrdiff_t e = 0; e < n_surf_elements; ++e) {
//             geom_t        element_sq_dist = 1000000;
//             element_idx_t sp              = surf_parents[e];
//             int16_t       s               = surf_idx[e];

// #pragma omp atomic update
//             emap_ptr[sp + 1]++;

//             geom_t barycenter[3] = {0., 0., 0.};
//             for (int d = 0; d < dim; ++d) {
//                 for (int n = 0; n < nnxs; n++) {
//                     idx_t node = elements[lst[s * nnxs + n]][sp];
//                     barycenter[d] += points[d][node];
//                 }
//                 barycenter[d] /= nnxs;
//             }

//             real_t sq_dist = 0.;
//             for (int d = 0; d < dim; ++d) {
//                 const real_t m_x   = barycenter[d];
//                 const real_t roi_x = roi[d];
//                 const real_t diff  = m_x - roi_x;
//                 sq_dist += diff * diff;
//             }

//             element_sq_dist = MIN(element_sq_dist, sq_dist);

// #pragma omp critical
//             {
//                 if (element_sq_dist < closest_sq_dist) {
//                     closest_sq_dist = element_sq_dist;
//                     closest_element = e;
//                     closest_side    = s;
//                 }
//             }
//         }

//         for (ptrdiff_t i = 0; i < n_elements; i++) {
//             emap_ptr[i + 1] += emap_ptr[i];
//         }

//         ptrdiff_t nmaps               = emap_ptr[n_elements];
//         auto      element_mapping_idx =
//         sfem::create_host_buffer<ptrdiff_t>(nmaps);

//         // assert(nmaps < table->size());

//         auto emap_idx = element_mapping_idx->data();

//         {
//             auto book_keeping =
//             sfem::create_host_buffer<element_idx_t>(n_elements); auto bk =
//             book_keeping->data(); for (ptrdiff_t e = 0; e < n_surf_elements;
//             ++e) {
//                 element_idx_t sp                  = surf_parents[e];
//                 emap_idx[emap_ptr[sp] + bk[sp]++] = e;
//             }
//         }

// #ifndef NDEBUG
//         for (ptrdiff_t i = 0; i < nmaps; i++) {
//             assert(emap_idx[i] < table->size());
//         }
// #endif

//         if (closest_element == SFEM_IDX_INVALID) {
//             SFEM_ERROR("Invalid set up! for mesh #nelements %ld #nodes
//             %ld\n", n_elements, n_nodes);
//         }

//         auto adj      = table->data();
//         auto selected = sfem::create_host_buffer<uint8_t>(n_surf_elements);
//         auto eselect  = selected->data();

//         ptrdiff_t size_queue    = (n_surf_elements + 1);
//         auto      element_queue =
//         sfem::create_host_buffer<ptrdiff_t>(size_queue); auto      equeue =
//         element_queue->data();

//         equeue[0] = closest_element;
//         for (ptrdiff_t e = 1; e < size_queue; ++e) {
//             equeue[e] = SFEM_PTRDIFF_INVALID;
//         }

//         // Next slot
//         ptrdiff_t next_slot = 1;

//         ///////////////////////////////////////////////////////////////////////////////
//         // Create marker for different faces based on dihedral angles
//         ///////////////////////////////////////////////////////////////////////////////

//         // Build node-to-element CSR once, to walk boundary by shared
//         nodes/edges count_t       *n2e_ptr  = nullptr; element_idx_t *n2e_el
//         = nullptr; const int      nxe_full =
//         elem_num_nodes(mesh->element_type()); if (build_n2e(n_elements,
//         n_nodes, nxe_full, elements, &n2e_ptr, &n2e_el) != SFEM_SUCCESS) {
//             SFEM_ERROR("Failed to build node->element incidence\n");
//         }

//         if (dim == 2) {
//             for (ptrdiff_t q = 0; equeue[q] >= 0; q = (q + 1) % size_queue) {
//                 const ptrdiff_t     e  = equeue[q];
//                 const element_idx_t sp = surf_parents[e];
//                 const int16_t       s  = surf_idx[e];
//                 if (eselect[e]) {
//                     equeue[q] = SFEM_PTRDIFF_INVALID;
//                     continue;
//                 }

//                 // face normal
//                 real_t      n2[2];
//                 const idx_t e_nodes[2] = {elements[lst[s * nnxs + 0]][sp],
//                 elements[lst[s * nnxs + 1]][sp]};
//                 {
//                     real_t p0[2], p1[2];
//                     for (int d = 0; d < 2; ++d) {
//                         p0[d] = points[d][e_nodes[0]];
//                         p1[d] = points[d][e_nodes[1]];
//                     }
//                     normal2(p0, p1, n2);
//                 }

//                 // candidate neighbors: faces on elements incident to either
//                 node for (int vn = 0; vn < 2; ++vn) {
//                     const idx_t   vtx = e_nodes[vn];
//                     const count_t beg = n2e_ptr[vtx];
//                     const count_t end = n2e_ptr[vtx + 1];
//                     for (count_t it = beg; it < end; ++it) {
//                         const element_idx_t esp = n2e_el[it];
//                         for (ptrdiff_t k = emap_ptr[esp]; k < emap_ptr[esp +
//                         1]; ++k) {
//                             const element_idx_t ne = emap_idx[k];
//                             if (ne == e || ne == SFEM_ELEMENT_IDX_INVALID ||
//                             eselect[ne]) continue; const int16_t ns         =
//                             surf_idx[ne]; const idx_t   n_nodes[2] =
//                             {elements[lst[ns * nnxs + 0]][esp],
//                             elements[lst[ns * nnxs + 1]][esp]};
//                             // share at least 1 node
//                             int shared = (n_nodes[0] == e_nodes[0]) ||
//                             (n_nodes[0] == e_nodes[1]) || (n_nodes[1] ==
//                             e_nodes[0]) ||
//                                          (n_nodes[1] == e_nodes[1]);
//                             if (!shared) continue;
//                             // angle check
//                             real_t p0a[2], p1a[2], na2[2];
//                             for (int d = 0; d < 2; ++d) {
//                                 p0a[d] = points[d][n_nodes[0]];
//                                 p1a[d] = points[d][n_nodes[1]];
//                             }
//                             normal2(p0a, p1a, na2);
//                             const real_t cos_angle = fabs(n2[0] * na2[0] +
//                             n2[1] * na2[1]); if (cos_angle > angle_threshold)
//                             {
//                                 equeue[next_slot++ % size_queue] = ne;
//                             }
//                         }
//                     }
//                 }

//                 eselect[e] = 1;
//                 equeue[q]  = SFEM_PTRDIFF_INVALID;
//             }
//         } else {
//             for (ptrdiff_t q = 0; equeue[q] >= 0; q = (q + 1) % size_queue) {
//                 const ptrdiff_t     e  = equeue[q];
//                 const element_idx_t sp = surf_parents[e];
//                 const int16_t       s  = surf_idx[e];
//                 if (eselect[e]) {
//                     equeue[q] = SFEM_PTRDIFF_INVALID;
//                     continue;
//                 }

//                 // face normal
//                 real_t      n3[3];
//                 const idx_t e_nodes[3] = {
//                         elements[lst[s * nnxs + 0]][sp], elements[lst[s *
//                         nnxs + 1]][sp], elements[lst[s * nnxs + 2]][sp]};
//                 {
//                     real_t u[3], v[3];
//                     for (int d = 0; d < 3; ++d) {
//                         u[d] = points[d][e_nodes[1]] - points[d][e_nodes[0]];
//                         v[d] = points[d][e_nodes[2]] - points[d][e_nodes[0]];
//                     }
//                     normal(u, v, n3);
//                 }

//                 // candidate neighbors: faces on elements incident to any
//                 face node for (int vn = 0; vn < nnxs; ++vn) {
//                     const idx_t   vtx = elements[lst[s * nnxs + vn]][sp];
//                     const count_t beg = n2e_ptr[vtx];
//                     const count_t end = n2e_ptr[vtx + 1];
//                     for (count_t it = beg; it < end; ++it) {
//                         const element_idx_t esp = n2e_el[it];
//                         for (ptrdiff_t k = emap_ptr[esp]; k < emap_ptr[esp +
//                         1]; ++k) {
//                             const element_idx_t ne = emap_idx[k];
//                             if (ne == e || ne == SFEM_ELEMENT_IDX_INVALID ||
//                             eselect[ne]) continue; const int16_t ns =
//                             surf_idx[ne];

//                             // count shared nodes (edge adjacency requires
//                             >=2) int shared = 0; for (int a = 0; a < nnxs;
//                             ++a) {
//                                 for (int b = 0; b < nnxs; ++b) {
//                                     shared += (elements[lst[s * nnxs +
//                                     a]][sp] == elements[lst[ns * nnxs +
//                                     b]][esp]);
//                                 }
//                             }
//                             if (shared < 2) continue;

//                             const idx_t n_nodes[3] = {elements[lst[ns * nnxs
//                             + 0]][esp],
//                                                       elements[lst[ns * nnxs
//                                                       + 1]][esp],
//                                                       elements[lst[ns * nnxs
//                                                       + 2]][esp]};

//                             // angle check
//                             real_t ua[3], va[3], na[3];
//                             for (int d = 0; d < 3; ++d) {
//                                 ua[d] = points[d][n_nodes[1]] -
//                                 points[d][n_nodes[0]]; va[d] =
//                                 points[d][n_nodes[2]] -
//                                 points[d][n_nodes[0]];
//                             }
//                             normal(ua, va, na);
//                             const real_t cos_angle = fabs(n3[0] * na[0] +
//                             n3[1] * na[1] + n3[2] * na[2]); if (cos_angle >
//                             angle_threshold) {
//                                 equeue[next_slot++ % size_queue] = ne;
//                             }
//                         }
//                     }
//                 }

//                 eselect[e] = 1;
//                 equeue[q]  = SFEM_PTRDIFF_INVALID;
//             }
//         }

//         free(n2e_ptr);
//         free(n2e_el);

//         ///////////////////////////////////////////////////////////////////////////////
//         // Create sidesets
//         ///////////////////////////////////////////////////////////////////////////////

//         ptrdiff_t n_selected = 0;
//         for (ptrdiff_t i = 0; i < n_surf_elements; i++) {
//             n_selected += eselect[i] == 1;
//         }

//         auto sideset_parent =
//         sfem::create_host_buffer<element_idx_t>(n_selected); auto sideset_lfi
//         = sfem::create_host_buffer<int16_t>(n_selected);

//         auto ssp   = sideset_parent->data();
//         auto sslfi = sideset_lfi->data();

//         for (ptrdiff_t e = 0, n_inserted = 0; e < n_surf_elements; ++e) {
//             if (eselect[e]) {
//                 sslfi[n_inserted] = surf_idx[e];
//                 ssp[n_inserted++] = surf_parents[e];
//             }
//         }

//         sfem::create_directory(output_folder.c_str());
//         sideset_parent->to_file((output_folder + "/parent.raw").c_str());
//         sideset_lfi->to_file((output_folder + "/lfi.int16.raw").c_str());

//         if (SFEM_DEBUG) {
//             printf("Extraced %ld/%ld surface elements\n",
//             long(sideset_parent->size()), long(parent->size())); auto
//             debug_elements =
//             sfem::create_host_buffer<idx_t>(elem_num_nodes(side_type(mesh->element_type())),
//             n_selected);

//             if (extract_surface_from_sideset(mesh->element_type(),
//                                              elements,
//                                              n_selected,
//                                              sideset_parent->data(),
//                                              sideset_lfi->data(),
//                                              debug_elements->data()) !=
//                                              SFEM_SUCCESS) {
//                 SFEM_ERROR("Unable to extract surface from sideset!\n");
//             }

//             sfem::create_directory((output_folder + "/surf").c_str());
//             debug_elements->to_files((output_folder +
//             "/surf/i%d.raw").c_str());
//         }

//         double tock = MPI_Wtime();

//         if (!rank) {
//             printf("----------------------------------------\n");
//             printf("TTS:\t\t\t%g seconds\n", tock - tick);
//         }
//     }

//     return SFEM_SUCCESS;
// }
