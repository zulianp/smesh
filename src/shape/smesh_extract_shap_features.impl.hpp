#ifndef SMESH_EXTRACT_SHAP_FEATURES_IMPL_HPP
#define SMESH_EXTRACT_SHAP_FEATURES_IMPL_HPP

#include "smesh_base.hpp"
#include "smesh_extract_shap_features.hpp"

#include "smesh_common.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_types.hpp"

namespace smesh {

template <typename idx_t, typename geom_t, typename count_t>
int extract_sharp_edges(
    const enum ElemType element_type, const ptrdiff_t nelements,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t nnodes,
    geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
    // CRS-graph (node to node)
    const count_t *const SMESH_RESTRICT rowptr,
    const idx_t *const SMESH_RESTRICT colidx, const geom_t angle_threshold,
    ptrdiff_t *out_n_sharp_edges, idx_t **out_e0, idx_t **out_e1) {
  const count_t nedges = rowptr[nnodes];

  geom_t *normal[3];
  for (int d = 0; d < 3; d++) {
    normal[d] = (geom_t *)calloc(nedges, sizeof(geom_t));
  }

  const int nxe = elem_num_nodes(element_type);
  count_t *opposite = (count_t *)malloc(nedges * sizeof(count_t));
  {
    // Opposite edge index
    // #pragma omp parallel for
    for (ptrdiff_t i = 0; i < nnodes; i++) {
      const count_t begin = rowptr[i];
      const count_t extent = rowptr[i + 1] - begin;
      const idx_t *cols = &colidx[begin];

      for (count_t k = 0; k < extent; k++) {
        const idx_t o = cols[k];
        if (i > o)
          continue;

        const count_t o_begin = rowptr[o];
        const count_t o_extent = rowptr[o + 1] - o_begin;
        const idx_t *o_cols = &colidx[o_begin];

        for (count_t o_k = 0; o_k < o_extent; o_k++) {
          if (i == o_cols[o_k]) {
            opposite[begin + k] = o_begin + o_k;
            opposite[o_begin + o_k] = begin + k;
            break;
          }
        }
      }
    }
  }

  {
    // Compute normals
     // #pragma omp parallel for
    for (ptrdiff_t e = 0; e < nelements; e++) {
      const idx_t i0 = elements[0][e];
      const idx_t i1 = elements[1][e];
      const idx_t i2 = elements[2][e];

      double ux = points[0][i1] - points[0][i0];
      double uy = points[1][i1] - points[1][i0];
      double uz = points[2][i1] - points[2][i0];

      double vx = points[0][i2] - points[0][i0];
      double vy = points[1][i2] - points[1][i0];
      double vz = points[2][i2] - points[2][i0];

      normalize3(&ux, &uy, &uz);
      normalize3(&vx, &vy, &vz);

      double nx = uy * vz - uz * vy;
      double ny = uz * vx - ux * vz;
      double nz = ux * vy - uy * vx;

      normalize3(&nx, &ny, &nz);

      for (int ln = 0; ln < nxe; ln++) {
        const int lnp1 = (ln + 1 == nxe) ? 0 : (ln + 1);

        const idx_t node_from = elements[ln][e];
        const idx_t node_to = elements[lnp1][e];

        const count_t extent = rowptr[node_from + 1] - rowptr[node_from];
        const idx_t *cols = &colidx[rowptr[node_from]];

        ptrdiff_t edge_id = invalid_idx<ptrdiff_t>();
        for (count_t k = 0; k < extent; k++) {
          if (cols[k] == node_to) {
            edge_id = rowptr[node_from] + k;
            break;
          }
        }

        SMESH_ASSERT(edge_id != invalid_idx<ptrdiff_t>());
        normal[0][edge_id] = nx;
        normal[1][edge_id] = ny;
        normal[2][edge_id] = nz;
      }
    }
  }

  ptrdiff_t n_sharp_edges = 0;
  idx_t *e0 = (idx_t *)malloc(nedges * sizeof(idx_t));
  idx_t *e1 = (idx_t *)malloc(nedges * sizeof(idx_t));

  {
    geom_t *dihedral_angle = (geom_t *)calloc(nedges, sizeof(geom_t));
    ptrdiff_t edge_count = 0;
    {
      for (ptrdiff_t i = 0; i < nnodes; i++) {
        const count_t begin = rowptr[i];
        const count_t extent = rowptr[i + 1] - begin;
        const idx_t *cols = &colidx[begin];

        for (count_t k = 0; k < extent; k++) {
          if (i >= cols[k])
            continue;

          ptrdiff_t edge_id = begin + k;
          ptrdiff_t o_edge_id = opposite[edge_id];

          SMESH_ASSERT(edge_id != o_edge_id);

          double da = dot3<double>(normal[0][edge_id], normal[1][edge_id],
                                   normal[2][edge_id], normal[0][o_edge_id],
                                   normal[1][o_edge_id], normal[2][o_edge_id]);

          // Store for minimum edge for exporting data
          dihedral_angle[edge_count] = (geom_t)da;
          e0[edge_count] = i;
          e1[edge_count] = cols[k];
          edge_count++;
        }
      }
    }

    // 1) select sharp edges create edge selection index
    // 2) create face islands index (for contact integral separation)
    // 3) export edge and face selection
    // TODO Future work: detect sharp corners

    {
      // Select edges
      for (ptrdiff_t i = 0; i < edge_count; i++) {
        if (dihedral_angle[i] <= angle_threshold) {
          e0[n_sharp_edges] = e0[i];
          e1[n_sharp_edges] = e1[i];
          dihedral_angle[n_sharp_edges] = dihedral_angle[i];
          n_sharp_edges++;
        }
      }
    }

    free(dihedral_angle);
  }

  *out_n_sharp_edges = n_sharp_edges;
  *out_e0 = e0;
  *out_e1 = e1;

  free(opposite);
  for (int d = 0; d < 3; d++) {
    free(normal[d]);
  }

  return SMESH_SUCCESS;
}

template <typename idx_t, typename count_t>
int extract_sharp_corners(const ptrdiff_t nnodes, const ptrdiff_t n_sharp_edges,
                          idx_t *const SMESH_RESTRICT e0,
                          idx_t *const SMESH_RESTRICT e1,
                          ptrdiff_t *out_ncorners, idx_t **out_corners,
                          int edge_clean_up) {
  ptrdiff_t n_corners = 0;
  idx_t *corners = 0;

  ptrdiff_t out_n_sharp_edges = n_sharp_edges;
  {
    int *incidence_count = (int *)calloc(nnodes, sizeof(int));

    for (ptrdiff_t i = 0; i < n_sharp_edges; i++) {
      incidence_count[e0[i]]++;
      incidence_count[e1[i]]++;
    }

    for (ptrdiff_t i = 0; i < nnodes; i++) {
      if (incidence_count[i] >= 3) {
        n_corners++;
      }
    }

    corners = (idx_t *)malloc(n_corners * sizeof(idx_t));
    for (ptrdiff_t i = 0, n_corners = 0; i < nnodes; i++) {
      if (incidence_count[i] >= 3) {
        corners[n_corners] = i;
        n_corners++;
      }
    }

    if (edge_clean_up) {
      out_n_sharp_edges = 0;
      for (ptrdiff_t i = 0; i < n_sharp_edges; i++) {
        if (incidence_count[e0[i]] < 3 && incidence_count[e1[i]] < 3) {
          e0[out_n_sharp_edges] = e0[i];
          e1[out_n_sharp_edges] = e1[i];

          out_n_sharp_edges++;
        }
      }
    }

    free(incidence_count);
  }

  *out_ncorners = n_corners;
  *out_corners = corners;

  return out_n_sharp_edges;
}

template <typename idx_t, typename count_t, typename element_idx_t>
int extract_disconnected_faces(
    const enum ElemType element_type, const ptrdiff_t nelements,
    const ptrdiff_t nnodes,
    idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
    const ptrdiff_t n_sharp_edges, const idx_t *const SMESH_RESTRICT e0,
    const idx_t *const SMESH_RESTRICT e1,
    ptrdiff_t *out_n_disconnected_elements,
    element_idx_t **out_disconnected_elements) {
  ptrdiff_t n_disconnected_elements = 0;
  element_idx_t *disconnected_elements = 0;
  {
    // Select unconnected faces
    short *checked = (short *)calloc(nnodes, sizeof(short));

    const int nxe = elem_num_nodes(element_type);
    for (ptrdiff_t i = 0; i < n_sharp_edges; i++) {
      checked[e0[i]] = 1;
      checked[e1[i]] = 1;
    }

    for (ptrdiff_t e = 0; e < nelements; e++) {
      short connected_to_sharp_edge = 0;
      for (int ln = 0; ln < nxe; ln++) {
        connected_to_sharp_edge += checked[elements[ln][e]];
      }

      n_disconnected_elements += connected_to_sharp_edge == 0;
    }

    disconnected_elements = (element_idx_t *)malloc(n_disconnected_elements *
                                                    sizeof(element_idx_t));

    ptrdiff_t eidx = 0;
    for (ptrdiff_t e = 0; e < nelements; e++) {
      short connected_to_sharp_edge = 0;
      for (int ln = 0; ln < nxe; ln++) {
        connected_to_sharp_edge += checked[elements[ln][e]];
      }

      if (connected_to_sharp_edge == 0) {
        disconnected_elements[eidx++] = e;
      }
    }

    free(checked);
  }

  *out_n_disconnected_elements = n_disconnected_elements;
  *out_disconnected_elements = disconnected_elements;

  return 0;
}

} // namespace smesh

#endif // SMESH_EXTRACT_SHAP_FEATURES_IMPL_HPP