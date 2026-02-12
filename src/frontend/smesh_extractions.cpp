#include "smesh_extractions.hpp"

#include "smesh_mesh.hpp"

namespace smesh {

// template <typename idx_t, typename count_t, typename element_idx_t>
// int extract_disconnected_faces(
//     const enum ElemType element_type, const ptrdiff_t nelements,
//     const ptrdiff_t nnodes,
//     idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elements,
//     const ptrdiff_t n_sharp_edges, const idx_t *const SMESH_RESTRICT e0,
//     const idx_t *const SMESH_RESTRICT e1,
//     ptrdiff_t *out_n_disconnected_elements,
//     element_idx_t **out_disconnected_elements);

SharedBuffer<idx_t *> extract_sharp_edges(Mesh &mesh,
                                          const geom_t cos_angle_threshold) {
  idx_t *out_e0 = nullptr;
  idx_t *out_e1 = nullptr;
  ptrdiff_t out_n_sharp_edges = 0;
  auto n2n = mesh.node_to_node_graph(); // TODO: see if we can use the upper
                                        // triangular graph instead
  extract_sharp_edges(mesh.element_type(), mesh.n_elements(),
                      mesh.elements()->data(), mesh.n_nodes(),
                      mesh.points()->data(), n2n->rowptr()->data(),
                      n2n->colidx()->data(), cos_angle_threshold,
                      &out_n_sharp_edges, &out_e0, &out_e1);

  idx_t **edges = (idx_t **)malloc(2 * sizeof(idx_t *));
  edges[0] = out_e0;
  edges[1] = out_e1;
  return manage_host_buffer<idx_t>(2, out_n_sharp_edges, edges);
}

SharedBuffer<idx_t> extract_sharp_corners(const ptrdiff_t n_nodes,
                                            SharedBuffer<idx_t *> &sharp_edges,
                                            const bool edge_clean_up) {
  idx_t *out_corners = nullptr;
  ptrdiff_t out_ncorners = 0;
  extract_sharp_corners(n_nodes, sharp_edges->extent(1), sharp_edges->data()[0],
                        sharp_edges->data()[1], &out_ncorners, &out_corners,
                        edge_clean_up);

  return manage_host_buffer<idx_t>(out_ncorners, out_corners);
}

SharedBuffer<element_idx_t> extract_disconnected_faces(Mesh &mesh,
                                                        SharedBuffer<idx_t *> &sharp_edges) {
  ptrdiff_t out_n_disconnected_faces = 0;
  element_idx_t *out_disconnected_faces = nullptr;
  extract_disconnected_faces(mesh.element_type(), mesh.n_elements(),
                             mesh.n_nodes(), mesh.elements()->data(),
                             sharp_edges->extent(1), sharp_edges->data()[0],
                             sharp_edges->data()[1], &out_n_disconnected_faces,
                             &out_disconnected_faces);
  return manage_host_buffer<element_idx_t>(out_n_disconnected_faces, out_disconnected_faces);
}

} // namespace smesh
