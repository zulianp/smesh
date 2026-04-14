#include "smesh_extractions.hpp"
#include "smesh_alloc.hpp"

#include "smesh_mesh.hpp"

namespace smesh {

SharedBuffer<idx_t *> extract_sharp_edges(Mesh &mesh,
                                          const geom_t cos_angle_threshold) {
  if (mesh.n_blocks() > 1) {
    SMESH_ERROR("extract_sharp_edges is not supported for multiblock meshes");
    return nullptr;
  }

  const block_idx_t block_id = 0;
  idx_t *out_e0 = nullptr;
  idx_t *out_e1 = nullptr;
  ptrdiff_t out_n_sharp_edges = 0;
  auto n2n = mesh.node_to_node_graph(); // TODO: see if we can use the upper
                                        // triangular graph instead
  extract_sharp_edges(mesh.element_type(block_id), mesh.n_elements(block_id),
                      mesh.elements(block_id)->data(), mesh.n_nodes(),
                      mesh.points()->data(), n2n->rowptr()->data(),
                      n2n->colidx()->data(), cos_angle_threshold,
                      &out_n_sharp_edges, &out_e0, &out_e1);

  idx_t **edges = (idx_t **)SMESH_ALLOC(2 * sizeof(idx_t *));
  edges[0] = out_e0;
  edges[1] = out_e1;
  return manage_host_buffer<idx_t>(2, out_n_sharp_edges, edges);
}

SharedBuffer<idx_t> extract_sharp_corners(const ptrdiff_t n_nodes,
                                          SharedBuffer<idx_t *> &sharp_edges,
                                          const bool edge_clean_up) {
  idx_t *out_corners = nullptr;
  ptrdiff_t out_ncorners = 0;
  ptrdiff_t out_n_sharp_edges =
      extract_sharp_corners(n_nodes, sharp_edges->extent(1),
                            sharp_edges->data()[0], sharp_edges->data()[1],
                            &out_ncorners, &out_corners, edge_clean_up);

  if (edge_clean_up && out_n_sharp_edges != sharp_edges->extent(1)) {
    auto cleaned_edges = create_host_buffer<idx_t>(2, out_n_sharp_edges);
    for (ptrdiff_t i = 0; i < out_n_sharp_edges; ++i) {
      cleaned_edges->data()[0][i] = sharp_edges->data()[0][i];
      cleaned_edges->data()[1][i] = sharp_edges->data()[1][i];
    }

    sharp_edges = cleaned_edges;
  }

  return manage_host_buffer<idx_t>(out_ncorners, out_corners);
}

SharedBuffer<element_idx_t>
extract_disconnected_faces(Mesh &mesh, SharedBuffer<idx_t *> &sharp_edges) {
  if (mesh.n_blocks() > 1) {
    SMESH_ERROR(
        "extract_disconnected_faces is not supported for multiblock meshes");
    return nullptr;
  }

  const block_idx_t block_id = 0;
  ptrdiff_t out_n_disconnected_faces = 0;
  element_idx_t *out_disconnected_faces = nullptr;
  extract_disconnected_faces(
      mesh.element_type(block_id), mesh.n_elements(block_id), mesh.n_nodes(),
      mesh.elements(block_id)->data(), sharp_edges->extent(1),
      sharp_edges->data()[0], sharp_edges->data()[1], &out_n_disconnected_faces,
      &out_disconnected_faces);
  return manage_host_buffer<element_idx_t>(out_n_disconnected_faces,
                                           out_disconnected_faces);
}

} // namespace smesh
