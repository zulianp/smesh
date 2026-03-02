#include "smesh_dual_graph.hpp"

#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_multiblock_graph.hpp"

namespace smesh {
class DualGraph::Impl {
public:
  SharedBuffer<count_t> adj_ptr;
  SharedBuffer<element_idx_t> adj_idx;
  SharedBuffer<block_idx_t> adj_block;
};

DualGraph::DualGraph() : impl_(std::make_unique<Impl>()) {}
DualGraph::~DualGraph() {}

SharedBuffer<count_t> DualGraph::adj_ptr() { return impl_->adj_ptr; }

SharedBuffer<element_idx_t> DualGraph::adj_idx() { return impl_->adj_idx; }

std::shared_ptr<DualGraph>
DualGraph::create(const std::shared_ptr<Mesh> &mesh) {
  auto ret = std::make_shared<DualGraph>();

  if (mesh->n_blocks() > 1) {
    const block_idx_t n_blocks = mesh->n_blocks();
    std::vector<enum ElemType> element_types(n_blocks);
    std::vector<ptrdiff_t> n_elements(n_blocks);
    std::vector<idx_t **> elements(n_blocks);
    for (block_idx_t i = 0; i < n_blocks; i++) {
      element_types[i] = mesh->block(i)->element_type();
      n_elements[i] = mesh->block(i)->n_elements();
      elements[i] = mesh->block(i)->elements()->data();
    }

    count_t *n2e_ptr = 0;
    element_idx_t *n2e_idx = 0;
    block_idx_t *n2e_block = 0;
    if (create_multiblock_n2e(n_blocks, element_types.data(), n_elements.data(),
                              elements.data(), mesh->n_nodes(), &n2e_block,
                              &n2e_ptr, &n2e_idx) != SMESH_SUCCESS) {
      return nullptr;
    }

    count_t *adj_ptr = nullptr;
    element_idx_t *adj_idx = nullptr;
    block_idx_t *adj_block = nullptr;
    create_multiblock_dual_graph_from_n2e(
        n_blocks, element_types.data(), n_elements.data(), mesh->n_nodes(),
        elements.data(), n2e_ptr, n2e_idx, n2e_block, &adj_ptr, &adj_idx,
        &adj_block);

    free(n2e_ptr);
    free(n2e_idx);
    free(n2e_block);

    const ptrdiff_t n_total_elements = mesh->n_elements();

    ret->impl_->adj_ptr =
        manage_host_buffer<count_t>(n_total_elements + 1, adj_ptr);
    ret->impl_->adj_idx =
        manage_host_buffer<element_idx_t>(adj_ptr[n_total_elements], adj_idx);
    ret->impl_->adj_block =
        manage_host_buffer<block_idx_t>(adj_ptr[n_total_elements], adj_block);
  } else {
    const ptrdiff_t n_elements = mesh->n_elements();
    const ptrdiff_t n_nodes = mesh->n_nodes();
    const enum ElemType element_type = mesh->element_type(0);
    enum ElemType element_type_for_algo = element_type;

    auto elems = mesh->elements(0)->data();
    if (element_type == TET10) {
      element_type_for_algo = TET4;
    } else if (element_type == TRI6) {
      element_type_for_algo = TRI3;
    }

    count_t *adj_ptr = 0;
    element_idx_t *adj_idx = 0;
    create_dual_graph(n_elements, n_nodes, element_type_for_algo, elems,
                      &adj_ptr, &adj_idx);

    ret->impl_->adj_ptr = manage_host_buffer<count_t>(n_elements + 1, adj_ptr);
    ret->impl_->adj_idx =
        manage_host_buffer<element_idx_t>(adj_ptr[n_elements], adj_idx);
  }

  return ret;
}
} // namespace smesh