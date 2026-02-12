#include "smesh_mesh_reorder.hpp"

#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_ops.hpp"
#include "smesh_reorder.hpp"
#include "smesh_sfc.hpp"

#include <algorithm>
#include <functional>
#include <map>

namespace smesh {

class SFC::Impl {
public:
  std::string ordering_type;
};

SFC::SFC() : impl_(std::make_unique<Impl>()) {
  impl_->ordering_type = "morton3";
}
SFC::~SFC() = default;

std::shared_ptr<SFC> SFC::create_from_env() {
  auto ret = std::make_shared<SFC>();
  ret->impl_->ordering_type =
      Env::read_string("SMESH_ORDERING_TYPE", ret->impl_->ordering_type);
  return ret;
}

int SFC::reorder(Mesh &mesh) {
  std::map<std::string,
           std::function<int(
               const ptrdiff_t, const geom_t *const SMESH_RESTRICT,
               const geom_t *const SMESH_RESTRICT,
               const geom_t *const SMESH_RESTRICT, u32 *const SMESH_RESTRICT)>>
      encode_functions = {
          {"morton3", encode_morton3<geom_t>},
          {"hilbert3", encode_hilbert3<geom_t>},
          {"cartesian3",
           [&](const ptrdiff_t n_points, const geom_t *const SMESH_RESTRICT x,
               const geom_t *const SMESH_RESTRICT y,
               const geom_t *const SMESH_RESTRICT z,
               u32 *const SMESH_RESTRICT encoding) {
             int fast = 0;
             int mid = 1;
             int slow = 2;
             return encode_cartesian3<geom_t>(n_points, x, y, z, fast, mid,
                                              slow, encoding);
           }}};

  auto iter = encode_functions.find(impl_->ordering_type);
  if (iter == encode_functions.end()) {
    SMESH_ERROR("Invalid ordering type");
    return SMESH_FAILURE;
  }

  int spatial_dim = mesh.spatial_dimension();
  int nxe = mesh.n_nodes_per_element();
  const ptrdiff_t n_elements = mesh.n_elements();
  const ptrdiff_t n_nodes = mesh.n_nodes();

  auto b = create_host_buffer<geom_t>(3, n_elements);
  barycenters(nxe, n_elements, mesh.elements()->data(), spatial_dim,
              mesh.points()->data(), b->data());

  auto encoding = create_host_buffer<u32>(n_elements);
  int err = iter->second(n_elements, b->data()[0], b->data()[1], b->data()[2],
                         encoding->data());

  auto buff =
      malloc(std::max(n_elements * sizeof(idx_t), n_nodes * sizeof(geom_t)));
  auto idx = (idx_t *)buff;
  for (ptrdiff_t i = 0; i < n_elements; i++) {
    idx[i] = i;
  }

  std::sort(idx, idx + n_elements,
            [key = encoding->data()](const idx_t l, const idx_t r) {
              return key[l] < key[r];
            });

  err |= mesh_block_reorder(nxe, n_elements, mesh.elements()->data(), idx,
                            mesh.elements()->data());

  auto n2n_scatter = create_host_buffer<idx_t>(n_nodes);

  idx_t next_node_idx = 0;
  err |= mesh_block_renumber_element_nodes<idx_t>(
      nxe, n_elements, mesh.elements()->data(), &next_node_idx,
      n2n_scatter->data());

  auto coords = (geom_t *)buff;
  memcpy(coords, mesh.points()->data()[0], n_nodes * sizeof(geom_t));
  err |= reorder_scatter(n_nodes, n2n_scatter->data(), coords,
                         mesh.points()->data()[0]);

  memcpy(coords, mesh.points()->data()[1], n_nodes * sizeof(geom_t));
  err |= reorder_scatter(n_nodes, n2n_scatter->data(), coords,
                         mesh.points()->data()[1]);

  memcpy(coords, mesh.points()->data()[2], n_nodes * sizeof(geom_t));
  err |= reorder_scatter(n_nodes, n2n_scatter->data(), coords,
                         mesh.points()->data()[2]);

  free(buff);
  return err;
}
} // namespace smesh