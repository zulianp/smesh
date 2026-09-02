#include "smesh_mesh_reorder.hpp"

#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_ops.hpp"
#include "smesh_reorder.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sfc.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sort.hpp"
#include "smesh_tracer.hpp"

#include <cstring>
#include <functional>
#include <map>
#include <vector>

namespace smesh {

class SFC::Impl {
public:
  std::string ordering_type;
};

SFC::SFC(const std::string &ordering_type) : impl_(std::make_unique<Impl>()) {
  impl_->ordering_type = ordering_type;
}
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

int SFC::reorder(Mesh &mesh, const std::vector<std::shared_ptr<Sideset>> &sidesets) {
  SMESH_TRACE_SCOPE("SFC::reorder");
  if (mesh.n_blocks() > 1) {
    SMESH_ERROR("SFC::reorder is not supported for multiblock meshes");
    return SMESH_FAILURE;
  }

  const block_idx_t block_id = 0;

  std::map<std::string,
           std::function<int(
               const ptrdiff_t, const geom_t *const SMESH_RESTRICT,
               const geom_t *const SMESH_RESTRICT,
               const geom_t *const SMESH_RESTRICT, u32 *const SMESH_RESTRICT)>>
      encode_functions = {
          {"morton3",
           [](const ptrdiff_t n_points, const geom_t *const SMESH_RESTRICT x,
              const geom_t *const SMESH_RESTRICT y,
              const geom_t *const SMESH_RESTRICT z,
              u32 *const SMESH_RESTRICT encoding) {
             return encode_morton3<geom_t>(n_points, x, y, z, encoding);
           }},
          {"hilbert3",
           [](const ptrdiff_t n_points, const geom_t *const SMESH_RESTRICT x,
              const geom_t *const SMESH_RESTRICT y,
              const geom_t *const SMESH_RESTRICT z,
              u32 *const SMESH_RESTRICT encoding) {
             return encode_hilbert3<geom_t>(n_points, x, y, z, encoding);
           }},
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
           }},
          {"random3", encode_random3<geom_t>}};

  auto iter = encode_functions.find(impl_->ordering_type);
  if (iter == encode_functions.end()) {
    SMESH_ERROR("Invalid ordering type");
    return SMESH_FAILURE;
  }

  int spatial_dim = mesh.spatial_dimension();
  int nxe = mesh.n_nodes_per_element(block_id);
  const ptrdiff_t n_elements = mesh.n_elements(block_id);
  const ptrdiff_t n_nodes = mesh.n_nodes();

  idx_t *const *const elems = mesh.elements(block_id)->data();
  geom_t *const *const pts  = mesh.points()->data();

  auto b = create_host_buffer<geom_t>(3, n_elements);
  geom_t **d_b = b->data();
  barycenters(nxe, n_elements, elems, spatial_dim, pts, d_b);

  auto encoding = create_host_buffer<u32>(n_elements);
  u32 *d_enc = encoding->data();
  SMESH_CATCH(iter->second(n_elements, d_b[0], d_b[1], d_b[2], d_enc));

  auto idx = create_host_buffer<idx_t>(n_elements);
  idx_t *d_idx = idx->data();
  argsort(n_elements, d_enc, d_idx);

  SMESH_CATCH(mesh_block_reorder(nxe, n_elements, elems, d_idx, elems));

  const bool remap_arg = !sidesets.empty();
  const bool remap_reg = !mesh.sidesets().empty();
  if (remap_arg || remap_reg) {
    auto old_to_new = create_host_buffer<element_idx_t>(n_elements);
    element_idx_t *d_otn = old_to_new->data();
    for (ptrdiff_t neu = 0; neu < n_elements; ++neu) {
      d_otn[d_idx[neu]] = static_cast<element_idx_t>(neu);
    }
    if (remap_arg) {
      if (remap_sidesets(sidesets, block_id, d_otn, n_elements) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
    }
    if (mesh.remap_registered_sidesets(block_id, d_otn, n_elements, sidesets) !=
        SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }

  if (is_semistructured_type(mesh.element_type(block_id))) {
    return semistructured_hierarchical_renumbering(mesh.element_type(block_id),
                                                   semistructured_level(mesh),
                                                   n_nodes,
                                                   mesh.elements(block_id),
                                                   mesh.points(),
                                                   false);
  }

  auto n2n_scatter = create_host_buffer<idx_t>(n_nodes);
  idx_t *d_n2n = n_nodes > 0 ? n2n_scatter->data() : nullptr;
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    d_n2n[i] = invalid_idx<idx_t>();
  }

  idx_t next_node_idx = 0;
  SMESH_CATCH(mesh_block_renumber_element_nodes<idx_t>(
      nxe, n_elements, elems, &next_node_idx, d_n2n));

  auto coords = create_host_buffer<geom_t>(n_nodes);
  geom_t *d_coords = n_nodes > 0 ? coords->data() : nullptr;
  if (n_nodes > 0) {
    memcpy(d_coords, pts[0], n_nodes * sizeof(geom_t));
    SMESH_CATCH(reorder_scatter(n_nodes, d_n2n, d_coords, pts[0]));

    if (spatial_dim > 1) {
      memcpy(d_coords, pts[1], n_nodes * sizeof(geom_t));
      SMESH_CATCH(reorder_scatter(n_nodes, d_n2n, d_coords, pts[1]));
    }

    if (spatial_dim > 2) {
      memcpy(d_coords, pts[2], n_nodes * sizeof(geom_t));
      SMESH_CATCH(reorder_scatter(n_nodes, d_n2n, d_coords, pts[2]));
    }
  }

  return SMESH_SUCCESS;
}
} // namespace smesh
