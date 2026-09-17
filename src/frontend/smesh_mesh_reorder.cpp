#include "smesh_mesh_reorder.hpp"

#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_ops.hpp"
#include "smesh_reorder.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sfc.hpp"
#include "smesh_sideset.hpp"
#include "smesh_edgeset.hpp"
#include "smesh_sort.hpp"
#include "smesh_tracer.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#endif

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
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

namespace {

int reorder_multiblock(
    const std::string &ordering_type, Mesh &mesh,
    const std::vector<std::shared_ptr<Sideset>> &sidesets) {
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

  auto iter = encode_functions.find(ordering_type);
  if (iter == encode_functions.end()) {
    SMESH_ERROR("Invalid ordering type");
    return SMESH_FAILURE;
  }

  const int spatial_dim = mesh.spatial_dimension();
  const ptrdiff_t n_nodes = mesh.n_nodes();
  const ptrdiff_t n_total = mesh.n_elements();
  const size_t n_blocks = mesh.n_blocks();
  geom_t *const *const pts = mesh.points()->data();

  std::vector<ptrdiff_t> offset(n_blocks + 1, 0);
  std::vector<int> nxe(n_blocks);
  std::vector<idx_t **> elems(n_blocks);
  for (size_t b = 0; b < n_blocks; ++b) {
    const block_idx_t bid = static_cast<block_idx_t>(b);
    nxe[b] = mesh.n_nodes_per_element(bid);
    elems[b] = mesh.elements(bid)->data();
    offset[b + 1] = offset[b] + mesh.n_elements(bid);
  }

  if (n_total == 0) {
    return SMESH_SUCCESS;
  }

  auto bary = create_host_buffer<geom_t>(3, static_cast<size_t>(n_total));
  geom_t **d_b = bary->data();
  for (size_t b = 0; b < n_blocks; ++b) {
    const ptrdiff_t n_e = offset[b + 1] - offset[b];
    if (n_e == 0) {
      continue;
    }
    geom_t *slice[3] = {d_b[0] + offset[b], d_b[1] + offset[b],
                        d_b[2] + offset[b]};
    barycenters(nxe[b], n_e, elems[b], spatial_dim, pts, slice);
  }
  if (spatial_dim < 3) {
    std::memset(d_b[2], 0, static_cast<size_t>(n_total) * sizeof(geom_t));
  }

  auto encoding = create_host_buffer<u32>(static_cast<size_t>(n_total));
  u32 *d_enc = encoding->data();
  SMESH_CATCH(iter->second(n_total, d_b[0], d_b[1], d_b[2], d_enc));

  auto idx = create_host_buffer<idx_t>(static_cast<size_t>(n_total));
  idx_t *d_idx = idx->data();
  argsort(n_total, d_enc, d_idx);

  auto concat_block_buf =
      create_host_buffer<block_idx_t>(static_cast<size_t>(n_total));
  block_idx_t *concat_block = concat_block_buf->data();
  for (size_t b = 0; b < n_blocks; ++b) {
    const block_idx_t bid = static_cast<block_idx_t>(b);
    for (ptrdiff_t e = offset[b]; e < offset[b + 1]; ++e) {
      concat_block[e] = bid;
    }
  }

  std::vector<SharedBuffer<idx_t>> gathers(n_blocks);
  std::vector<ptrdiff_t> cursor(n_blocks, 0);
  for (size_t b = 0; b < n_blocks; ++b) {
    const ptrdiff_t n_e = offset[b + 1] - offset[b];
    if (n_e > 0) {
      gathers[b] = create_host_buffer<idx_t>(static_cast<size_t>(n_e));
    }
  }
  for (ptrdiff_t i = 0; i < n_total; ++i) {
    const ptrdiff_t orig = static_cast<ptrdiff_t>(d_idx[i]);
    const size_t b = static_cast<size_t>(concat_block[orig]);
    gathers[b]->data()[cursor[b]++] = static_cast<idx_t>(orig - offset[b]);
  }

  const bool remap_arg = !sidesets.empty();
  const bool remap_reg = !mesh.sidesets().empty() || !mesh.edgesets().empty();
  for (size_t b = 0; b < n_blocks; ++b) {
    const ptrdiff_t n_e = offset[b + 1] - offset[b];
    if (n_e == 0) {
      continue;
    }
    idx_t *d_gather = gathers[b]->data();
    SMESH_CATCH(mesh_block_reorder(nxe[b], n_e, elems[b], d_gather, elems[b]));

    if (remap_arg || remap_reg) {
      auto old_to_new =
          create_host_buffer<element_idx_t>(static_cast<size_t>(n_e));
      element_idx_t *d_otn = old_to_new->data();
      for (ptrdiff_t neu = 0; neu < n_e; ++neu) {
        d_otn[d_gather[neu]] = static_cast<element_idx_t>(neu);
      }
      const block_idx_t bid = static_cast<block_idx_t>(b);
      if (remap_arg) {
        if (remap_sidesets(sidesets, bid, d_otn, n_e) != SMESH_SUCCESS) {
          return SMESH_FAILURE;
        }
      }
      if (mesh.remap_registered_sidesets(bid, d_otn, n_e, sidesets) !=
          SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
      if (mesh.remap_registered_edgesets(bid, d_otn, n_e) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
    }
  }

  auto n2n_scatter = create_host_buffer<idx_t>(static_cast<size_t>(n_nodes));
  idx_t *d_n2n = n_nodes > 0 ? n2n_scatter->data() : nullptr;
  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    d_n2n[i] = invalid_idx<idx_t>();
  }

  idx_t next_node_idx = 0;
  std::fill(cursor.begin(), cursor.end(), 0);
  for (ptrdiff_t i = 0; i < n_total; ++i) {
    const ptrdiff_t orig = static_cast<ptrdiff_t>(d_idx[i]);
    const size_t b = static_cast<size_t>(concat_block[orig]);
    const ptrdiff_t new_local = cursor[b]++;
    idx_t *const *const be = elems[b];
    const int bn = nxe[b];
    for (int d = 0; d < bn; ++d) {
      const idx_t ii = be[d][new_local];
      if (d_n2n[ii] == invalid_idx<idx_t>()) {
        d_n2n[ii] = next_node_idx++;
      }
      be[d][new_local] = d_n2n[ii];
    }
  }

  auto coords = create_host_buffer<geom_t>(static_cast<size_t>(n_nodes));
  geom_t *d_coords = n_nodes > 0 ? coords->data() : nullptr;
  if (n_nodes > 0) {
    memcpy(d_coords, pts[0], static_cast<size_t>(n_nodes) * sizeof(geom_t));
    SMESH_CATCH(reorder_scatter(n_nodes, d_n2n, d_coords, pts[0]));

    if (spatial_dim > 1) {
      memcpy(d_coords, pts[1], static_cast<size_t>(n_nodes) * sizeof(geom_t));
      SMESH_CATCH(reorder_scatter(n_nodes, d_n2n, d_coords, pts[1]));
    }

    if (spatial_dim > 2) {
      memcpy(d_coords, pts[2], static_cast<size_t>(n_nodes) * sizeof(geom_t));
      SMESH_CATCH(reorder_scatter(n_nodes, d_n2n, d_coords, pts[2]));
    }
  }

  if (!mesh.nodesets().empty()) {
    if (mesh.remap_registered_nodesets(d_n2n, n_nodes) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }

  return SMESH_SUCCESS;
}

#ifdef SMESH_ENABLE_MPI
int encode_with_bbox(
    const std::string &ordering_type, const ptrdiff_t n,
    const geom_t *const x, const geom_t *const y, const geom_t *const z,
    const geom_t x_min, const geom_t x_max, const geom_t y_min,
    const geom_t y_max, const geom_t z_min, const geom_t z_max,
    u32 *const encoding) {
  if (ordering_type == "morton3") {
    return encode_morton3<geom_t>(n, x, y, z, x_min, x_max, y_min, y_max,
                                  z_min, z_max, encoding);
  }
  if (ordering_type == "hilbert3") {
    return encode_hilbert3<geom_t>(n, x, y, z, x_min, x_max, y_min, y_max,
                                   z_min, z_max, encoding);
  }
  if (ordering_type == "cartesian3") {
    return encode_cartesian3<geom_t>(n, x, y, z, x_min, x_max, y_min, y_max,
                                     z_min, z_max, 0, 1, 2, encoding);
  }
  if (ordering_type == "random3") {
    return encode_random3<geom_t>(n, x, y, z, encoding);
  }
  SMESH_ERROR("Invalid ordering type");
  return SMESH_FAILURE;
}

void argsort_segment(const u32 *const keys, const ptrdiff_t lo,
                     const ptrdiff_t hi, idx_t *const gather) {
  const ptrdiff_t n = hi - lo;
  if (n <= 0) {
    return;
  }
  auto tmp = create_host_buffer<idx_t>(static_cast<size_t>(n));
  idx_t *d_tmp = tmp->data();
  argsort(n, keys + lo, d_tmp);
  for (ptrdiff_t i = 0; i < n; ++i) {
    gather[lo + i] = static_cast<idx_t>(lo + d_tmp[i]);
  }
}

void permute_large_idx(const ptrdiff_t n, large_idx_t *const arr,
                       const idx_t *const gather, const ptrdiff_t gather_off,
                       const ptrdiff_t gather_base) {
  if (n <= 0 || !arr) {
    return;
  }
  auto tmp = create_host_buffer<large_idx_t>(static_cast<size_t>(n));
  large_idx_t *d_tmp = tmp->data();
  memcpy(d_tmp, arr, static_cast<size_t>(n) * sizeof(large_idx_t));
  for (ptrdiff_t i = 0; i < n; ++i) {
    arr[i] = d_tmp[static_cast<ptrdiff_t>(gather[gather_off + i]) - gather_base];
  }
}

int reorder_distributed(
    const std::string &ordering_type, Mesh &mesh,
    const std::vector<std::shared_ptr<Sideset>> &sidesets) {
  if (ordering_type != "morton3" && ordering_type != "hilbert3" &&
      ordering_type != "cartesian3" && ordering_type != "random3") {
    SMESH_ERROR("Invalid ordering type");
    return SMESH_FAILURE;
  }

  const int spatial_dim = mesh.spatial_dimension();
  const ptrdiff_t n_total = mesh.n_elements();
  const size_t n_blocks = mesh.n_blocks();
  geom_t *const *const pts = mesh.points()->data();
  MPI_Comm comm = mesh.comm()->get();

  std::vector<ptrdiff_t> offset(n_blocks + 1, 0);
  std::vector<int> nxe(n_blocks);
  std::vector<idx_t **> elems(n_blocks);
  std::vector<ptrdiff_t> n_ons(n_blocks, 0);
  std::vector<ptrdiff_t> n_owned(n_blocks, 0);
  std::vector<ptrdiff_t> n_ghosts(n_blocks, 0);
  for (size_t b = 0; b < n_blocks; ++b) {
    const block_idx_t bid = static_cast<block_idx_t>(b);
    auto block = mesh.block(b);
    nxe[b] = mesh.n_nodes_per_element(bid);
    elems[b] = mesh.elements(bid)->data();
    const ptrdiff_t n_e = mesh.n_elements(bid);
    offset[b + 1] = offset[b] + n_e;
    n_owned[b] = block->n_elements_owned();
    n_ghosts[b] = block->n_elements_ghosts();
    if (n_owned[b] == 0 && n_ghosts[b] == 0) {
      n_owned[b] = n_e;
    }
    n_ons[b] = n_owned[b] - block->n_elements_shared();
  }

  auto bary = n_total > 0
                  ? create_host_buffer<geom_t>(3, static_cast<size_t>(n_total))
                  : nullptr;
  geom_t **d_b = bary ? bary->data() : nullptr;
  if (n_total > 0) {
    for (size_t b = 0; b < n_blocks; ++b) {
      const ptrdiff_t n_e = offset[b + 1] - offset[b];
      if (n_e == 0) {
        continue;
      }
      geom_t *slice[3] = {d_b[0] + offset[b], d_b[1] + offset[b],
                          d_b[2] + offset[b]};
      barycenters(nxe[b], n_e, elems[b], spatial_dim, pts, slice);
    }
    if (spatial_dim < 3) {
      std::memset(d_b[2], 0, static_cast<size_t>(n_total) * sizeof(geom_t));
    }
  }

  geom_t local_min[3] = {std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max(),
                         std::numeric_limits<geom_t>::max()};
  geom_t local_max[3] = {std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest(),
                         std::numeric_limits<geom_t>::lowest()};
  if (n_total > 0) {
    for (int d = 0; d < 3; ++d) {
      for (ptrdiff_t i = 0; i < n_total; ++i) {
        local_min[d] = std::min(local_min[d], d_b[d][i]);
        local_max[d] = std::max(local_max[d], d_b[d][i]);
      }
    }
  }
  geom_t global_min[3];
  geom_t global_max[3];
  SMESH_MPI_CATCH(MPI_Allreduce(local_min, global_min, 3, mpi_type<geom_t>(),
                                MPI_MIN, comm));
  SMESH_MPI_CATCH(MPI_Allreduce(local_max, global_max, 3, mpi_type<geom_t>(),
                                MPI_MAX, comm));

  if (n_total == 0) {
    return SMESH_SUCCESS;
  }

  auto encoding = create_host_buffer<u32>(static_cast<size_t>(n_total));
  u32 *d_enc = encoding->data();
  SMESH_CATCH(encode_with_bbox(ordering_type, n_total, d_b[0], d_b[1], d_b[2],
                               global_min[0], global_max[0], global_min[1],
                               global_max[1], global_min[2], global_max[2],
                               d_enc));

  const bool remap_arg = !sidesets.empty();
  const bool remap_reg = !mesh.sidesets().empty() || !mesh.edgesets().empty();

  for (size_t b = 0; b < n_blocks; ++b) {
    const ptrdiff_t n_e = offset[b + 1] - offset[b];
    if (n_e == 0) {
      continue;
    }
    auto gather_buf = create_host_buffer<idx_t>(static_cast<size_t>(n_e));
    idx_t *d_gather = gather_buf->data();
    const u32 *keys = d_enc + offset[b];
    argsort_segment(keys, 0, n_ons[b], d_gather);
    argsort_segment(keys, n_ons[b], n_owned[b], d_gather);
    argsort_segment(keys, n_owned[b], n_e, d_gather);

    SMESH_CATCH(mesh_block_reorder(nxe[b], n_e, elems[b], d_gather, elems[b]));

    auto block = mesh.block(b);
    if (block->element_mapping() && n_owned[b] > 0) {
      permute_large_idx(n_owned[b], block->element_mapping()->data(), d_gather,
                        0, 0);
    }
    if (block->aura_element_mapping() && n_ghosts[b] > 0) {
      permute_large_idx(n_ghosts[b], block->aura_element_mapping()->data(),
                        d_gather, n_owned[b], n_owned[b]);
    }

    if (remap_arg || remap_reg) {
      auto old_to_new =
          create_host_buffer<element_idx_t>(static_cast<size_t>(n_e));
      element_idx_t *d_otn = old_to_new->data();
      for (ptrdiff_t neu = 0; neu < n_e; ++neu) {
        d_otn[d_gather[neu]] = static_cast<element_idx_t>(neu);
      }
      const block_idx_t bid = static_cast<block_idx_t>(b);
      if (remap_arg) {
        if (remap_sidesets(sidesets, bid, d_otn, n_e) != SMESH_SUCCESS) {
          return SMESH_FAILURE;
        }
      }
      if (mesh.remap_registered_sidesets(bid, d_otn, n_e, sidesets) !=
          SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
      if (mesh.remap_registered_edgesets(bid, d_otn, n_e) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
    }
  }

  if (n_blocks > 1 && mesh.distributed()) {
    auto dist = mesh.distributed();
    auto emap = dist->element_mapping();
    auto amap = dist->aura_element_mapping();
    ptrdiff_t wo = 0;
    ptrdiff_t wa = 0;
    for (size_t b = 0; b < n_blocks; ++b) {
      auto block = mesh.block(b);
      if (n_owned[b] > 0 && block->element_mapping() && emap) {
        memcpy(emap->data() + wo, block->element_mapping()->data(),
               static_cast<size_t>(n_owned[b]) * sizeof(large_idx_t));
        wo += n_owned[b];
      }
      if (n_ghosts[b] > 0 && block->aura_element_mapping() && amap) {
        memcpy(amap->data() + wa, block->aura_element_mapping()->data(),
               static_cast<size_t>(n_ghosts[b]) * sizeof(large_idx_t));
        wa += n_ghosts[b];
      }
    }
  }

  return SMESH_SUCCESS;
}
#endif

} // namespace

int SFC::reorder(Mesh &mesh, const std::vector<std::shared_ptr<Sideset>> &sidesets) {
  SMESH_TRACE_SCOPE("SFC::reorder");
#ifdef SMESH_ENABLE_MPI
  if (mesh.is_distributed() && mesh.comm() && mesh.comm()->size() > 1) {
    return reorder_distributed(impl_->ordering_type, mesh, sidesets);
  }
#endif
  if (mesh.n_blocks() > 1) {
    return reorder_multiblock(impl_->ordering_type, mesh, sidesets);
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
  if (spatial_dim < 3) {
    std::memset(d_b[2], 0, static_cast<size_t>(n_elements) * sizeof(geom_t));
  }

  auto encoding = create_host_buffer<u32>(n_elements);
  u32 *d_enc = encoding->data();
  SMESH_CATCH(iter->second(n_elements, d_b[0], d_b[1], d_b[2], d_enc));

  auto idx = create_host_buffer<idx_t>(n_elements);
  idx_t *d_idx = idx->data();
  argsort(n_elements, d_enc, d_idx);

  SMESH_CATCH(mesh_block_reorder(nxe, n_elements, elems, d_idx, elems));

  const bool remap_arg = !sidesets.empty();
  const bool remap_reg = !mesh.sidesets().empty() || !mesh.edgesets().empty();
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
    if (mesh.remap_registered_edgesets(block_id, d_otn, n_elements) != SMESH_SUCCESS) {
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

  if (!mesh.nodesets().empty()) {
    if (mesh.remap_registered_nodesets(d_n2n, n_nodes) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }

  return SMESH_SUCCESS;
}
} // namespace smesh
