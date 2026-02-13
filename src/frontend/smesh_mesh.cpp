#include "smesh_mesh.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_build.hpp"
#include "smesh_conversion.hpp"
#include "smesh_glob.hpp"
#include "smesh_graph.hpp"
#include "smesh_mask.hpp"
#include "smesh_multiblock_graph.hpp"
#include "smesh_path.hpp"
#include "smesh_promotions.hpp"
#include "smesh_read.hpp"
#include "smesh_refine.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_tracer.hpp"
#include "smesh_write.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_read.hpp"
#include "smesh_distributed_write.hpp"
#endif

#include <algorithm>
#include <list>
#include <map>
#include <vector>

namespace smesh {

static ptrdiff_t
max_node_id(const enum ElemType element_type, const ptrdiff_t n_elements,
            const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT elems) {
  const int nxe = elem_num_nodes(element_type);
  idx_t max_id{-1};
  for (int v = 0; v < nxe; ++v) {
    const idx_t *const col = elems[v];
    for (ptrdiff_t e = 0; e < n_elements; ++e) {
      max_id = std::max(max_id, col[e]);
    }
  }
  return static_cast<ptrdiff_t>(max_id);
}

class Mesh::Block::Impl {
public:
  std::string name;
  enum ElemType element_type;
  SharedBuffer<idx_t *> elements;
};

Mesh::Block::Block() : impl_(std::make_unique<Impl>()) {}

Mesh::Block::~Block() = default;

const std::string &Mesh::Block::name() const { return impl_->name; }
enum ElemType Mesh::Block::element_type() const { return impl_->element_type; }
int Mesh::Block::n_nodes_per_element() const {
  return elem_num_nodes(impl_->element_type);
}
const SharedBuffer<idx_t *> &Mesh::Block::elements() const {
  return impl_->elements;
}

void Mesh::Block::set_name(const std::string &name) { impl_->name = name; }
void Mesh::Block::set_element_type(enum ElemType element_type) {
  impl_->element_type = element_type;
}
void Mesh::Block::set_elements(SharedBuffer<idx_t *> elements) {
  impl_->elements = elements;
}

ptrdiff_t Mesh::Block::n_elements() const { return impl_->elements->extent(1); }

class Mesh::Impl {
public:
  std::shared_ptr<Communicator> comm;
  std::vector<std::shared_ptr<Block>> blocks;

  SharedBuffer<geom_t *> points; // Node coordinates

  // MPI-related data using Buffers
  SharedBuffer<idx_t> node_mapping;
  SharedBuffer<int> node_owner;
  SharedBuffer<element_idx_t> element_mapping;
  SharedBuffer<idx_t> node_offsets;
  SharedBuffer<idx_t> ghosts;

  // Metadata
  int spatial_dim;
  ptrdiff_t nnodes;

  // MPI ownership info
  ptrdiff_t n_owned_nodes;
  ptrdiff_t n_owned_nodes_with_ghosts;
  ptrdiff_t n_owned_elements;
  ptrdiff_t n_owned_elements_with_ghosts;
  ptrdiff_t n_shared_elements;

  std::shared_ptr<NodeToNodeGraph> crs_graph;
  std::shared_ptr<NodeToNodeGraph> crs_graph_upper_triangular;
  std::shared_ptr<NodeToElementGraph> node_to_element_graph;

  ~Impl() {}

  void clear() {
    comm = nullptr;
    blocks.clear();
    spatial_dim = 0;
    nnodes = 0;
    points = nullptr;
    node_mapping = nullptr;
    node_owner = nullptr;
    element_mapping = nullptr;
    node_offsets = nullptr;
    ghosts = nullptr;
    n_owned_nodes = 0;
    n_owned_nodes_with_ghosts = 0;
    n_owned_elements = 0;
    n_owned_elements_with_ghosts = 0;
    n_shared_elements = 0;
    crs_graph = nullptr;
    crs_graph_upper_triangular = nullptr;
  }

  // Helper methods for backward compatibility
  ptrdiff_t total_elements() const {
    ptrdiff_t total = 0;
    for (const auto &block : blocks) {
      if (block && block->elements()) {
        // For Buffer<T*>, extent(1) gives the number of elements
        total += block->elements()->extent(1);
      }
    }
    return total;
  }

  enum ElemType default_element_type() const {
    if (blocks.empty() || !blocks[0]) {
      return INVALID;
    }
    return blocks[0]->element_type();
  }

  SharedBuffer<idx_t *> default_elements() const {
    if (blocks.empty() || !blocks[0]) {
      return nullptr;
    }
    return blocks[0]->elements();
  }

  void create_node_to_element_graph() {
    if (node_to_element_graph) {
      return;
    }
    node_to_element_graph = std::make_shared<NodeToElementGraph>();

    count_t *rowptr{nullptr};
    element_idx_t *colidx{nullptr};

    create_n2e(default_elements()->extent(1), nnodes,
               default_elements()->extent(0), default_elements()->data(),
               &rowptr, &colidx);

    node_to_element_graph = std::make_shared<Mesh::NodeToElementGraph>(
        Buffer<count_t>::own(nnodes + 1, rowptr, free, MEMORY_SPACE_HOST),
        Buffer<element_idx_t>::own(rowptr[nnodes], colidx, free,
                                   MEMORY_SPACE_HOST));
  }
};

std::shared_ptr<Communicator> Mesh::comm() const { return impl_->comm; }

Mesh::Mesh(const std::shared_ptr<Communicator> &comm,
           enum ElemType element_type, SharedBuffer<idx_t *> elements,
           SharedBuffer<geom_t *> points)
    : impl_(std::make_unique<Impl>()) {
  impl_->comm = comm;
  impl_->spatial_dim = points->extent(0);
  impl_->nnodes = points->extent(1);
  impl_->points = points;

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(element_type);
  default_block->set_elements(elements);
  impl_->blocks.push_back(default_block);
}

Mesh::Mesh(const std::shared_ptr<Communicator> &comm,
           const std::vector<std::shared_ptr<Block>> &blocks,
           SharedBuffer<geom_t *> points)
    : impl_(std::make_unique<Impl>()) {
  impl_->comm = comm;
  impl_->spatial_dim = points->extent(0);
  impl_->nnodes = points->extent(1);
  impl_->points = points;
  impl_->blocks = blocks;
}

Mesh::Mesh() : impl_(std::make_unique<Impl>()) {
  impl_->clear();
  impl_->comm = Communicator::world();
}

Mesh::Mesh(const std::shared_ptr<Communicator> &comm)
    : impl_(std::make_unique<Impl>()) {
  impl_->clear();
  impl_->comm = comm;
}

Mesh::~Mesh() = default;

// Block-related methods
size_t Mesh::n_blocks() const { return impl_->blocks.size(); }

std::shared_ptr<const Mesh::Block> Mesh::block(size_t index) const {
  if (index >= impl_->blocks.size() || !impl_->blocks[index]) {
    SMESH_ERROR("Block index out of range");
  }
  return impl_->blocks[index];
}

std::shared_ptr<Mesh::Block> Mesh::block(size_t index) {
  if (index >= impl_->blocks.size() || !impl_->blocks[index]) {
    SMESH_ERROR("Block index out of range");
  }
  return impl_->blocks[index];
}

void Mesh::add_block(const std::string &name, enum ElemType element_type,
                     SharedBuffer<idx_t *> elements) {
  auto new_block = std::make_shared<Block>();
  new_block->set_name(name);
  new_block->set_element_type(element_type);
  new_block->set_elements(elements);
  impl_->blocks.push_back(new_block);
}

void Mesh::remove_block(size_t index) {
  if (index >= impl_->blocks.size()) {
    SMESH_ERROR("Block index out of range");
  }

  impl_->blocks.erase(impl_->blocks.begin() + index);
}

int Mesh::read(const Path &path) {
  SMESH_TRACE_SCOPE("Mesh::read");

#ifdef SMESH_ENABLE_MPI
  int comm_size = impl_->comm->size();

  if (comm_size == 1)
#endif
  {
    idx_t **elements = nullptr;
    geom_t **points = nullptr;
    int nnodesxelem;
    int spatial_dim;
    ptrdiff_t nelements;

    if (mesh_from_folder(path, &nnodesxelem, &nelements, &elements,
                         &spatial_dim, &impl_->nnodes,
                         &points) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    auto elements_buffer =
        manage_host_buffer<idx_t>(nnodesxelem, nelements, elements);
    impl_->points =
        manage_host_buffer<geom_t>(spatial_dim, impl_->nnodes, points);
    impl_->spatial_dim = spatial_dim;
    impl_->nnodes = impl_->nnodes;
    impl_->n_owned_nodes = impl_->nnodes;
    impl_->n_owned_elements = nelements;
    impl_->n_owned_elements_with_ghosts = 0;
    impl_->n_shared_elements = 0;
    impl_->n_owned_nodes_with_ghosts = 0;

    // Create default block
    auto default_block = std::make_shared<Block>();
    default_block->set_name("default");
    default_block->set_element_type((enum ElemType)nnodesxelem);
    default_block->set_elements(elements_buffer);
    impl_->blocks.push_back(default_block);
  }
#ifdef SMESH_ENABLE_MPI
  else {
    int nnodesxelem;
    ptrdiff_t nelements;
    idx_t **elements;
    int spatial_dim;
    ptrdiff_t nnodes;
    geom_t **points;
    ptrdiff_t n_owned_nodes;
    ptrdiff_t n_owned_elements;
    element_idx_t *element_mapping;
    idx_t *node_mapping;
    int *node_owner;
    idx_t *node_offsets;
    idx_t *ghosts;
    ptrdiff_t n_owned_nodes_with_ghosts;
    ptrdiff_t n_shared_elements;
    ptrdiff_t n_owned_elements_with_ghosts;

    if (mesh_from_folder(impl_->comm->get(), path, &nnodesxelem, &nelements,
                         &elements, &spatial_dim, &nnodes, &points,
                         &n_owned_nodes, &n_owned_elements, &element_mapping,
                         &node_mapping, &node_owner, &node_offsets, &ghosts,
                         &n_owned_nodes_with_ghosts, &n_shared_elements,
                         &n_owned_elements_with_ghosts) != SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }

    auto elements_buffer =
        manage_host_buffer<idx_t>(nnodesxelem, nelements, elements);
    impl_->points = manage_host_buffer<geom_t>(spatial_dim, nnodes, points);
    impl_->node_mapping = manage_host_buffer<idx_t>(nnodes, node_mapping);
    impl_->node_owner = manage_host_buffer<int>(nnodes, node_owner);
    impl_->element_mapping =
        manage_host_buffer<element_idx_t>(nelements, element_mapping);

    int comm_size;
    MPI_Comm_size(impl_->comm->get(), &comm_size);
    impl_->node_offsets =
        manage_host_buffer<idx_t>(comm_size + 1, node_offsets);

    ptrdiff_t n_ghost_nodes = nnodes - n_owned_nodes;
    impl_->ghosts = manage_host_buffer<idx_t>(n_ghost_nodes, ghosts);

    impl_->n_owned_nodes = n_owned_nodes;
    impl_->n_owned_nodes_with_ghosts = n_owned_nodes_with_ghosts;
    impl_->n_owned_elements = n_owned_elements;
    impl_->n_owned_elements_with_ghosts = n_owned_elements_with_ghosts;
    impl_->n_shared_elements = n_shared_elements;

    // Create default block
    auto default_block = std::make_shared<Block>();
    default_block->set_name("default");
    default_block->set_element_type((enum ElemType)nnodesxelem);
    default_block->set_elements(elements_buffer);
    impl_->blocks.push_back(default_block);
  }
#endif // SMESH_ENABLE_MPI

  int SMESH_USE_MACRO = 0;
  SMESH_READ_ENV(SMESH_USE_MACRO, atoi);

  if (SMESH_USE_MACRO) {
    for (auto &block : blocks()) {
      if (block) {
        block->set_element_type(macro_type_variant(block->element_type()));
      }
    }
  }

  return SMESH_SUCCESS;
}

const std::vector<std::shared_ptr<Mesh::Block>> &Mesh::blocks() const {
  return impl_->blocks;
}

int Mesh::write(const Path &path) const {
  SMESH_TRACE_SCOPE("Mesh::write");

  create_directory(path);

  if (impl_->node_mapping) {
    Path path_node_mapping =
        path / ("node_mapping." + std::string(TypeToString<idx_t>::value()));
    impl_->node_mapping->to_file(path_node_mapping);
  }

  // Write the default block (block 0)
  if (impl_->blocks.empty() || !impl_->blocks[0]) {
    return SMESH_FAILURE;
  }

  if (impl_->blocks.size() == 1) {
    return mesh_to_folder(path, impl_->blocks[0]->element_type(),
                          impl_->blocks[0]->elements()->extent(1),
                          impl_->blocks[0]->elements()->data(),
                          impl_->spatial_dim, impl_->nnodes,
                          impl_->points->data());
  } else {
    std::vector<ptrdiff_t> n_elements;
    std::vector<enum ElemType> element_types;
    std::vector<idx_t **> elements;
    std::vector<std::string> block_names;

    for (auto &block : impl_->blocks) {
      n_elements.push_back(block->elements()->extent(1));
      element_types.push_back(block->element_type());
      elements.push_back(block->elements()->data());
      block_names.push_back(block->name());
    }
    return mesh_multiblock_to_folder(
        path, block_names, element_types, n_elements, elements.data(),
        impl_->spatial_dim, impl_->nnodes, impl_->points->data());
  }
}

const geom_t *Mesh::points(const int coord) const {
  assert(coord < spatial_dimension());
  assert(coord >= 0);
  return impl_->points->data()[coord];
}

const idx_t *Mesh::idx(const int node_num) const {
  assert(node_num < n_nodes_per_element());
  assert(node_num >= 0);
  return impl_->default_elements()->data()[node_num];
}

std::shared_ptr<Mesh::NodeToNodeGraph> Mesh::node_to_node_graph() {
  initialize_node_to_node_graph();
  return impl_->crs_graph;
}

std::shared_ptr<Mesh::NodeToElementGraph> Mesh::node_to_element_graph() {
  impl_->create_node_to_element_graph();
  return impl_->node_to_element_graph;
}

SharedBuffer<element_idx_t> Mesh::half_face_table() {
  // FIXME it should be allocated outisde
  element_idx_t *table{nullptr};
  create_element_adj_table(n_elements(), n_nodes(), element_type(),
                           default_elements()->data(), &table);

  int nsxe = elem_num_sides(element_type());
  return manage_host_buffer<element_idx_t>(n_elements() * nsxe, table);
}

std::shared_ptr<Mesh::NodeToNodeGraph>
Mesh::create_node_to_node_graph(const enum ElemType element_type) {
  if (impl_->default_element_type() == element_type) {
    return node_to_node_graph();
  }

  const ptrdiff_t n_nodes = max_node_id(element_type, impl_->total_elements(),
                                        impl_->default_elements()->data()) +
                            1;

  count_t *rowptr{nullptr};
  idx_t *colidx{nullptr};
  create_crs_graph_for_elem_type(element_type, impl_->total_elements(), n_nodes,
                                 impl_->default_elements()->data(), &rowptr,
                                 &colidx);

  auto crs_graph = std::make_shared<Mesh::NodeToNodeGraph>(
      Buffer<count_t>::own(n_nodes + 1, rowptr, free, MEMORY_SPACE_HOST),
      Buffer<idx_t>::own(rowptr[n_nodes], colidx, free, MEMORY_SPACE_HOST));

  return crs_graph;
}

int Mesh::initialize_node_to_node_graph() {
  if (impl_->crs_graph) {
    return SMESH_SUCCESS;
  }

  SMESH_TRACE_SCOPE("Mesh::initialize_node_to_node_graph");

  impl_->crs_graph = std::make_shared<NodeToNodeGraph>();

  count_t *rowptr{nullptr};
  idx_t *colidx{nullptr};

  if (impl_->blocks.size() == 1) {
    create_crs_graph_for_elem_type(
        impl_->default_element_type(), impl_->total_elements(), impl_->nnodes,
        impl_->default_elements()->data(), &rowptr, &colidx);
  } else {
    // AoS to SoA
    std::vector<enum ElemType> element_types;
    std::vector<ptrdiff_t> n_elements;
    std::vector<idx_t **> elements;

    for (auto &block : impl_->blocks) {
      element_types.push_back(block->element_type());
      n_elements.push_back(block->elements()->extent(1));
      elements.push_back(block->elements()->data());
    }

    create_multiblock_crs_graph(impl_->blocks.size(), element_types.data(),
                                n_elements.data(), elements.data(),
                                impl_->nnodes, &rowptr, &colidx);
  }

  impl_->crs_graph = std::make_shared<Mesh::NodeToNodeGraph>(
      Buffer<count_t>::own(impl_->nnodes + 1, rowptr, free, MEMORY_SPACE_HOST),
      Buffer<idx_t>::own(rowptr[impl_->nnodes], colidx, free,
                         MEMORY_SPACE_HOST));

  return SMESH_SUCCESS;
}

std::shared_ptr<Mesh::NodeToNodeGraph>
Mesh::node_to_node_graph_upper_triangular() {
  if (impl_->crs_graph_upper_triangular)
    return impl_->crs_graph_upper_triangular;
  SMESH_TRACE_SCOPE("Mesh::node_to_node_graph_upper_triangular");

  count_t *rowptr{nullptr};
  idx_t *colidx{nullptr};

  if (impl_->blocks.size() == 1) {
    create_crs_graph_upper_triangular_from_element(
        impl_->total_elements(), impl_->nnodes,
        elem_num_nodes(impl_->default_element_type()),
        impl_->default_elements()->data(), &rowptr, &colidx);
  } else {
    // AoS to SoA
    std::vector<enum ElemType> element_types;
    std::vector<ptrdiff_t> n_elements;
    std::vector<idx_t **> elements;

    for (auto &block : impl_->blocks) {
      element_types.push_back(block->element_type());
      n_elements.push_back(block->elements()->extent(1));
      elements.push_back(block->elements()->data());
    }

    create_multiblock_crs_graph_upper_triangular(
        impl_->blocks.size(), element_types.data(), n_elements.data(),
        elements.data(), impl_->nnodes, &rowptr, &colidx);
  }

  impl_->crs_graph_upper_triangular = std::make_shared<Mesh::NodeToNodeGraph>(
      Buffer<count_t>::own(impl_->nnodes + 1, rowptr, free, MEMORY_SPACE_HOST),
      Buffer<idx_t>::own(rowptr[impl_->nnodes], colidx, free,
                         MEMORY_SPACE_HOST));

  return impl_->crs_graph_upper_triangular;
}

int Mesh::convert_to_macro_element_mesh() {
  if (!impl_->blocks.empty() && impl_->blocks[0]) {
    impl_->blocks[0]->set_element_type(
        macro_type_variant(impl_->blocks[0]->element_type()));
  }
  return SMESH_SUCCESS;
}

SharedBuffer<count_t> Mesh::node_to_node_rowptr() const {
  return impl_->crs_graph->rowptr();
}
SharedBuffer<idx_t> Mesh::node_to_node_colidx() const {
  return impl_->crs_graph->colidx();
}

std::shared_ptr<Mesh>
Mesh::create_hex8_cube(const std::shared_ptr<Communicator> &comm, const int nx,
                       const int ny, const int nz, const geom_t xmin,
                       const geom_t ymin, const geom_t zmin, const geom_t xmax,
                       const geom_t ymax, const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);

  ret->impl_->spatial_dim = 3;
  ret->impl_->nnodes = nnodes;
  ret->impl_->points = create_host_buffer<geom_t>(3, nnodes);
  auto elements_buffer = create_host_buffer<idx_t>(8, nelements);

  ret->impl_->n_owned_nodes = nnodes;
  ret->impl_->n_owned_elements = nelements;

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_hex8_cube<idx_t, geom_t>(nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                                     zmax, elements, points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(HEX8);
  default_block->set_elements(elements_buffer);
  ret->impl_->blocks.push_back(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_tri3_square(const std::shared_ptr<Communicator> &comm,
                         const int nx, const int ny, const geom_t xmin,
                         const geom_t ymin, const geom_t xmax,
                         const geom_t ymax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = 2 * nx * ny;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1);

  ret->impl_->spatial_dim = 2;
  ret->impl_->nnodes = nnodes;
  ret->impl_->points = create_host_buffer<geom_t>(2, nnodes);
  auto elements_buffer = create_host_buffer<idx_t>(3, nelements);

  ret->impl_->n_owned_nodes = nnodes;
  ret->impl_->n_owned_elements = nelements;

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_tri3_square<idx_t, geom_t>(nx, ny, xmin, ymin, xmax, ymax, elements,
                                       points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(TRI3);
  default_block->set_elements(elements_buffer);
  ret->impl_->blocks.push_back(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_quad4_square(const std::shared_ptr<Communicator> &comm,
                          const int nx, const int ny, const geom_t xmin,
                          const geom_t ymin, const geom_t xmax,
                          const geom_t ymax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1);

  ret->impl_->spatial_dim = 2;
  ret->impl_->nnodes = nnodes;
  ret->impl_->points = create_host_buffer<geom_t>(2, nnodes);
  auto elements_buffer = create_host_buffer<idx_t>(4, nelements);

  ret->impl_->n_owned_nodes = nnodes;
  ret->impl_->n_owned_elements = nelements;

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_quad4_square<idx_t, geom_t>(nx, ny, xmin, ymin, xmax, ymax,
                                        elements, points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(QUAD4);
  default_block->set_elements(elements_buffer);
  ret->impl_->blocks.push_back(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_square(const std::shared_ptr<Communicator> &comm,
                    const enum ElemType element_type, const int nx,
                    const int ny, const geom_t xmin, const geom_t ymin,
                    const geom_t xmax, const geom_t ymax) {
  switch (element_type) {
  case QUAD4:
    return create_quad4_square(comm, nx, ny, xmin, ymin, xmax, ymax);
  case TRI3:
    return create_tri3_square(comm, nx, ny, xmin, ymin, xmax, ymax);
  default:
    SMESH_ERROR("Invalid element type: %d\n", element_type);
    return nullptr;
  }
}

std::shared_ptr<Mesh>
Mesh::create_quad4_ring(const std::shared_ptr<Communicator> &comm,
                        const geom_t inner_radius, const geom_t outer_radius,
                        const int nlayers, const int nelements) {
  auto elements = create_host_buffer<idx_t>(4, nlayers * nelements);
  auto points = create_host_buffer<geom_t>(3, (nlayers + 1) * nelements);
  mesh_fill_quad4_ring<idx_t, geom_t>(inner_radius, outer_radius, nlayers,
                                      nelements, elements->data(),
                                      points->data());

  auto ret = std::make_shared<Mesh>(comm, QUAD4, elements, points);
  return ret;
}

std::shared_ptr<Mesh> Mesh::create_hex8_checkerboard_cube(
    const std::shared_ptr<Communicator> &comm, const int nx, const int ny,
    const int nz, const geom_t xmin, const geom_t ymin, const geom_t zmin,
    const geom_t xmax, const geom_t ymax, const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);

  if (nx % 2 != 0 || ny % 2 != 0 || nz % 2 != 0) {
    SMESH_ERROR("nx, ny, and nz must be even");
  }

  ret->impl_->spatial_dim = 3;
  ret->impl_->nnodes = nnodes;
  ret->impl_->points = create_host_buffer<geom_t>(3, nnodes);
  auto white_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);
  auto black_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);

  ret->impl_->n_owned_nodes = nnodes;
  ret->impl_->n_owned_elements = nelements;

  auto points = ret->impl_->points->data();
  auto white_elements = white_elements_buffer->data();
  auto black_elements = black_elements_buffer->data();

  mesh_fill_hex8_checkerboard_cube<idx_t, geom_t>(
      nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax, white_elements,
      black_elements, points);
  // Create white and black blocks
  auto white_block = std::make_shared<Block>();
  white_block->set_name("white");
  white_block->set_element_type(HEX8);
  white_block->set_elements(white_elements_buffer);
  ret->impl_->blocks.push_back(white_block);

  auto black_block = std::make_shared<Block>();
  black_block->set_name("black");
  black_block->set_element_type(HEX8);
  black_block->set_elements(black_elements_buffer);
  ret->impl_->blocks.push_back(black_block);
  return ret;
}

std::shared_ptr<Mesh> Mesh::create_hex8_bidomain_cube(
    const std::shared_ptr<Communicator> &comm, const int nx, const int ny,
    const int nz, const geom_t xmin, const geom_t ymin, const geom_t zmin,
    const geom_t xmax, const geom_t ymax, const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);

  ret->impl_->spatial_dim = 3;
  ret->impl_->nnodes = nnodes;
  ret->impl_->points = create_host_buffer<geom_t>(3, nnodes);
  auto left_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);
  auto right_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);

  ret->impl_->n_owned_nodes = nnodes;
  ret->impl_->n_owned_elements = nelements;

  auto points = ret->impl_->points->data();
  auto left_elements = left_elements_buffer->data();
  auto right_elements = right_elements_buffer->data();

  mesh_fill_hex8_bidomain_cube<idx_t, geom_t>(
      nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax, 0, nx / 2, left_elements,
      right_elements, points);

  // Create left and right blocks
  auto left_block = std::make_shared<Block>();
  left_block->set_name("left");
  left_block->set_element_type(HEX8);
  left_block->set_elements(left_elements_buffer);
  ret->impl_->blocks.push_back(left_block);

  auto right_block = std::make_shared<Block>();
  right_block->set_name("right");
  right_block->set_element_type(HEX8);
  right_block->set_elements(right_elements_buffer);
  ret->impl_->blocks.push_back(right_block);

  return ret;
}

int Mesh::spatial_dimension() const { return impl_->spatial_dim; }
int Mesh::n_nodes_per_element() const {
  return elem_num_nodes(impl_->default_element_type());
}

ptrdiff_t Mesh::n_nodes() const { return impl_->nnodes; }
ptrdiff_t Mesh::n_elements() const { return impl_->total_elements(); }

enum ElemType Mesh::element_type() const {
  return impl_->default_element_type();
}

ptrdiff_t Mesh::n_owned_nodes() const { return impl_->n_owned_nodes; }
ptrdiff_t Mesh::n_owned_nodes_with_ghosts() const {
  return impl_->n_owned_nodes_with_ghosts;
}
ptrdiff_t Mesh::n_owned_elements() const { return impl_->n_owned_elements; }
ptrdiff_t Mesh::n_owned_elements_with_ghosts() const {
  return impl_->n_owned_elements_with_ghosts;
}
ptrdiff_t Mesh::n_shared_elements() const { return impl_->n_shared_elements; }

SharedBuffer<idx_t> Mesh::node_mapping() const { return impl_->node_mapping; }
SharedBuffer<idx_t> Mesh::element_mapping() const {
  return impl_->element_mapping;
}

SharedBuffer<idx_t> Mesh::node_offsets() const { return impl_->node_offsets; }
SharedBuffer<idx_t> Mesh::ghosts() const { return impl_->ghosts; }
SharedBuffer<int> Mesh::node_owner() const { return impl_->node_owner; }

SharedBuffer<geom_t *> Mesh::points() { return impl_->points; }

SharedBuffer<idx_t *> Mesh::elements() { return impl_->default_elements(); }
SharedBuffer<idx_t *> Mesh::default_elements() {
  return impl_->default_elements();
}

void Mesh::set_node_mapping(const SharedBuffer<idx_t> &node_mapping) {
  impl_->node_mapping = node_mapping;
}

void Mesh::set_comm(const std::shared_ptr<Communicator> &comm) {
  impl_->comm = comm;
}

void Mesh::set_element_type(const enum ElemType element_type) {
  if (!impl_->blocks.empty() && impl_->blocks[0]) {
    impl_->blocks[0]->set_element_type(element_type);
  }
}

std::vector<std::shared_ptr<Mesh::Block>>
Mesh::blocks(const std::vector<std::string> &block_names) const {
  if (block_names.empty()) {
    return impl_->blocks;
  }

  std::vector<std::shared_ptr<Mesh::Block>> ret;
  for (auto &block : impl_->blocks) {
    if (std::find(block_names.begin(), block_names.end(), block->name()) !=
        block_names.end()) {
      ret.push_back(block);
    }
  }

  return ret;
}

std::shared_ptr<Mesh> Mesh::create_hex8_reference_cube() {
  auto ret = std::make_shared<Mesh>(Communicator::null());
  ret->impl_->spatial_dim = 3;
  ret->impl_->nnodes = 8;
  ret->impl_->points = create_host_buffer<geom_t>(3, 8);
  auto elements_buffer = create_host_buffer<idx_t>(8, 1);

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_hex8_reference_cube<idx_t, geom_t>(elements, points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(HEX8);
  default_block->set_elements(elements_buffer);
  ret->impl_->blocks.push_back(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_tet4_cube(const std::shared_ptr<Communicator> &comm, const int nx,
                       const int ny, const int nz, const geom_t xmin,
                       const geom_t ymin, const geom_t zmin, const geom_t xmax,
                       const geom_t ymax, const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes_vertices = (nx + 1) * (ny + 1) * (nz + 1);
  const ptrdiff_t nnodes_total = nnodes_vertices + nelements;

  ret->impl_->spatial_dim = 3;
  ret->impl_->nnodes = nnodes_total;
  ret->impl_->points = create_host_buffer<geom_t>(3, nnodes_total);
  auto elements_buffer = create_host_buffer<idx_t>(4, nelements * 12);

  ret->impl_->n_owned_nodes = nnodes_total;
  ret->impl_->n_owned_elements = nelements * 12;

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_tet4_cube<idx_t, geom_t>(nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                                     zmax, elements, points);

  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(TET4);
  default_block->set_elements(elements_buffer);
  ret->impl_->blocks.push_back(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_cube(const std::shared_ptr<Communicator> &comm,
                  const enum ElemType element_type, const int nx, const int ny,
                  const int nz, const geom_t xmin, const geom_t ymin,
                  const geom_t zmin, const geom_t xmax, const geom_t ymax,
                  const geom_t zmax) {
  switch (element_type) {
  case HEX8:
    return create_hex8_cube(comm, nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                            zmax);
  case TET4:
    return create_tet4_cube(comm, nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                            zmax);
  default:
    SMESH_ERROR("Invalid element type: %d\n", element_type);
    return nullptr;
  }
}

std::pair<SharedBuffer<geom_t>, SharedBuffer<geom_t>>
Mesh::compute_bounding_box() {
  auto points = impl_->points->data();

  int dim = spatial_dimension();
  auto min = create_host_buffer<geom_t>(dim);
  auto max = create_host_buffer<geom_t>(dim);

  auto d_min = min->data();
  auto d_max = max->data();

  for (int d = 0; d < dim; d++) {
    d_min[d] = points[d][0];
    d_max[d] = points[d][0];
  }

  ptrdiff_t n_nodes = this->n_nodes();

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    for (int d = 0; d < dim; d++) {
      d_min[d] = std::min(d_min[d], points[d][i]);
      d_max[d] = std::max(d_max[d], points[d][i]);
    }
  }

  return {min, max};
}

std::shared_ptr<Mesh::Block> Mesh::find_block(const std::string &name) const {
  for (auto &block : impl_->blocks) {
    if (block->name() == name) {
      return block;
    }
  }
  return nullptr;
}

int Mesh::split_block(const SharedBuffer<element_idx_t> &elements,
                      const std::string &name) {
  if (n_blocks() != 1) {
    SMESH_ERROR("Mesh must have exactly one block to split boundary layer!\n");
    return SMESH_FAILURE;
  }

  {
    const int nxe = n_nodes_per_element();
    const ptrdiff_t n_elements = this->n_elements();

    auto bdry_mask = create_host_buffer<mask_t>(mask_count(n_elements));

    auto d_parent = elements->data();
    auto d_bdry_mask = bdry_mask->data();
    const ptrdiff_t size_sideset = elements->size();

    auto default_block = impl_->blocks[0];

    auto d_elements = default_block->elements()->data();

    ptrdiff_t n_bdry_elements = 0;
    for (ptrdiff_t i = 0; i < size_sideset; i++) {
      if (mask_get(d_parent[i], d_bdry_mask) == 0) {
        n_bdry_elements++;
        mask_set(d_parent[i], d_bdry_mask);
      }
    }

    memset(d_bdry_mask, 0, mask_count(n_elements) * sizeof(mask_t));

    auto bdry_elements = create_host_buffer<idx_t>(nxe, n_bdry_elements);
    ptrdiff_t n_bdry_elements_count = 0;

    auto d_bdry_elements = bdry_elements->data();

    for (int e = 0; e < size_sideset; e++) {
      if (mask_get(d_parent[e], d_bdry_mask) == 0) {
        for (int v = 0; v < nxe; v++) {
          d_bdry_elements[v][n_bdry_elements_count] =
              d_elements[v][d_parent[e]];
        }
        mask_set(d_parent[e], d_bdry_mask);
        n_bdry_elements_count++;
      }
    }

    auto interior_elements =
        create_host_buffer<idx_t>(nxe, n_elements - n_bdry_elements);
    ptrdiff_t n_interior_elements_count = 0;
    auto d_interior_elements = interior_elements->data();

    for (ptrdiff_t i = 0; i < n_elements; i++) {
      if (mask_get(i, d_bdry_mask) == 0) {
        for (int v = 0; v < nxe; v++) {
          assert(n_interior_elements_count <
                 static_cast<ptrdiff_t>(interior_elements->extent(1)));
          d_interior_elements[v][n_interior_elements_count] = d_elements[v][i];
        }
        n_interior_elements_count++;
      }
    }

    // !!!!
    remove_block(0);

    { // Boundary block
      auto block = std::make_shared<Block>();
      block->set_name(name);
      block->set_element_type(default_block->element_type());
      block->set_elements(bdry_elements);
      impl_->blocks.push_back(block);
    }

    { // Interior block
      auto block = std::make_shared<Block>();
      block->set_name(default_block->name());
      block->set_element_type(default_block->element_type());
      block->set_elements(interior_elements);
      impl_->blocks.push_back(block);
    }
  }

  return SMESH_SUCCESS;
}

int Mesh::split_boundary_layer() {
  if (n_blocks() != 1) {
    SMESH_ERROR("Mesh must have exactly one block to split boundary layer!\n");
    return SMESH_FAILURE;
  }

  std::shared_ptr<Sideset> sideset;

  {
    ptrdiff_t n_surf_elements = 0;
    element_idx_t *parent = 0;
    int16_t *side_idx = 0;

    if (extract_skin_sideset(this->n_elements(), this->n_nodes(),
                             this->element_type(), this->elements()->data(),
                             &n_surf_elements, &parent,
                             &side_idx) != SMESH_SUCCESS) {
      SMESH_ERROR("Failed to extract skin!\n");
    }

    sideset = std::make_shared<Sideset>(
        this->comm(), manage_host_buffer(n_surf_elements, parent),
        manage_host_buffer(n_surf_elements, side_idx));
  }

  return split_block(sideset->parent(), "boundary_layer");
}

int Mesh::renumber_nodes() {
  const int nxe = n_nodes_per_element();
  auto n_nodes = this->n_nodes();

  auto new_idx_buff = create_host_buffer<idx_t>(n_nodes);

  auto new_idx = new_idx_buff->data();
  for (ptrdiff_t i = 0; i < n_nodes; i++) {
    new_idx[i] = -1;
  }

  idx_t next_node_id = 0;
  for (auto &b : impl_->blocks) {
    auto elements = b->elements()->data();
    auto n_elements = b->n_elements();

    for (ptrdiff_t e = 0; e < n_elements; e++) {
      for (int v = 0; v < nxe; v++) {
        auto node = elements[v][e];
        if (new_idx[node] == -1) {
          new_idx[node] = next_node_id++;
        }
      }
    }
  }

  return renumber_nodes(new_idx_buff);

  // const int dim = spatial_dimension();

  // auto new_points_buff = create_host_buffer<geom_t>(dim, n_nodes);
  // auto new_points      = new_points_buff->data();

  // for (int d = 0; d < dim; d++) {
  //     for (ptrdiff_t i = 0; i < n_nodes; i++) {
  //         new_points[d][new_idx[i]] = points[d][i];
  //     }
  // }

  // impl_->points = new_points_buff;

  // for (auto &b : impl_->blocks) {
  //     auto elements   = b->elements()->data();
  //     auto n_elements = b->n_elements();

  //     for (ptrdiff_t e = 0; e < n_elements; e++) {
  //         for (int v = 0; v < nxe; v++) {
  //             elements[v][e] = new_idx[elements[v][e]];
  //         }
  //     }
  // }

  // return SMESH_SUCCESS;
}

int Mesh::renumber_nodes(const SharedBuffer<idx_t> &node_mapping) {
  const int dim = spatial_dimension();
  const int nxe = n_nodes_per_element();
  const ptrdiff_t n_nodes = this->n_nodes();

  auto points = this->points()->data();
  auto new_points_buff = create_host_buffer<geom_t>(dim, n_nodes);
  auto new_points = new_points_buff->data();

  auto d_node_mapping = node_mapping->data();

  for (int d = 0; d < dim; d++) {
    for (ptrdiff_t i = 0; i < n_nodes; i++) {
      assert(d_node_mapping[i] < n_nodes);
      assert(d_node_mapping[i] >= 0);
      new_points[d][d_node_mapping[i]] = points[d][i];
    }
  }

  impl_->points = new_points_buff;

  for (auto &b : impl_->blocks) {
    auto elements = b->elements()->data();
    auto n_elements = b->n_elements();

    for (ptrdiff_t e = 0; e < n_elements; e++) {
      for (int v = 0; v < nxe; v++) {
        elements[v][e] = d_node_mapping[elements[v][e]];
      }
    }
  }

  return SMESH_SUCCESS;
}

std::vector<std::pair<block_idx_t, SharedBuffer<element_idx_t>>>
Mesh::select_elements(const std::function<bool(const geom_t, const geom_t,
                                               const geom_t)> &selector,
                      const std::vector<std::string> &block_names) {
  SMESH_TRACE_SCOPE("Sideset::create_from_selector");

  // const ptrdiff_t nelements = mesh->n_elements();
  const int dim = spatial_dimension();

  auto points = this->points()->data();
  int nxe = n_nodes_per_element();

  size_t n_blocks = this->n_blocks();
  std::vector<std::pair<block_idx_t, SharedBuffer<element_idx_t>>>
      selected_elements;

  for (size_t b = 0; b < n_blocks; b++) {
    auto block = this->block(b);
    if (!block_names.empty() && //
        std::find(block_names.begin(), block_names.end(), block->name()) ==
            block_names.end()) {
      continue;
    }

    const ptrdiff_t nelements = block->n_elements();
    auto elements = block->elements()->data();

    std::list<element_idx_t> selected_element_list;
    for (ptrdiff_t e = 0; e < nelements; e++) {
      // Barycenter of element
      double p[3] = {0, 0, 0};

      for (int v = 0; v < nxe; v++) {
        const idx_t node = elements[v][e];

        for (int d = 0; d < dim; d++) {
          p[d] += points[d][node];
        }
      }

      for (int d = 0; d < dim; d++) {
        p[d] /= nxe;
      }

      if (selector(p[0], p[1], p[2])) {
        selected_element_list.push_back(e);
      }
    }

    const ptrdiff_t nselected_elements = selected_element_list.size();
    auto selected_element =
        create_host_buffer<element_idx_t>(nselected_elements);

    {
      ptrdiff_t idx = 0;
      for (auto p : selected_element_list) {
        selected_element->data()[idx++] = p;
      }
    }

    selected_elements.push_back(std::make_pair(b, selected_element));
  }

  return selected_elements;
}

void Mesh::reorder_elements_from_tags(const SharedBuffer<idx_t> &tags) {
  const ptrdiff_t nelems = n_elements();
  auto temp = create_host_buffer<idx_t>(nelems);
  auto d_temp = temp->data();
  auto d_tags = tags->data();

  auto d_elements = elements()->data();

  idx_t ntags = 0;
  for (ptrdiff_t i = 0; i < nelems; i++) {
    ntags = std::max(ntags, d_tags[i]);
  }

  if (!ntags)
    return;

  ntags += 1;

  auto bookkeeping = create_host_buffer<ptrdiff_t>(ntags);
  auto d_bk = bookkeeping->data();

  int nxe = n_nodes_per_element();
  for (int d = 0; d < nxe; d++) {
    memcpy(d_temp, d_elements[d], nelems * sizeof(idx_t));

    for (ptrdiff_t i = 0; i < nelems; i++) {
      auto t = d_tags[i];
      d_elements[d][d_bk[t]++] = d_temp[i];
    }
  }
}

std::shared_ptr<Mesh> Mesh::clone() const {
  auto ret = std::make_shared<Mesh>();
  SMESH_IMPLEMENT_ME();
  return ret;
}

std::shared_ptr<Mesh> convert_to(const enum ElemType element_type,
                                 const std::shared_ptr<Mesh> &mesh) {

  // FIXME the multiblock case is not really supported yet
  if (mesh->n_blocks() > 1) {
    SMESH_ERROR(
        "Conversion from %d to %d is not supported for multiblock meshes\n",
        mesh->element_type(), element_type);
    return nullptr;
  }

  std::map<std::pair<enum ElemType, enum ElemType>,
           std::function<void(Mesh::Block &, Mesh::Block &)>>
      cmap;

  cmap[std::make_pair(HEX8, TET4)] = [](Mesh::Block &block,
                                        Mesh::Block &new_block) {
    new_block.set_element_type(TET4);
    new_block.set_elements(
        create_host_buffer<idx_t>(4, block.n_elements() * 6));
    mesh_hex8_to_6x_tet4(block.n_elements(), block.elements()->data(),
                         new_block.elements()->data());
  };

  cmap[std::make_pair(TET15, HEX8)] = [](Mesh::Block &block,
                                         Mesh::Block &new_block) {
    new_block.set_element_type(HEX8);
    new_block.set_elements(
        create_host_buffer<idx_t>(8, block.n_elements() * 4));
    mesh_tet15_to_4x_hex8(block.n_elements(), block.elements()->data(),
                          new_block.elements()->data());
  };

  cmap[std::make_pair(WEDGE6, TET4)] = [](Mesh::Block &block,
                                          Mesh::Block &new_block) {
    new_block.set_element_type(TET4);
    new_block.set_elements(
        create_host_buffer<idx_t>(4, block.n_elements() * 3));
    mesh_wedge6_to_3x_tet4(block.n_elements(), block.elements()->data(),
                           new_block.elements()->data());
  };

  std::vector<std::shared_ptr<Mesh::Block>> blocks;
  for (auto &block : mesh->blocks()) {
    auto new_block = std::make_shared<Mesh::Block>();
    new_block->set_name(block->name());
    new_block->set_element_type(element_type);

    if (block->element_type() == element_type) {
      new_block->set_elements(block->elements());
    } else {
      auto it = cmap.find(std::make_pair(block->element_type(), element_type));
      if (it != cmap.end()) {
        it->second(*block, *new_block);
      } else {
        SMESH_ERROR("Conversion from %d to %d is not supported\n",
                    block->element_type(), element_type);
        return nullptr;
      }
    }

    blocks.push_back(new_block);
  }

  return std::make_shared<Mesh>(mesh->comm(), blocks, mesh->points());
}

std::shared_ptr<Mesh> promote_to(const enum ElemType element_type,
                                 const std::shared_ptr<Mesh> &mesh) {
  // FIXME the multiblock case is not really supported yet
  if (mesh->n_blocks() > 1) {
    SMESH_ERROR(
        "Promotion from %d to %d is not supported for multiblock meshes\n",
        mesh->element_type(), element_type);
    return nullptr;
  }

  std::map<std::pair<enum ElemType, enum ElemType>,
           std::function<std::shared_ptr<Mesh>(Mesh &)>>
      cmap;

  cmap[std::make_pair(TET4, TET15)] = [](Mesh &mesh) -> std::shared_ptr<Mesh> {
    auto elements = create_host_buffer<idx_t>(15, mesh.n_elements());
    auto n2n_upper_triangular = mesh.node_to_node_graph_upper_triangular();
    auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
    auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

    auto hft = mesh.half_face_table();
    auto e2e_table = hft->data();

    ptrdiff_t n_new_nodes = 0;
    mesh_tet4_to_tet15(mesh.n_elements(), mesh.n_nodes(),
                       mesh.elements()->data(), n2n_upper_triangular_ptr,
                       n2n_upper_triangular_idx, e2e_table, elements->data(),
                       &n_new_nodes);

    auto points = create_host_buffer<geom_t>(3, n_new_nodes);
    mesh_tet4_to_tet15_points(mesh.n_elements(), mesh.n_nodes(),
                              mesh.points()->data(), n2n_upper_triangular_ptr,
                              n2n_upper_triangular_idx, elements->data(),
                              points->data());

    return std::make_shared<Mesh>(mesh.comm(), TET15, elements, points);
  };

  cmap[std::make_pair(TET4, TET10)] = [](Mesh &mesh) {
    auto n2n_upper_triangular = mesh.node_to_node_graph_upper_triangular();
    auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
    auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

    auto elements = create_host_buffer<idx_t>(10, mesh.n_elements());
    auto points = create_host_buffer<geom_t>(
        mesh.spatial_dimension(),
        n2n_upper_triangular->colidx()->size() + mesh.n_nodes());

    p1_to_p2(TET4, mesh.n_elements(), mesh.elements()->data(),
             mesh.spatial_dimension(), mesh.n_nodes(), mesh.points()->data(),
             n2n_upper_triangular_ptr, n2n_upper_triangular_idx,
             elements->data(), points->data());
    return std::make_shared<Mesh>(mesh.comm(), TET10, elements, points);
  };

  cmap[std::make_pair(TRI3, TRI6)] = [](Mesh &mesh) {
    auto n2n_upper_triangular = mesh.node_to_node_graph_upper_triangular();
    auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
    auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

    auto elements = create_host_buffer<idx_t>(6, mesh.n_elements());
    auto points = create_host_buffer<geom_t>(
        mesh.spatial_dimension(), n2n_upper_triangular->colidx()->size());

    p1_to_p2(TRI3, mesh.n_elements(), mesh.elements()->data(),
             mesh.spatial_dimension(), mesh.n_nodes(), mesh.points()->data(),
             n2n_upper_triangular_ptr, n2n_upper_triangular_idx,
             elements->data(), points->data());
    return std::make_shared<Mesh>(mesh.comm(), TRI6, elements, points);
  };

  auto it = cmap.find(std::make_pair(mesh->element_type(), element_type));
  if (it != cmap.end()) {
    return it->second(*mesh);
  } else {
    SMESH_ERROR("Promotion from %d to %d is not supported\n",
                mesh->element_type(), element_type);
    return nullptr;
  }
}

std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                             const int levels) {
  const int refine_factor = [](const ElemType element_type) {
    switch (element_type) {
    case HEX8:
      return 8;
    case TET4:
      return 8;
    case TRI3:
      return 4;
    default:
      SMESH_ERROR("Refinement factor not supported for element type %d\n",
                  element_type);
      return 0;
    }
  }(mesh->element_type());

  auto out = mesh;
  if (mesh->element_type() == HEX8) {
    const ptrdiff_t n_elements = mesh->n_elements();

    const int ss_levels = pow(2, levels);
    const int nxe = sshex8_nxe(ss_levels);
    const int txe = sshex8_txe(ss_levels);

    auto sshex8_elements = create_host_buffer<idx_t>(nxe, mesh->n_elements());
    auto d_sshex8_elements = sshex8_elements->data();

    ptrdiff_t n_unique_nodes = 0;
    ptrdiff_t interior_start = 0;

    sshex8_generate_elements(ss_levels, n_elements, mesh->n_nodes(),
                             mesh->elements()->data(), d_sshex8_elements,
                             &n_unique_nodes, &interior_start);

    ptrdiff_t n_micro_elements = n_elements * txe;
    auto hex8_elements = create_host_buffer<idx_t>(8, n_micro_elements);
    auto d_hex8_elements = hex8_elements->data();

    sshex8_to_standard_hex8_mesh(ss_levels, n_elements, d_sshex8_elements,
                                 d_hex8_elements);

    auto hex8_points = create_host_buffer<geom_t>(3, n_unique_nodes);
    auto d_hex8_points = hex8_points->data();

    sshex8_fill_points(ss_levels, n_elements, d_sshex8_elements,
                       mesh->points()->data(), d_hex8_points);
    out =
        std::make_shared<Mesh>(mesh->comm(), HEX8, hex8_elements, hex8_points);

  } else {
    for (int i = 0; i < levels; i++) {
      auto n2n_upper_triangular = out->node_to_node_graph_upper_triangular();
      auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
      auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

      auto refined_elements = create_host_buffer<idx_t>(
          out->n_nodes_per_element(), out->n_elements() * refine_factor);

      auto refined_points = create_host_buffer<geom_t>(
          out->spatial_dimension(),
          n2n_upper_triangular->colidx()->size() + out->n_nodes());

      int err = mesh_refine(out->element_type(), out->n_elements(),
                            out->elements()->data(), out->spatial_dimension(),
                            out->n_nodes(), out->points()->data(),
                            n2n_upper_triangular_ptr, n2n_upper_triangular_idx,
                            refined_elements->data(), refined_points->data());

      if (err != SMESH_SUCCESS) {
        SMESH_ERROR("Refinement failed\n");
        return nullptr;
      }

      out = std::make_shared<Mesh>(out->comm(), out->element_type(),
                                   refined_elements, refined_points);
    }
  }

  return out;
}

std::shared_ptr<Sideset> skin_sideset(const std::shared_ptr<Mesh> &mesh) {
  auto hft = mesh->half_face_table();
  auto e2e_table = hft->data();

  ptrdiff_t n_surf_elements = 0;
  element_idx_t *parent_element = 0;
  i16 *side_idx = 0;

  int err = extract_sideset_from_adj_table(
      mesh->element_type(), mesh->n_elements(), e2e_table, &n_surf_elements,
      &parent_element, &side_idx);

  if (err != SMESH_SUCCESS) {
    SMESH_ERROR("Unable to extract skin sideset!\n");
    return nullptr;
  }

  return std::make_shared<Sideset>(
      mesh->comm(),
      manage_host_buffer<element_idx_t>(n_surf_elements, parent_element),
      manage_host_buffer<i16>(n_surf_elements, side_idx), 0);
}

std::shared_ptr<Mesh>
mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                  const std::shared_ptr<Sideset> &sideset) {

  auto [surface_type, surface_elements] =
      create_surface_from_sideset(mesh, sideset);

  const ptrdiff_t n_nodes = mesh->n_nodes();
  auto vol2surf = create_host_buffer<idx_t>(n_nodes);
  auto b_vol2surf = vol2surf->data();
  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    b_vol2surf[i] = invalid_idx<idx_t>();
  }

  const int nnxs = surface_elements->extent(0);
  ptrdiff_t n_surf_elements = surface_elements->extent(1);
  auto b_surface_elements = surface_elements->data();

  ptrdiff_t n_surf_nodes = 0;
  for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
    for (int d = 0; d < nnxs; ++d) {
      idx_t idx = b_surface_elements[d][i];
      if (b_vol2surf[idx] == invalid_idx<idx_t>()) {
        b_vol2surf[idx] = n_surf_nodes++;
      }
    }
  }

  auto b_points = mesh->points()->data();
  auto surf_points =
      create_host_buffer<geom_t>(mesh->spatial_dimension(), n_surf_nodes);

  auto mapping = create_host_buffer<idx_t>(n_surf_nodes);
  auto b_surf_points = surf_points->data();
  auto b_mapping = mapping->data();

  int spatial_dim = mesh->spatial_dimension();
  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    if (b_vol2surf[i] == invalid_idx<idx_t>())
      continue;

    b_mapping[b_vol2surf[i]] = i;
    for (int d = 0; d < spatial_dim; ++d) {
      b_surf_points[d][b_vol2surf[i]] = b_points[d][i];
    }
  }

  auto ret = std::make_shared<Mesh>(mesh->comm(), surface_type,
                                    surface_elements, surf_points);
  ret->set_node_mapping(mapping);
  return ret;
}

std::shared_ptr<Mesh> skin(const std::shared_ptr<Mesh> &mesh) {
  auto sideset = skin_sideset(mesh);
  return mesh_from_sideset(mesh, sideset);
}

std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                              const geom_t height, const ptrdiff_t nlayers) {

  // This is a hack
  if (mesh->n_nodes_per_element() == 4) {
    auto hex8_elements =
        create_host_buffer<idx_t>(8, mesh->n_elements() * nlayers);

    auto hex8_points =
        create_host_buffer<geom_t>(3, mesh->n_nodes() * (nlayers + 1));

    quad4_to_hex8_extrude(mesh->n_elements(), mesh->n_nodes(),
                          mesh->elements()->data(), mesh->points()->data(),
                          nlayers, height, hex8_elements->data(),
                          hex8_points->data());

    return std::make_shared<Mesh>(mesh->comm(), HEX8, hex8_elements,
                                  hex8_points);
  } else if (mesh->element_type() == TRI3) {
    auto wedge6_elements =
        create_host_buffer<idx_t>(6, mesh->n_elements() * nlayers);

    auto wedge6_points =
        create_host_buffer<geom_t>(3, mesh->n_nodes() * (nlayers + 1));

    tri3_to_wedge6_extrude(mesh->n_elements(), mesh->n_nodes(),
                           mesh->elements()->data(), mesh->points()->data(),
                           nlayers, height, wedge6_elements->data(),
                           wedge6_points->data());

    return std::make_shared<Mesh>(mesh->comm(), WEDGE6, wedge6_elements,
                                  wedge6_points);
  } else {
    SMESH_ERROR("Extrusion not supported for element type %d\n",
                mesh->element_type());
    return nullptr;
  }
}
} // namespace smesh
