#include "smesh_mesh.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_build.hpp"
#include "smesh_conversion.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_glob.hpp"
#include "smesh_graph.hpp"
#include "smesh_mask.hpp"
#include "smesh_multiblock_graph.hpp"
#include "smesh_path.hpp"
#include "smesh_promotions.hpp"
#include "smesh_read.hpp"
#include "smesh_refine.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_tracer.hpp"
#include "smesh_write.hpp"
#include "smesh_volume_to_surface.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_distributed_write.hpp"
#include "smesh_decompose.hpp"
#endif

#ifdef SMESH_ENABLE_RYAML
#include <ryml.hpp>
#include <ryml_std.hpp>
#endif

#include <algorithm>
#include <array>
#include <fstream>
#include <list>
#include <map>
#include <math.h>
#include <unordered_map>
#include <vector>

namespace smesh {

class Distributed::Impl {
public:
  ptrdiff_t n_nodes_global = 0;
  ptrdiff_t n_nodes_owned = 0;
  ptrdiff_t n_nodes_shared = 0;
  ptrdiff_t n_nodes_ghosts = 0;
  ptrdiff_t n_nodes_aura = 0;

  ptrdiff_t n_elements_global = 0;
  ptrdiff_t n_elements_owned = 0;
  ptrdiff_t n_elements_shared = 0;
  ptrdiff_t n_elements_ghosts = 0;

  SharedBuffer<large_idx_t> node_mapping;
  SharedBuffer<large_idx_t> element_mapping;

  SharedBuffer<int> node_owner;
  SharedBuffer<ptrdiff_t> node_offsets;
  SharedBuffer<idx_t> ghosts_and_aura;

  // Only exists for aura elements (global id)
  SharedBuffer<large_idx_t> element_id_aura;
};

SharedBuffer<large_idx_t> Distributed::node_mapping() const {
  SMESH_ASSERT(impl_->node_mapping);
  return impl_->node_mapping;
}

SharedBuffer<large_idx_t> Distributed::element_mapping() const {
  return impl_->element_mapping;
}
SharedBuffer<int> Distributed::node_owner() const { return impl_->node_owner; }
SharedBuffer<ptrdiff_t> Distributed::node_offsets() const {
  return impl_->node_offsets;
}
SharedBuffer<idx_t> Distributed::ghosts() const {
  return view(impl_->ghosts_and_aura, 0, impl_->n_nodes_ghosts);
}
SharedBuffer<idx_t> Distributed::ghosts_and_aura() const {
  return impl_->ghosts_and_aura;
}

Distributed::Distributed() : impl_(std::make_unique<Impl>()) {}
Distributed::~Distributed() = default;

ptrdiff_t Distributed::n_nodes_global() const { return impl_->n_nodes_global; }
ptrdiff_t Distributed::n_elements_global() const {
  return impl_->n_elements_global;
}
ptrdiff_t Distributed::n_nodes_local() const {
  return impl_->n_nodes_owned + impl_->n_nodes_ghosts + impl_->n_nodes_aura;
}
ptrdiff_t Distributed::n_nodes_owned_not_shared() const {
  return impl_->n_nodes_owned - impl_->n_nodes_shared;
}
ptrdiff_t Distributed::n_nodes_owned() const { return impl_->n_nodes_owned; }
ptrdiff_t Distributed::n_nodes_shared() const { return impl_->n_nodes_shared; }
ptrdiff_t Distributed::n_nodes_ghosts() const { return impl_->n_nodes_ghosts; }
ptrdiff_t Distributed::n_nodes_aura() const { return impl_->n_nodes_aura; }
ptrdiff_t Distributed::n_elements_local() const {
  return impl_->n_elements_owned + impl_->n_elements_ghosts;
}
ptrdiff_t Distributed::n_elements_owned_not_shared() const {
  return impl_->n_elements_owned - impl_->n_elements_shared;
}
ptrdiff_t Distributed::n_elements_owned() const {
  return impl_->n_elements_owned;
}
ptrdiff_t Distributed::n_elements_shared() const {
  return impl_->n_elements_shared;
}
ptrdiff_t Distributed::n_elements_ghosts() const {
  return impl_->n_elements_ghosts;
}

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

Mesh::Block::Block(const std::string &name, enum ElemType element_type,
                   SharedBuffer<idx_t *> elements)
    : impl_(std::make_unique<Impl>()) {
  impl_->name = name;
  impl_->element_type = element_type;
  impl_->elements = elements;
}
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
  SharedBuffer<geom_t *> points;
  SharedBuffer<idx_t> node_mapping;

  std::shared_ptr<Distributed> distributed;
  std::shared_ptr<NodeToNodeGraph> crs_graph;
  std::shared_ptr<NodeToNodeGraph> crs_graph_upper_triangular;
  std::shared_ptr<NodeToElementGraph> node_to_element_graph;

  ~Impl() {}

  void clear() {
    comm = nullptr;
    blocks.clear();
    points = nullptr;
    distributed = nullptr;
    crs_graph = nullptr;
    crs_graph_upper_triangular = nullptr;
    distributed = nullptr;
    node_to_element_graph = nullptr;
  }

  // Helper methods for backward compatibility
  ptrdiff_t total_elements() const {
    ptrdiff_t total = 0;
    for (const auto &block : blocks) {
      if (block && block->elements()) {
        total += block->n_elements();
      }
    }
    return total;
  }

  void create_node_to_element_graph() {
    if (node_to_element_graph) {
      return;
    }
    node_to_element_graph = std::make_shared<NodeToElementGraph>();

    count_t *rowptr{nullptr};
    element_idx_t *colidx{nullptr};

    const ptrdiff_t nnodes = points->extent(1);
    if (blocks.size() == 1) {
      auto block0 = blocks[0];
      create_n2e(block0->n_elements(), nnodes, block0->n_nodes_per_element(),
                 block0->elements()->data(), &rowptr, &colidx);
    } else {
      // Multiblock: build node-to-element graph over all blocks.
      std::vector<enum ElemType> element_types;
      std::vector<ptrdiff_t> n_elements;
      std::vector<idx_t **> elements;

      for (auto &block : blocks) {
        if (!block || !block->elements()) {
          continue;
        }
        element_types.push_back(block->element_type());
        n_elements.push_back(block->elements()->extent(1));
        elements.push_back(block->elements()->data());
      }

      create_multiblock_n2e<idx_t, count_t, element_idx_t>(
          static_cast<block_idx_t>(element_types.size()), element_types.data(),
          n_elements.data(), elements.data(), nnodes, nullptr, &rowptr,
          &colidx);
    }

    node_to_element_graph = std::make_shared<Mesh::NodeToElementGraph>(
        Buffer<count_t>::own(nnodes + 1, rowptr, free, MEMORY_SPACE_HOST),
        Buffer<element_idx_t>::own(rowptr[nnodes], colidx, free,
                                   MEMORY_SPACE_HOST));
  }
};

std::shared_ptr<Communicator> Mesh::comm() const { return impl_->comm; }

std::shared_ptr<Distributed> Mesh::distributed() const {
  SMESH_ASSERT(impl_->distributed);
  return impl_->distributed;
}

Mesh::Mesh(const std::shared_ptr<Communicator> &comm,
           enum ElemType element_type, SharedBuffer<idx_t *> elements,
           SharedBuffer<geom_t *> points)
    : impl_(std::make_unique<Impl>()) {
  impl_->comm = comm;
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

void Mesh::add_block(const std::shared_ptr<Block> &block) {
  impl_->blocks.push_back(block);
}

void Mesh::remove_block(size_t index) {
  if (index >= impl_->blocks.size()) {
    SMESH_ERROR("Block index out of range");
  }

  impl_->blocks.erase(impl_->blocks.begin() + index);
}

void read_meta(const std::shared_ptr<Communicator> &comm, const Path &path,
               enum ElemType &element_type) {

  if (!comm->rank()) {
    auto meta_file = Path(path) / "meta.yaml";
    if (meta_file.exists()) {
#if defined(SMESH_ENABLE_RYAML)
      std::ifstream ifs(meta_file.c_str(), std::ios::binary);
      if (ifs.good()) {
        std::string yaml((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
        if (!yaml.empty()) {
          ryml::Tree tree =
              ryml::parse_in_arena(ryml::to_csubstr(yaml)); // modifies input
          auto root = tree.rootref();
          if (root.has_child("element_type")) {
            auto v = root["element_type"].val();
            std::string s(v.str, v.len);
            element_type = type_from_string(s.c_str());
          }
        }
      }
#else
      std::ifstream ifs(meta_file.c_str());
      while (ifs.good()) {
        std::string line;
        std::getline(ifs, line);
        if (line.find("element_type:") != std::string::npos) {
          auto element_type_str = trim(line.substr(line.find(":") + 1));
          element_type = type_from_string(element_type_str.c_str());
          break;
        }
      }
#endif
    }
  }

  if (comm->size() > 1) {
    int element_type_int = (int)element_type;
    comm->broadcast(&element_type_int, 1, 0);
    element_type = (enum ElemType)element_type_int;
  }
}

static bool read_blocks_meta(const Path &path,
                             std::vector<std::string> &block_names,
                             std::vector<enum ElemType> &element_types) {
  block_names.clear();
  element_types.clear();

  auto meta_file = Path(path) / "meta.yaml";
  if (!meta_file.exists()) {
    return false;
  }

#if defined(SMESH_ENABLE_RYAML)
  std::ifstream ifs(meta_file.c_str(), std::ios::binary);
  if (!ifs.good()) {
    return false;
  }

  std::string yaml((std::istreambuf_iterator<char>(ifs)),
                   std::istreambuf_iterator<char>());
  if (yaml.empty()) {
    return false;
  }

  ryml::Tree tree = ryml::parse_in_arena(ryml::to_csubstr(yaml));
  auto root = tree.rootref();
  if (!root.has_child("blocks")) {
    return false;
  }

  auto blocks = root["blocks"];
  const size_t n = blocks.num_children();
  block_names.reserve(n);
  element_types.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    auto blk = blocks[i];
    if (!blk.has_child("name")) {
      continue;
    }

    auto name_val = blk["name"].val();
    std::string name(name_val.str, name_val.len);

    enum ElemType et = INVALID;
    if (blk.has_child("element_type")) {
      auto et_val = blk["element_type"].val();
      std::string et_str(et_val.str, et_val.len);
      et = type_from_string(et_str.c_str());
    }

    block_names.push_back(name);
    element_types.push_back(et);
  }

  return !block_names.empty();
#else
  (void)path;
  (void)block_names;
  (void)element_types;
  return false;
#endif
}

int Mesh::read(const Path &path) {
  SMESH_TRACE_SCOPE("Mesh::read");

  if (impl_->comm->size() == 1) {
    std::vector<std::string> block_names;
    std::vector<enum ElemType> element_types;
    const bool has_blocks = read_blocks_meta(path, block_names, element_types);

    impl_->blocks.clear();

    if (has_blocks) {
      // Shared coordinates for all blocks
      geom_t **points = nullptr;
      int spatial_dim = 0;
      ptrdiff_t nnodes = 0;
      if (mesh_coordinates_from_folder(path, &spatial_dim, &points, &nnodes) !=
          SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }
      impl_->points = manage_host_buffer<geom_t>(spatial_dim, nnodes, points);

      const size_t n_blocks = block_names.size();
      for (size_t b = 0; b < n_blocks; ++b) {
        int nnodesxelem = 0;
        ptrdiff_t nelements = 0;
        idx_t **elements = nullptr;

        Path block_folder = Path(path) / "blocks" / block_names[b];
        if (mesh_block_from_folder(block_folder, &nnodesxelem, &elements,
                                   &nelements) != SMESH_SUCCESS) {
          return SMESH_FAILURE;
        }

        auto elements_buffer =
            manage_host_buffer<idx_t>(nnodesxelem, nelements, elements);

        enum ElemType et = element_types[b];
        if (et == INVALID) {
          // Fallback: infer from number of nodes per element when metadata is
          // missing.
          et = (enum ElemType)nnodesxelem;
        }

        auto block = std::make_shared<Block>();
        block->set_name(block_names[b]);
        block->set_element_type(et);
        block->set_elements(elements_buffer);
        this->add_block(block);
      }
    } else {
      // Legacy single-block layout: connectivity and points live directly
      // under path.
      idx_t **elements = nullptr;
      geom_t **points = nullptr;
      int nnodesxelem;
      int spatial_dim;
      ptrdiff_t nnodes;
      ptrdiff_t nelements;

      if (mesh_from_folder(path, &nnodesxelem, &nelements, &elements,
                           &spatial_dim, &nnodes, &points) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
      }

      auto elements_buffer =
          manage_host_buffer<idx_t>(nnodesxelem, nelements, elements);
      impl_->points = manage_host_buffer<geom_t>(spatial_dim, nnodes, points);

      enum ElemType element_type = (enum ElemType)nnodesxelem;
      read_meta(impl_->comm, path, element_type);

      auto default_block = std::make_shared<Block>();
      default_block->set_name("default");
      default_block->set_element_type(element_type);
      default_block->set_elements(elements_buffer);
      this->add_block(default_block);
    }
  }
#ifdef SMESH_ENABLE_MPI
  else {

    auto dist = std::make_shared<Distributed>();
    int nnodesxelem;
    large_idx_t *element_mapping = nullptr;
    large_idx_t *node_mapping = nullptr;
    idx_t **elements = nullptr;
    geom_t **points = nullptr;
    ptrdiff_t *node_offsets = nullptr;
    int *node_owner = nullptr;
    idx_t *ghosts = nullptr;
    ptrdiff_t n_nodes_aura = 0;

    int spatial_dim;
    if (mesh_from_folder(
            impl_->comm->get(), path,
            // Elements
            &nnodesxelem, &dist->impl_->n_elements_global,
            &dist->impl_->n_elements_owned, &dist->impl_->n_elements_shared,
            &dist->impl_->n_elements_ghosts, &element_mapping, &elements,
            // Nodes
            &spatial_dim, &dist->impl_->n_nodes_global,
            &dist->impl_->n_nodes_owned, &dist->impl_->n_nodes_shared,
            &dist->impl_->n_nodes_ghosts, &n_nodes_aura, &node_mapping, &points,
            // Distributed connectivities
            &node_owner, &node_offsets, &ghosts) != SMESH_SUCCESS) {
      SMESH_ERROR("Failed to read mesh from folder %s\n", path.c_str());
      return SMESH_FAILURE;
    }
    dist->impl_->n_nodes_aura = n_nodes_aura;

    auto elements_buffer = manage_host_buffer<idx_t>(
        nnodesxelem, dist->n_elements_local(), elements);
    impl_->points =
        manage_host_buffer<geom_t>(spatial_dim, dist->n_nodes_local(), points);
    dist->impl_->node_mapping =
        manage_host_buffer<large_idx_t>(dist->n_nodes_local(), node_mapping);
    dist->impl_->node_owner =
        manage_host_buffer<int>(dist->n_nodes_local(), node_owner);
    dist->impl_->element_mapping = manage_host_buffer<large_idx_t>(
        dist->n_elements_owned(), element_mapping);

    int comm_size;
    MPI_Comm_size(impl_->comm->get(), &comm_size);
    dist->impl_->node_offsets =
        manage_host_buffer<ptrdiff_t>(comm_size + 1, node_offsets);

    dist->impl_->ghosts_and_aura = manage_host_buffer<idx_t>(
        dist->impl_->n_nodes_ghosts + dist->impl_->n_nodes_aura, ghosts);

    // Best effort for basic types
    enum ElemType element_type = (enum ElemType)nnodesxelem;
    read_meta(impl_->comm, path, element_type);

    // Create default block
    auto default_block = std::make_shared<Block>();
    default_block->set_name("default");
    default_block->set_element_type(element_type);
    default_block->set_elements(elements_buffer);
    this->add_block(default_block);

    impl_->distributed = dist;
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

  if (impl_->comm->size() == 1) {
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
                            this->spatial_dimension(), this->n_nodes(),
                            this->points()->data());
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
          this->spatial_dimension(), this->n_nodes(), this->points()->data());
    }
  }
#ifdef SMESH_ENABLE_MPI
  else {
    // Parallel topology write for single-block meshes.
    if (!impl_->distributed) {
      SMESH_ERROR("Mesh::write (MPI) requires a Distributed object. "
                  "Did you create the mesh via distributed read?\n");
    }

    if (impl_->blocks.size() != 1 || !impl_->blocks[0]) {
      SMESH_ERROR(
          "Mesh::write (MPI) currently supports only single-block meshes.\n");
    }

    auto dist = impl_->distributed;

    const enum ElemType et = impl_->blocks[0]->element_type();
    const int nxe = elem_num_nodes(et);

    // Write coordinates and connectivity in parallel.
    int err = write_distributed_mesh_topology(
        impl_->comm->get(), path, et, this->spatial_dimension(),
        dist->n_elements_global(), dist->n_elements_owned(),
        dist->impl_->element_mapping->data(), nxe,
        impl_->blocks[0]->elements()->data(), dist->n_nodes_global(),
        dist->n_nodes_owned(), dist->impl_->node_mapping->data(),
        impl_->points->data());

    if (err != SMESH_SUCCESS) {
      return err;
    }

    return SMESH_SUCCESS;
  }
#endif

  return SMESH_FAILURE;
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
  if (n_blocks() > 1) {
    SMESH_ERROR("half_face_table is not supported for multiblock meshes");
    return nullptr;
  }

  const block_idx_t block_id = 0;
  // FIXME it should be allocated outisde
  element_idx_t *table{nullptr};
  create_element_adj_table(n_elements(block_id), n_nodes(),
                           element_type(block_id), elements(block_id)->data(),
                           &table);

  int nsxe = elem_num_sides(element_type(block_id));
  return manage_host_buffer<element_idx_t>(n_elements(block_id) * nsxe, table);
}

// FIXME: redesign to support multiblock meshes
std::shared_ptr<Mesh::NodeToNodeGraph>
Mesh::create_node_to_node_graph(const enum ElemType element_type) {
  if (n_blocks() != 1) {
    SMESH_ERROR(
        "create_node_to_node_graph is not supported for multi-block meshes!\n");
    return nullptr;
  }

  if (this->element_type(0) == element_type) {
    return node_to_node_graph();
  }

  const ptrdiff_t n_nodes =
      max_node_id(element_type, n_elements(0), elements(0)->data()) + 1;

  count_t *rowptr{nullptr};
  idx_t *colidx{nullptr};
  if (is_semistructured_type(this->element_type(0))) {

    // TODO: check of it works
    SMESH_ERROR("Semistructured meshes by create_node_to_node_graph for "
                "different element type!\n");
    // for other semistructured elements
    sshex8_crs_graph<element_idx_t, count_t, idx_t>(
        proteus_hex_micro_elements_per_dim(element_type), this->n_elements(0),
        this->n_nodes(), this->elements(0)->data(), &rowptr, &colidx);

  } else {

    create_crs_graph_for_elem_type(element_type, n_elements(0), n_nodes,
                                   elements(0)->data(), &rowptr, &colidx);
  }

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
    if (is_semistructured_type(this->element_type(0))) {

      sshex8_crs_graph<element_idx_t, count_t, idx_t>(
          proteus_hex_micro_elements_per_dim(this->element_type(0)),
          this->n_elements(0), this->n_nodes(), this->elements(0)->data(),
          &rowptr, &colidx);

    } else {

      create_crs_graph_for_elem_type(this->element_type(0), this->n_elements(0),
                                     this->n_nodes(), this->elements(0)->data(),
                                     &rowptr, &colidx);
    }
  } else {

    if (is_semistructured_type(this->element_type(0))) {
      SMESH_ERROR(
          "Semistructured meshes are not supported for multi-block meshes");
      return SMESH_FAILURE;
    }
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
                                this->n_nodes(), &rowptr, &colidx);
  }

  impl_->crs_graph = std::make_shared<Mesh::NodeToNodeGraph>(
      Buffer<count_t>::own(this->n_nodes() + 1, rowptr, free,
                           MEMORY_SPACE_HOST),
      Buffer<idx_t>::own(rowptr[this->n_nodes()], colidx, free,
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
        impl_->total_elements(), this->n_nodes(),
        elem_num_nodes(this->element_type(0)), this->elements(0)->data(),
        &rowptr, &colidx);
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
        elements.data(), this->n_nodes(), &rowptr, &colidx);
  }

  impl_->crs_graph_upper_triangular = std::make_shared<Mesh::NodeToNodeGraph>(
      Buffer<count_t>::own(this->n_nodes() + 1, rowptr, free,
                           MEMORY_SPACE_HOST),
      Buffer<idx_t>::own(rowptr[this->n_nodes()], colidx, free,
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
Mesh::create_hex8_cube(const std::shared_ptr<Communicator> &comm,
                       const ptrdiff_t nx, const ptrdiff_t ny,
                       const ptrdiff_t nz, const geom_t xmin, const geom_t ymin,
                       const geom_t zmin, const geom_t xmax, const geom_t ymax,
                       const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);

  ret->impl_->points = create_host_buffer<geom_t>(3, nnodes);
  auto elements_buffer = create_host_buffer<idx_t>(8, nelements);

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_hex8_cube<idx_t, geom_t>(nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                                     zmax, elements, points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(HEX8);
  default_block->set_elements(elements_buffer);
  ret->add_block(default_block);

  return ret;
}

std::shared_ptr<Mesh> Mesh::create_semistructured_hex_cube(
    const std::shared_ptr<Communicator> &comm, const int micro_elements_per_dim,
    const ptrdiff_t nx, const ptrdiff_t ny, const ptrdiff_t nz,
    const geom_t xmin, const geom_t ymin, const geom_t zmin, const geom_t xmax,
    const geom_t ymax, const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = (nx) * (ny) * (nz);

  const int micro_nodes_per_dim = (micro_elements_per_dim + 1);
  const ptrdiff_t nnodes = (nx * micro_nodes_per_dim) *
                           (ny * micro_nodes_per_dim) *
                           (nz * micro_nodes_per_dim);

  const int nxme =
      micro_nodes_per_dim * micro_nodes_per_dim * micro_nodes_per_dim;
  ret->impl_->points = create_host_buffer<geom_t>(3, nnodes);
  auto elements = create_host_buffer<idx_t>(nxme, nelements);

  mesh_fill_proteus_hex_cube<idx_t, geom_t>(
      micro_elements_per_dim, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax,
      elements->data(), ret->impl_->points->data());

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(proteus_hex_type(micro_elements_per_dim));
  default_block->set_elements(elements);
  ret->add_block(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_tri3_square(const std::shared_ptr<Communicator> &comm,
                         const ptrdiff_t nx, const ptrdiff_t ny,
                         const geom_t xmin, const geom_t ymin,
                         const geom_t xmax, const geom_t ymax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = 2 * nx * ny;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1);

  ret->impl_->points = create_host_buffer<geom_t>(2, nnodes);
  auto elements_buffer = create_host_buffer<idx_t>(3, nelements);

  auto points = ret->impl_->points->data();
  auto elements = elements_buffer->data();

  mesh_fill_tri3_square<idx_t, geom_t>(nx, ny, xmin, ymin, xmax, ymax, elements,
                                       points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(TRI3);
  default_block->set_elements(elements_buffer);
  ret->add_block(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_quad4_square(const std::shared_ptr<Communicator> &comm,
                          const ptrdiff_t nx, const ptrdiff_t ny,
                          const geom_t xmin, const geom_t ymin,
                          const geom_t xmax, const geom_t ymax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1);

  ret->impl_->points = create_host_buffer<geom_t>(2, nnodes);
  auto elements_buffer = create_host_buffer<idx_t>(4, nelements);

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
                    const enum ElemType element_type, const ptrdiff_t nx,
                    const ptrdiff_t ny, const geom_t xmin, const geom_t ymin,
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
                        const ptrdiff_t nlayers, const ptrdiff_t nelements) {
  auto elements = create_host_buffer<idx_t>(4, nlayers * nelements);
  auto points = create_host_buffer<geom_t>(3, (nlayers + 1) * nelements);
  mesh_fill_quad4_ring<idx_t, geom_t>(inner_radius, outer_radius, nlayers,
                                      nelements, elements->data(),
                                      points->data());

  auto ret = std::make_shared<Mesh>(comm, QUAD4, elements, points);
  return ret;
}

std::shared_ptr<Mesh> Mesh::create_hex8_checkerboard_cube(
    const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx,
    const ptrdiff_t ny, const ptrdiff_t nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);

  if (nx % 2 != 0 || ny % 2 != 0 || nz % 2 != 0) {
    SMESH_ERROR("nx, ny, and nz must be even");
  }

  ret->set_points(create_host_buffer<geom_t>(3, nnodes));
  auto white_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);
  auto black_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);

  auto points = ret->points()->data();
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
  ret->add_block(white_block);

  auto black_block = std::make_shared<Block>();
  black_block->set_name("black");
  black_block->set_element_type(HEX8);
  black_block->set_elements(black_elements_buffer);
  ret->add_block(black_block);
  return ret;
}

std::shared_ptr<Mesh> Mesh::create_hex8_bidomain_cube(
    const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx,
    const ptrdiff_t ny, const ptrdiff_t nz, const geom_t xmin,
    const geom_t ymin, const geom_t zmin, const geom_t xmax, const geom_t ymax,
    const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);

  ret->set_points(create_host_buffer<geom_t>(3, nnodes));
  auto left_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);
  auto right_elements_buffer = create_host_buffer<idx_t>(8, nelements / 2);

  auto points = ret->points()->data();
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
  ret->add_block(left_block);

  auto right_block = std::make_shared<Block>();
  right_block->set_name("right");
  right_block->set_element_type(HEX8);
  right_block->set_elements(right_elements_buffer);
  ret->add_block(right_block);

  return ret;
}

int Mesh::spatial_dimension() const { return points()->extent(0); }

ptrdiff_t Mesh::n_nodes() const { return points()->extent(1); }
ptrdiff_t Mesh::n_elements() const { return impl_->total_elements(); }

int Mesh::n_nodes_per_element(block_idx_t block_id) const {
  auto blk = this->block(block_id);
  SMESH_ASSERT(blk);
  return blk->n_nodes_per_element();
}

ptrdiff_t Mesh::n_elements(block_idx_t block_id) const {
  auto blk = this->block(block_id);
  SMESH_ASSERT(blk);
  return blk->n_elements();
}

enum ElemType Mesh::element_type(block_idx_t block_id) const {
  auto blk = this->block(block_id);
  SMESH_ASSERT(blk);
  return blk->element_type();
}

SharedBuffer<geom_t *> Mesh::points() { return impl_->points; }
SharedBuffer<geom_t *> Mesh::points() const { return impl_->points; }

void Mesh::set_points(const SharedBuffer<geom_t *> &points) {
  impl_->points = points;
}

SharedBuffer<idx_t *> Mesh::elements(block_idx_t block_id) {
  auto blk = this->block(block_id);
  SMESH_ASSERT(blk);
  return blk->elements();
}

SharedBuffer<idx_t *> Mesh::elements(block_idx_t block_id) const {
  auto blk = this->block(block_id);
  SMESH_ASSERT(blk);
  return blk->elements();
}

void Mesh::set_node_mapping(const SharedBuffer<idx_t> &node_mapping) {
  impl_->node_mapping = node_mapping;
}

void Mesh::set_comm(const std::shared_ptr<Communicator> &comm) {
  impl_->comm = comm;
}

void Mesh::set_element_type(const block_idx_t block_id,
                            const enum ElemType element_type) {
  auto blk = this->block(block_id);
  SMESH_ASSERT(blk);
  blk->set_element_type(element_type);
}

std::vector<std::shared_ptr<Mesh::Block>>
Mesh::blocks(const std::vector<std::string> &block_names) const {
  if (block_names.empty()) {
    return this->blocks();
  }

  std::vector<std::shared_ptr<Mesh::Block>> ret;
  for (auto &block : this->blocks()) {
    if (std::find(block_names.begin(), block_names.end(), block->name()) !=
        block_names.end()) {
      ret.push_back(block);
    }
  }

  return ret;
}

std::shared_ptr<Mesh> Mesh::create_hex8_reference_cube() {
  auto ret = std::make_shared<Mesh>(Communicator::null());
  ret->set_points(create_host_buffer<geom_t>(3, 8));
  auto elements_buffer = create_host_buffer<idx_t>(8, 1);

  auto points = ret->points()->data();
  auto elements = elements_buffer->data();

  mesh_fill_hex8_reference_cube<idx_t, geom_t>(elements, points);

  // Create default block
  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(HEX8);
  default_block->set_elements(elements_buffer);
  ret->add_block(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_tet4_cube(const std::shared_ptr<Communicator> &comm,
                       const ptrdiff_t nx, const ptrdiff_t ny,
                       const ptrdiff_t nz, const geom_t xmin, const geom_t ymin,
                       const geom_t zmin, const geom_t xmax, const geom_t ymax,
                       const geom_t zmax) {
  auto ret = std::make_shared<Mesh>(comm);
  const ptrdiff_t nelements = nx * ny * nz;
  const ptrdiff_t nnodes_vertices = (nx + 1) * (ny + 1) * (nz + 1);
  const ptrdiff_t nnodes_total = nnodes_vertices + nelements;

  ret->set_points(create_host_buffer<geom_t>(3, nnodes_total));
  auto elements_buffer = create_host_buffer<idx_t>(4, nelements * 12);

  auto points = ret->points()->data();
  auto elements = elements_buffer->data();

  mesh_fill_tet4_cube<idx_t, geom_t>(nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                                     zmax, elements, points);

  auto default_block = std::make_shared<Block>();
  default_block->set_name("default");
  default_block->set_element_type(TET4);
  default_block->set_elements(elements_buffer);
  ret->add_block(default_block);

  return ret;
}

std::shared_ptr<Mesh>
Mesh::create_cube(const std::shared_ptr<Communicator> &comm,
                  const enum ElemType element_type, const ptrdiff_t nx,
                  const ptrdiff_t ny, const ptrdiff_t nz, const geom_t xmin,
                  const geom_t ymin, const geom_t zmin, const geom_t xmax,
                  const geom_t ymax, const geom_t zmax) {
  SMESH_TRACE_SCOPE("Mesh::create_cube");
  switch (element_type) {
  case HEX8:
    return create_hex8_cube(comm, nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                            zmax);
  case TET4:
    return create_tet4_cube(comm, nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                            zmax);
  case PROTEUS_HEX8:
    return create_semistructured_hex_cube(comm, 1, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX27:
    return create_semistructured_hex_cube(comm, 2, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX64:
    return create_semistructured_hex_cube(comm, 3, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX125:
    return create_semistructured_hex_cube(comm, 4, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX216:
    return create_semistructured_hex_cube(comm, 5, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX343:
    return create_semistructured_hex_cube(comm, 6, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX512:
    return create_semistructured_hex_cube(comm, 7, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
  case PROTEUS_HEX729:
    return create_semistructured_hex_cube(comm, 8, nx, ny, nz, xmin, ymin, zmin,
                                          xmax, ymax, zmax);
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
    const int nxe = n_nodes_per_element(0);
    const ptrdiff_t n_elements = this->n_elements(0);

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
          SMESH_ASSERT(n_interior_elements_count <
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
      this->add_block(block);
    }

    { // Interior block
      auto block = std::make_shared<Block>();
      block->set_name(default_block->name());
      block->set_element_type(default_block->element_type());
      block->set_elements(interior_elements);
      this->add_block(block);
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

    if (extract_skin_sideset(this->n_elements(0), this->n_nodes(),
                             this->element_type(0), this->elements(0)->data(),
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
    auto nxe = b->n_nodes_per_element();

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
}

int Mesh::renumber_nodes(const SharedBuffer<idx_t> &node_mapping) {
  const int dim = spatial_dimension();
  const ptrdiff_t n_nodes = this->n_nodes();

  auto points = this->points()->data();
  auto new_points_buff = create_host_buffer<geom_t>(dim, n_nodes);
  auto new_points = new_points_buff->data();

  auto d_node_mapping = node_mapping->data();

  for (int d = 0; d < dim; d++) {
    for (ptrdiff_t i = 0; i < n_nodes; i++) {
      SMESH_ASSERT(d_node_mapping[i] < n_nodes);
      SMESH_ASSERT(d_node_mapping[i] >= 0);
      new_points[d][d_node_mapping[i]] = points[d][i];
    }
  }

  impl_->points = new_points_buff;

  for (auto &b : impl_->blocks) {
    auto elements = b->elements()->data();
    auto n_elements = b->n_elements();
    auto nxe = b->n_nodes_per_element();

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

  const int dim = spatial_dimension();
  auto points = this->points()->data();

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

    int nxe = block->n_nodes_per_element();
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

void Mesh::reorder_elements_from_tags(const block_idx_t block_id,
                                      const SharedBuffer<idx_t> &tags) {
  const ptrdiff_t nelems = n_elements(block_id);
  auto temp = create_host_buffer<idx_t>(nelems);
  auto d_temp = temp->data();
  auto d_tags = tags->data();

  auto d_elements = elements(block_id)->data();

  idx_t ntags = 0;
  for (ptrdiff_t i = 0; i < nelems; i++) {
    ntags = std::max(ntags, d_tags[i]);
  }

  if (!ntags)
    return;

  ntags += 1;

  auto bookkeeping = create_host_buffer<ptrdiff_t>(ntags);
  auto d_bk = bookkeeping->data();

  int nxe = n_nodes_per_element(block_id);
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
    SMESH_ERROR("Conversion from %s is not supported for multiblock meshes\n",
                type_to_string(element_type));
    return nullptr;
  }

  std::map<std::pair<enum ElemType, enum ElemType>,
           std::function<void(const Mesh::Block &, Mesh::Block &)>>
      cmap;

  cmap[std::make_pair(HEX8, TET4)] = [](const Mesh::Block &block,
                                        Mesh::Block &new_block) {
    new_block.set_element_type(TET4);
    new_block.set_elements(
        create_host_buffer<idx_t>(4, block.n_elements() * 6));
    mesh_hex8_to_6x_tet4(block.n_elements(), block.elements()->data(),
                         new_block.elements()->data());
  };

  cmap[std::make_pair(TET15, HEX8)] = [](const Mesh::Block &block,
                                         Mesh::Block &new_block) {
    new_block.set_element_type(HEX8);
    new_block.set_elements(
        create_host_buffer<idx_t>(8, block.n_elements() * 4));
    mesh_tet15_to_4x_hex8(block.n_elements(), block.elements()->data(),
                          new_block.elements()->data());
  };

  cmap[std::make_pair(WEDGE6, TET4)] = [](const Mesh::Block &block,
                                          Mesh::Block &new_block) {
    new_block.set_element_type(TET4);
    new_block.set_elements(
        create_host_buffer<idx_t>(4, block.n_elements() * 3));
    mesh_wedge6_to_3x_tet4(block.n_elements(), block.elements()->data(),
                           new_block.elements()->data());
  };

  cmap[std::make_pair(HEX8, PROTEUS_HEX8)] = [](const Mesh::Block &block,
                                                Mesh::Block &new_block) {
    new_block.set_element_type(PROTEUS_HEX8);
    auto elements = block.elements();

    auto view = std::make_shared<Buffer<idx_t *>>(
        8, block.n_elements(), (idx_t **)malloc(8 * sizeof(idx_t *)),
        [keep_alive = elements](int, void **v) {
          (void)keep_alive;
          free(v);
        },
        elements->mem_space());

    const int pts[8] = {// Bottom
                        sshex8_lidx(1, 0, 0, 0), sshex8_lidx(1, 1, 0, 0),
                        sshex8_lidx(1, 1, 1, 0), sshex8_lidx(1, 0, 1, 0),

                        // Top
                        sshex8_lidx(1, 0, 0, 1), sshex8_lidx(1, 1, 0, 1),
                        sshex8_lidx(1, 1, 1, 1), sshex8_lidx(1, 0, 1, 1)};

    view->data()[0] = elements->data()[pts[0]];
    view->data()[1] = elements->data()[pts[1]];
    view->data()[2] = elements->data()[pts[2]];
    view->data()[3] = elements->data()[pts[3]];
    view->data()[4] = elements->data()[pts[4]];
    view->data()[5] = elements->data()[pts[5]];
    view->data()[6] = elements->data()[pts[6]];
    view->data()[7] = elements->data()[pts[7]];

    new_block.set_elements(view);
  };

  cmap[std::make_pair(PROTEUS_HEX8, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX27, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX64, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX125, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX216, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX343, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX512, HEX8)] = sshex_block_to_hex8_block;
  cmap[std::make_pair(PROTEUS_HEX729, HEX8)] = sshex_block_to_hex8_block;

  cmap[std::make_pair(PROTEUS_QUAD4, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD9, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD16, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD25, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD36, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD49, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD64, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUAD81, QUAD4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL4, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL9, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL16, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL25, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL36, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL49, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL64, QUADSHELL4)] = ssquad_block_to_quad4_block;
  cmap[std::make_pair(PROTEUS_QUADSHELL81, QUADSHELL4)] = ssquad_block_to_quad4_block;

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
    SMESH_ERROR("Promotion to %s is not supported for multiblock meshes\n",
                type_to_string(element_type));
    return nullptr;
  }

  std::map<std::pair<enum ElemType, enum ElemType>,
           std::function<std::shared_ptr<Mesh>(Mesh &)>>
      cmap;

  cmap[std::make_pair(TET4, TET15)] = [](Mesh &mesh) -> std::shared_ptr<Mesh> {
    auto elements = create_host_buffer<idx_t>(15, mesh.n_elements(0));
    auto n2n_upper_triangular = mesh.node_to_node_graph_upper_triangular();
    auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
    auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

    auto hft = mesh.half_face_table();
    auto e2e_table = hft->data();

    ptrdiff_t n_new_nodes = 0;
    mesh_tet4_to_tet15(mesh.n_elements(0), mesh.n_nodes(),
                       mesh.elements(0)->data(), n2n_upper_triangular_ptr,
                       n2n_upper_triangular_idx, e2e_table, elements->data(),
                       &n_new_nodes);

    auto points = create_host_buffer<geom_t>(3, n_new_nodes);
    mesh_tet4_to_tet15_points(mesh.n_elements(0), mesh.n_nodes(),
                              mesh.points()->data(), n2n_upper_triangular_ptr,
                              n2n_upper_triangular_idx, elements->data(),
                              points->data());

    return std::make_shared<Mesh>(mesh.comm(), TET15, elements, points);
  };

  cmap[std::make_pair(TET4, TET10)] = [](Mesh &mesh) {
    auto n2n_upper_triangular = mesh.node_to_node_graph_upper_triangular();
    auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
    auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

    auto elements = create_host_buffer<idx_t>(10, mesh.n_elements(0));
    auto points = create_host_buffer<geom_t>(
        mesh.spatial_dimension(),
        n2n_upper_triangular->colidx()->size() + mesh.n_nodes());

    p1_to_p2(TET4, mesh.n_elements(0), mesh.elements(0)->data(),
             mesh.spatial_dimension(), mesh.n_nodes(), mesh.points()->data(),
             n2n_upper_triangular_ptr, n2n_upper_triangular_idx,
             elements->data(), points->data());
    return std::make_shared<Mesh>(mesh.comm(), TET10, elements, points);
  };

  cmap[std::make_pair(TRI3, TRI6)] = [](Mesh &mesh) {
    auto n2n_upper_triangular = mesh.node_to_node_graph_upper_triangular();
    auto n2n_upper_triangular_ptr = n2n_upper_triangular->rowptr()->data();
    auto n2n_upper_triangular_idx = n2n_upper_triangular->colidx()->data();

    auto elements = create_host_buffer<idx_t>(6, mesh.n_elements(0));
    auto points = create_host_buffer<geom_t>(
        mesh.spatial_dimension(), n2n_upper_triangular->colidx()->size());

    p1_to_p2(TRI3, mesh.n_elements(0), mesh.elements(0)->data(),
             mesh.spatial_dimension(), mesh.n_nodes(), mesh.points()->data(),
             n2n_upper_triangular_ptr, n2n_upper_triangular_idx,
             elements->data(), points->data());
    return std::make_shared<Mesh>(mesh.comm(), TRI6, elements, points);
  };

  auto it = cmap.find(std::make_pair(mesh->element_type(0), element_type));
  if (it != cmap.end()) {
    return it->second(*mesh);
  } else {
    SMESH_ERROR("Promotion from %s to %s is not supported\n",
                type_to_string(mesh->element_type(0)),
                type_to_string(element_type));
    return nullptr;
  }
}

std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                             const int levels) {
  if (mesh->n_blocks() != 1) {
    SMESH_ERROR("Refinement is not supported for multiblock meshes\n");
    return nullptr;
  }

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
  }(mesh->element_type(0));

  auto out = mesh;
  if (mesh->element_type(0) == HEX8) {
    const ptrdiff_t n_elements = mesh->n_elements(0);

    const int ss_levels = pow(2, levels);
    const int nxe = sshex8_nxe(ss_levels);
    const int txe = sshex8_txe(ss_levels);

    auto sshex8_elements = create_host_buffer<idx_t>(nxe, mesh->n_elements(0));
    auto d_sshex8_elements = sshex8_elements->data();

    ptrdiff_t n_unique_nodes = 0;
    ptrdiff_t interior_start = 0;

    sshex8_generate_elements(ss_levels, n_elements, mesh->n_nodes(),
                             mesh->elements(0)->data(), d_sshex8_elements,
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
          out->n_nodes_per_element(0), out->n_elements(0) * refine_factor);

      auto refined_points = create_host_buffer<geom_t>(
          out->spatial_dimension(),
          n2n_upper_triangular->colidx()->size() + out->n_nodes());

      int err = mesh_refine(out->element_type(0), out->n_elements(0),
                            out->elements(0)->data(), out->spatial_dimension(),
                            out->n_nodes(), out->points()->data(),
                            n2n_upper_triangular_ptr, n2n_upper_triangular_idx,
                            refined_elements->data(), refined_points->data());

      if (err != SMESH_SUCCESS) {
        SMESH_ERROR("Refinement failed\n");
        return nullptr;
      }

      out = std::make_shared<Mesh>(out->comm(), out->element_type(0),
                                   refined_elements, refined_points);
    }
  }

  return out;
}

std::shared_ptr<Sideset> skin_sideset(const std::shared_ptr<Mesh> &mesh) {
  SMESH_TRACE_SCOPE("skin_sideset");
  if (mesh->n_blocks() != 1) {
    SMESH_ERROR("Skin sideset is not supported for multiblock meshes\n");
    return nullptr;
  }

  if(is_semistructured_type(mesh->element_type(0))) {
    return skin_sideset(derefine(mesh, 1));
  }

  auto n2e_graph = mesh->node_to_element_graph();
  auto n2e_graph_ptr = n2e_graph->rowptr()->data();
  auto n2e_graph_idx = n2e_graph->colidx()->data();

  ptrdiff_t n_surf_elements = 0;
  element_idx_t *parent_element = 0;
  i16 *side_idx = 0;

  int err = extract_skin_sideset_from_n2e(mesh->n_elements(0), mesh->n_nodes(),
                                          mesh->element_type(0),
                                          mesh->elements(0)->data(),
                                          n2e_graph_ptr, n2e_graph_idx,
                                          &n_surf_elements, &parent_element,
                                          &side_idx);

  if (err != SMESH_SUCCESS) {
    SMESH_ERROR("Unable to extract skin sideset!\n");
    return nullptr;
  }

  if(mesh->comm()->size() > 1) {
    const auto dist = mesh->distributed();
    const auto n_owned_elements = dist->n_elements_owned();

    ptrdiff_t write_pos = 0;
    for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
      const element_idx_t parent = parent_element[i];
      if (parent >= n_owned_elements) {
        continue;
      }

      parent_element[write_pos] = parent;
      side_idx[write_pos] = side_idx[i];
      ++write_pos;
    }
    n_surf_elements = write_pos;

#if 0
    // FIXME: AI Slop to be clean-up
    if (n_surf_elements > 0) {
      LocalSideTable lst;
      lst.fill(mesh->element_type(0));
      const int ns = elem_num_sides(mesh->element_type(0));
      const int nnxs = elem_num_nodes(side_type(mesh->element_type(0)));
      auto elems = mesh->elements(0)->data();
      auto node_mapping = dist->node_mapping()->data();

      using FaceKey =
          std::array<large_idx_t, LocalSideTable::MAX_NUM_NODES_PER_SIDE>;

      // Build a canonical key for each candidate surface face on this rank by
      // mapping its local node ids to distributed/global node ids and sorting
      // them so the key is independent of side orientation.
      std::vector<FaceKey> local_face_keys((size_t)n_surf_elements);
      for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
        auto &key = local_face_keys[(size_t)i];
        const element_idx_t parent = parent_element[i];
        const int side = side_idx[i];

        for (int n = 0; n < nnxs; ++n) {
          const idx_t node = elems[lst(side, n)][parent];
          key[(size_t)n] = node_mapping[node];
        }

        std::sort(key.begin(), key.begin() + nnxs);
      }

      // Enumerate every side of every owned volume element on this rank using
      // the same canonical representation. These are the faces this rank is
      // allowed to contribute to the final distributed sideset.
      const ptrdiff_t n_owned_faces =
          static_cast<ptrdiff_t>(n_owned_elements) * ns;
      std::vector<large_idx_t> local_owned_face_key_data(
          (size_t)n_owned_faces * (size_t)nnxs);
      for (ptrdiff_t parent = 0; parent < n_owned_elements; ++parent) {
        for (int side = 0; side < ns; ++side) {
          std::array<large_idx_t, LocalSideTable::MAX_NUM_NODES_PER_SIDE> key;
          for (int n = 0; n < nnxs; ++n) {
            const idx_t node = elems[lst(side, n)][parent];
            key[(size_t)n] = node_mapping[node];
          }

          std::sort(key.begin(), key.begin() + nnxs);

          const ptrdiff_t face_idx = parent * ns + side;
          for (int n = 0; n < nnxs; ++n) {
            local_owned_face_key_data[(size_t)face_idx * (size_t)nnxs +
                                      (size_t)n] = key[(size_t)n];
          }
        }
      }

      // Gather the owned-face keys from all ranks so we can detect whether a
      // candidate surface face corresponds to exactly one globally owned face.
      std::vector<ptrdiff_t> local_counts(mesh->comm()->size());
      SMESH_MPI_CATCH(MPI_Allgather(&n_owned_faces, 1, mpi_type<ptrdiff_t>(),
                                    local_counts.data(), 1,
                                    mpi_type<ptrdiff_t>(),
                                    mesh->comm()->get()));

      std::vector<int> local_counts_i(mesh->comm()->size());
      std::vector<int> local_displs_i(mesh->comm()->size());
      ptrdiff_t total_faces = 0;
      ptrdiff_t total_face_entries = 0;
      for (int r = 0; r < mesh->comm()->size(); ++r) {
        local_counts_i[r] = static_cast<int>(local_counts[r] * nnxs);
        local_displs_i[r] = static_cast<int>(total_face_entries);
        total_faces += local_counts[r];
        total_face_entries += local_counts[r] * nnxs;
      }

      std::vector<large_idx_t> global_face_key_data(
          (size_t)total_face_entries);
      SMESH_MPI_CATCH(MPI_Allgatherv(
          local_owned_face_key_data.empty() ? nullptr
                                            : local_owned_face_key_data.data(),
          static_cast<int>(local_owned_face_key_data.size()),
          mpi_type<large_idx_t>(),
          global_face_key_data.empty() ? nullptr : global_face_key_data.data(),
          local_counts_i.data(), local_displs_i.data(),
          mpi_type<large_idx_t>(), mesh->comm()->get()));

      std::vector<FaceKey> global_face_keys((size_t)total_faces);
      for (ptrdiff_t i = 0; i < total_faces; ++i) {
        auto &key = global_face_keys[(size_t)i];
        for (int n = 0; n < nnxs; ++n) {
          key[(size_t)n] =
              global_face_key_data[(size_t)i * (size_t)nnxs + (size_t)n];
        }
      }

      const auto face_key_less = [nnxs](const FaceKey &lhs,
                                        const FaceKey &rhs) -> bool {
        return std::lexicographical_compare(lhs.begin(), lhs.begin() + nnxs,
                                            rhs.begin(), rhs.begin() + nnxs);
      };

      std::sort(global_face_keys.begin(), global_face_keys.end(), face_key_less);

      write_pos = 0;
      for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
        // Keep only faces that match a unique owned face globally. If the key
        // is missing or duplicated, the face belongs to another rank or is
        // ambiguous, so we drop it from this local sideset.
        const auto range = std::equal_range(
            global_face_keys.begin(), global_face_keys.end(),
            local_face_keys[(size_t)i], face_key_less);
        if (range.second - range.first != 1) {
          continue;
        }

        parent_element[write_pos] = parent_element[i];
        side_idx[write_pos] = side_idx[i];
        ++write_pos;
      }

      n_surf_elements = write_pos;
    }
#endif
  }

  return std::make_shared<Sideset>(
      mesh->comm(),
      manage_host_buffer<element_idx_t>(n_surf_elements, parent_element),
      manage_host_buffer<i16>(n_surf_elements, side_idx), 0, mesh->comm()->size() > 1 ? mesh->distributed()->element_mapping() : nullptr);
}

#ifdef SMESH_ENABLE_MPI

// FIXME: This is AI Slop code, it should be simplified and optimized
std::shared_ptr<Mesh>
mesh_from_sideset_parallel(const std::shared_ptr<Mesh> &mesh,
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

  const auto parent_dist = mesh->distributed();
  const ptrdiff_t n_parent_owned = parent_dist->n_nodes_owned();
  const ptrdiff_t n_parent_ghosts = parent_dist->n_nodes_ghosts();

  for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
    for (int d = 0; d < nnxs; ++d) {
      idx_t idx = b_surface_elements[d][i];
      if (b_vol2surf[idx] == invalid_idx<idx_t>()) {
        b_vol2surf[idx] = 0;
      }
    }
  }

  ptrdiff_t n_surf_nodes = 0;
  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    if (b_vol2surf[i] == invalid_idx<idx_t>()) {
      continue;
    }

    b_vol2surf[i] = n_surf_nodes++;
  }

  auto local_parent = create_host_buffer<idx_t>(n_surf_nodes);
  auto local_parent_global = create_host_buffer<large_idx_t>(n_surf_nodes);
  auto b_local_parent = local_parent->data();
  auto b_local_parent_global = local_parent_global->data();
  auto b_parent_node_mapping = parent_dist->node_mapping()->data();

  // Compact the sparse volume-to-surface map into dense surface-node arrays while
  // keeping both the parent local index and its globally unique parent id.
  for (ptrdiff_t i = 0; i < n_nodes; ++i) {
    const idx_t local_idx = b_vol2surf[i];
    if (local_idx == invalid_idx<idx_t>()) {
      continue;
    }

    b_local_parent[local_idx] = i;
    b_local_parent_global[local_idx] = b_parent_node_mapping[i];
  }

  std::vector<ptrdiff_t> local_counts(mesh->comm()->size());
  SMESH_MPI_CATCH(MPI_Allgather(&n_surf_nodes, 1, mpi_type<ptrdiff_t>(),
                                local_counts.data(), 1,
                                mpi_type<ptrdiff_t>(), mesh->comm()->get()));

  std::vector<int> local_counts_i(mesh->comm()->size());
  std::vector<int> local_displs_i(mesh->comm()->size());
  ptrdiff_t total_surface_nodes = 0;
  // Build displacements for the allgatherv so every rank can index the flattened
  // list of parent global ids contributed by all other ranks.
  for (int r = 0; r < mesh->comm()->size(); ++r) {
    local_counts_i[r] = static_cast<int>(local_counts[r]);
    local_displs_i[r] = static_cast<int>(total_surface_nodes);
    total_surface_nodes += local_counts[r];
  }

  std::vector<large_idx_t> global_parent_ids((size_t)total_surface_nodes);
  SMESH_MPI_CATCH(MPI_Allgatherv(
      b_local_parent_global, static_cast<int>(n_surf_nodes),
      mpi_type<large_idx_t>(), global_parent_ids.data(), local_counts_i.data(),
      local_displs_i.data(), mpi_type<large_idx_t>(), mesh->comm()->get()));

  std::unordered_map<large_idx_t, int> surface_owner;
  std::unordered_map<large_idx_t, bool> surface_shared;
  surface_owner.reserve(global_parent_ids.size());
  surface_shared.reserve(global_parent_ids.size());
  // The first rank that reports a parent global id becomes its owner; seeing the
  // same id on another rank marks that surface node as shared.
  for (int r = 0; r < mesh->comm()->size(); ++r) {
    const ptrdiff_t begin = local_displs_i[r];
    const ptrdiff_t end = begin + local_counts[r];
    for (ptrdiff_t i = begin; i < end; ++i) {
      const large_idx_t gid = global_parent_ids[(size_t)i];
      const auto inserted = surface_owner.emplace(gid, r);
      if (!inserted.second && inserted.first->second != r) {
        surface_shared[gid] = true;
      }
    }
  }

  std::vector<idx_t> owned_nodes;
  std::vector<idx_t> ghost_nodes;
  std::vector<idx_t> aura_nodes;
  owned_nodes.reserve((size_t)n_surf_nodes);
  ghost_nodes.reserve((size_t)n_surf_nodes);
  aura_nodes.reserve((size_t)n_surf_nodes);

  ptrdiff_t n_surf_shared = 0;
  const int rank = mesh->comm()->rank();
  // Classify each local surface node from the current rank's perspective:
  // owned if this rank won ownership, ghost if the parent node is already a
  // ghost in the volume mesh, otherwise aura.
  for (idx_t i = 0; i < n_surf_nodes; ++i) {
    const large_idx_t gid = b_local_parent_global[i];
    const int owner = surface_owner[gid];
    if (owner == rank) {
      owned_nodes.push_back(i);
      if (surface_shared[gid]) {
        ++n_surf_shared;
      }
    } else if (b_local_parent[i] < n_parent_owned + n_parent_ghosts) {
      ghost_nodes.push_back(i);
    } else {
      aura_nodes.push_back(i);
    }
  }

  auto sort_by_owner = [&](std::vector<idx_t> &nodes) {
    std::stable_sort(nodes.begin(), nodes.end(), [&](const idx_t a, const idx_t b) {
      return surface_owner[b_local_parent_global[a]] <
             surface_owner[b_local_parent_global[b]];
    });
  };
  sort_by_owner(ghost_nodes);
  sort_by_owner(aura_nodes);

  ptrdiff_t n_surf_owned = owned_nodes.size();
  ptrdiff_t n_surf_ghosts = ghost_nodes.size();
  ptrdiff_t n_surf_aura = aura_nodes.size();

  auto surf_points =
      create_host_buffer<geom_t>(mesh->spatial_dimension(), n_surf_nodes);
  auto mapping = create_host_buffer<idx_t>(n_surf_nodes);
  auto old_to_new = create_host_buffer<idx_t>(n_surf_nodes);
  auto surf_node_owner = create_host_buffer<int>(n_surf_nodes);

  auto b_points = mesh->points()->data();
  auto b_surf_points = surf_points->data();
  auto b_mapping = mapping->data();
  auto b_old_to_new = old_to_new->data();
  auto b_surf_node_owner = surf_node_owner->data();

  // Reorder nodes into owned/ghost/aura blocks and carry over geometry, parent
  // mapping, and owning rank in one pass.
  auto assign_nodes = [&](const std::vector<idx_t> &nodes, ptrdiff_t offset) {
    const int spatial_dim = mesh->spatial_dimension();
    for (ptrdiff_t k = 0; k < (ptrdiff_t)nodes.size(); ++k) {
      const idx_t old_idx = nodes[(size_t)k];
      const ptrdiff_t new_idx = offset + k;
      const ptrdiff_t parent_local_idx = b_local_parent[old_idx];
      b_old_to_new[old_idx] = static_cast<idx_t>(new_idx);
      b_mapping[new_idx] = parent_local_idx;
      b_surf_node_owner[new_idx] = surface_owner[b_local_parent_global[old_idx]];
      for (int d = 0; d < spatial_dim; ++d) {
        b_surf_points[d][new_idx] = b_points[d][parent_local_idx];
      }
    }
  };

  assign_nodes(owned_nodes, 0);
  assign_nodes(ghost_nodes, n_surf_owned);
  assign_nodes(aura_nodes, n_surf_owned + n_surf_ghosts);

  for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
    for (int d = 0; d < nnxs; ++d) {
      const idx_t old_idx = b_vol2surf[b_surface_elements[d][i]];
      b_surface_elements[d][i] = b_old_to_new[old_idx];
    }
  }

  auto ret = std::make_shared<Mesh>(mesh->comm(), surface_type,
                                    surface_elements, surf_points);

  ret->set_node_mapping(mapping);

  auto surf_dist = std::make_shared<Distributed>();
  auto surf_node_mapping = create_host_buffer<large_idx_t>(n_surf_nodes);
  auto surf_node_offsets =
      create_host_buffer<ptrdiff_t>(mesh->comm()->size() + 1);
  auto surf_ghosts_and_aura =
      create_host_buffer<idx_t>(n_surf_ghosts + n_surf_aura);

  auto b_surf_node_mapping = surf_node_mapping->data();
  auto b_surf_node_offsets = surf_node_offsets->data();
  auto b_surf_ghosts_and_aura = surf_ghosts_and_aura->data();

  std::vector<ptrdiff_t> owned_counts(mesh->comm()->size());
  SMESH_MPI_CATCH(MPI_Allgather(&n_surf_owned, 1, mpi_type<ptrdiff_t>(),
                                owned_counts.data(), 1,
                                mpi_type<ptrdiff_t>(), mesh->comm()->get()));

  b_surf_node_offsets[0] = 0;
  for (int r = 0; r < mesh->comm()->size(); ++r) {
    b_surf_node_offsets[r + 1] = b_surf_node_offsets[r] + owned_counts[r];
  }

  std::vector<int> owned_counts_i(mesh->comm()->size());
  std::vector<int> owned_displs_i(mesh->comm()->size());
  for (int r = 0; r < mesh->comm()->size(); ++r) {
    owned_counts_i[r] = static_cast<int>(owned_counts[r]);
    owned_displs_i[r] = static_cast<int>(b_surf_node_offsets[r]);
  }

  std::vector<large_idx_t> owned_parent_global_ids(n_surf_owned);
  for (ptrdiff_t i = 0; i < n_surf_owned; ++i) {
    owned_parent_global_ids[i] = b_parent_node_mapping[b_mapping[i]];
    b_surf_node_mapping[i] = b_surf_node_offsets[mesh->comm()->rank()] + i;
  }

  std::vector<large_idx_t> global_owned_parent_ids(
      (size_t)b_surf_node_offsets[mesh->comm()->size()]);
  SMESH_MPI_CATCH(MPI_Allgatherv(
      owned_parent_global_ids.data(), static_cast<int>(n_surf_owned),
      mpi_type<large_idx_t>(), global_owned_parent_ids.data(),
      owned_counts_i.data(), owned_displs_i.data(), mpi_type<large_idx_t>(),
      mesh->comm()->get()));

  std::unordered_map<large_idx_t, large_idx_t> global_surface_node_ids;
  global_surface_node_ids.reserve(global_owned_parent_ids.size());
  for (ptrdiff_t i = 0; i < b_surf_node_offsets[mesh->comm()->size()]; ++i) {
    global_surface_node_ids.emplace(global_owned_parent_ids[(size_t)i], i);
  }

  ptrdiff_t import_idx = 0;
  for (ptrdiff_t i = n_surf_owned; i < n_surf_nodes; ++i) {
    const auto it = global_surface_node_ids.find(b_parent_node_mapping[b_mapping[i]]);
    SMESH_ASSERT(it != global_surface_node_ids.end());
    b_surf_node_mapping[i] = it->second;
    b_surf_ghosts_and_aura[import_idx++] =
        static_cast<idx_t>(b_surf_node_mapping[i]);
  }

  const ptrdiff_t n_parent_owned_elements = parent_dist->n_elements_owned();
  const auto b_parent_surface_elements = sideset->parent()->data();
  ptrdiff_t n_surf_owned_elements = 0;
  ptrdiff_t n_surf_shared_elements = 0;
  for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
    if (b_parent_surface_elements[i] >= n_parent_owned_elements) {
      continue;
    }

    ++n_surf_owned_elements;
    for (int d = 0; d < nnxs; ++d) {
      if (b_surface_elements[d][i] >= n_surf_owned) {
        ++n_surf_shared_elements;
        break;
      }
    }
  }

  ptrdiff_t element_offset = 0;
  ptrdiff_t n_surf_global_elements = 0;
  SMESH_MPI_CATCH(MPI_Exscan(&n_surf_owned_elements, &element_offset, 1,
                             mpi_type<ptrdiff_t>(), MPI_SUM,
                             mesh->comm()->get()));
  if (mesh->comm()->rank() == 0) {
    element_offset = 0;
  }
  SMESH_MPI_CATCH(MPI_Allreduce(&n_surf_owned_elements, &n_surf_global_elements,
                                1, mpi_type<ptrdiff_t>(), MPI_SUM,
                                mesh->comm()->get()));

  auto surf_element_mapping =
      create_host_buffer<large_idx_t>(n_surf_owned_elements);
  auto b_surf_element_mapping = surf_element_mapping->data();
  for (ptrdiff_t i = 0; i < n_surf_owned_elements; ++i) {
    b_surf_element_mapping[i] = element_offset + i;
  }

  surf_dist->impl_->n_nodes_global = b_surf_node_offsets[mesh->comm()->size()];
  surf_dist->impl_->n_nodes_owned = n_surf_owned;
  surf_dist->impl_->n_nodes_shared = n_surf_shared;
  surf_dist->impl_->n_nodes_ghosts = n_surf_ghosts;
  surf_dist->impl_->n_nodes_aura = n_surf_aura;
  surf_dist->impl_->n_elements_global = n_surf_global_elements;
  surf_dist->impl_->n_elements_owned = n_surf_owned_elements;
  surf_dist->impl_->n_elements_shared = n_surf_shared_elements;
  surf_dist->impl_->n_elements_ghosts = n_surf_elements - n_surf_owned_elements;
  surf_dist->impl_->node_mapping = surf_node_mapping;
  surf_dist->impl_->element_mapping = surf_element_mapping;
  surf_dist->impl_->node_owner = surf_node_owner;
  surf_dist->impl_->node_offsets = surf_node_offsets;
  surf_dist->impl_->ghosts_and_aura = surf_ghosts_and_aura;

  ret->impl_->distributed = surf_dist;
  return ret;
}
#endif

std::shared_ptr<Mesh>
mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                  const std::shared_ptr<Sideset> &sideset) {
  SMESH_TRACE_SCOPE("mesh_from_sideset");
#ifdef SMESH_ENABLE_MPI
  if (mesh->comm()->size() > 1) {
    return mesh_from_sideset_parallel(mesh, sideset);
  }
#endif

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

  for (ptrdiff_t i = 0; i < n_surf_elements; ++i) {
    for (int d = 0; d < nnxs; ++d) {
      b_surface_elements[d][i] = b_vol2surf[b_surface_elements[d][i]];
    }
  }

  ret->set_node_mapping(mapping);
  return ret;
}

std::shared_ptr<Mesh> skin(const std::shared_ptr<Mesh> &mesh) {
  auto sideset = skin_sideset(mesh);
  return mesh_from_sideset(mesh, sideset);
}

std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                              const geom_t height, const ptrdiff_t nlayers) {
  if (mesh->n_blocks() != 1) {
    SMESH_ERROR("Extrusion is not supported for multiblock meshes\n");
    return nullptr;
  }
  // This is a hack
  if (mesh->n_nodes_per_element(0) == 4) {
    auto hex8_elements =
        create_host_buffer<idx_t>(8, mesh->n_elements(0) * nlayers);

    auto hex8_points =
        create_host_buffer<geom_t>(3, mesh->n_nodes() * (nlayers + 1));

    quad4_to_hex8_extrude(mesh->n_elements(0), mesh->n_nodes(),
                          mesh->elements(0)->data(), mesh->points()->data(),
                          nlayers, height, hex8_elements->data(),
                          hex8_points->data());

    return std::make_shared<Mesh>(mesh->comm(), HEX8, hex8_elements,
                                  hex8_points);
  } else if (mesh->element_type(0) == TRI3) {
    auto wedge6_elements =
        create_host_buffer<idx_t>(6, mesh->n_elements(0) * nlayers);

    auto wedge6_points =
        create_host_buffer<geom_t>(3, mesh->n_nodes() * (nlayers + 1));

    tri3_to_wedge6_extrude(mesh->n_elements(0), mesh->n_nodes(),
                           mesh->elements(0)->data(), mesh->points()->data(),
                           nlayers, height, wedge6_elements->data(),
                           wedge6_points->data());

    return std::make_shared<Mesh>(mesh->comm(), WEDGE6, wedge6_elements,
                                  wedge6_points);
  } else {
    SMESH_ERROR("Extrusion not supported for element type %s\n",
                type_to_string(mesh->element_type(0)));
    return nullptr;
  }
}

void Mesh::print(std::ostream &os) const {
  os << "n_blocks: " << n_blocks() << "\n";

  for (size_t i = 0; i < n_blocks(); i++) {
    os << i << ")\n";
    block(i)->elements()->print(os);
  }

  os << "n_nodes: " << n_nodes() << "\n";
  points()->print(os);
}
} // namespace smesh
