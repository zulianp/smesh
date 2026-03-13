#ifndef SMESH_MESH_HPP
#define SMESH_MESH_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_crs_graph.hpp"
#include "smesh_forward_declarations.hpp"

// STL
#include <functional>

namespace smesh {

class Distributed {
public:
  Distributed();
  ~Distributed();

  ptrdiff_t n_nodes_global() const;
  ptrdiff_t n_elements_global() const;

  ptrdiff_t n_nodes_local() const;
  ptrdiff_t n_nodes_owned_not_shared() const;
  ptrdiff_t n_nodes_owned() const;
  ptrdiff_t n_nodes_shared() const;
  ptrdiff_t n_nodes_ghosts() const;
  ptrdiff_t n_nodes_aura() const;

  ptrdiff_t n_elements_local() const;
  ptrdiff_t n_elements_owned_not_shared() const;
  ptrdiff_t n_elements_owned() const;
  ptrdiff_t n_elements_shared() const;
  ptrdiff_t n_elements_ghosts() const;

  SharedBuffer<large_idx_t> node_mapping() const;
  SharedBuffer<large_idx_t> element_mapping() const;
  SharedBuffer<large_idx_t> aura_element_mapping() const;

  SharedBuffer<int> node_owner() const;
  SharedBuffer<ptrdiff_t> node_offsets() const;
  SharedBuffer<idx_t> ghosts() const;
  SharedBuffer<idx_t> ghosts_and_aura() const;

  friend class Mesh;
  friend std::shared_ptr<Mesh>
  mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                    const std::shared_ptr<Sideset> &sideset);
  friend std::shared_ptr<Mesh>
  mesh_from_sideset_parallel(const std::shared_ptr<Mesh> &mesh,
                             const std::shared_ptr<Sideset> &sideset);

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

enum ComputeDataFlags {
  DD_NONE = 0,
  DD_ELEMENT_SOA = 1 << 0,
  DD_ELEMENT_AOS = 1 << 1,
  DD_POINT_AOS = 1 << 2,
  DD_POINT_SOA = 1 << 3,
  DD_JACOBIAN_SOA = 1 << 0,
  DD_JACOBIAN_AOS = 1 << 1,
  DD_JACOBIAN_ADJUGATE_SOA = 1 << 2,
  DD_JACOBIAN_ADJUGATE_AOS = 1 << 3,
  DD_JACOBIAN_DETERMINANT = 1 << 4,
  DD_ALL_SOA = DD_ELEMENT_SOA | DD_POINT_SOA | DD_JACOBIAN_SOA |
               DD_JACOBIAN_ADJUGATE_SOA | DD_JACOBIAN_DETERMINANT,
  DD_ALL_AOS = DD_ELEMENT_AOS | DD_POINT_AOS | DD_JACOBIAN_AOS |
               DD_JACOBIAN_ADJUGATE_AOS,
  DD_ALL = DD_ALL_SOA | DD_ALL_AOS
};

class ComputeData {
public:
  ComputeData();
  ~ComputeData();
  SharedBuffer<idx_t *> elements_SoA(const block_idx_t block_id);
  SharedBuffer<idx_t> elements_AoS(const block_idx_t block_id);

  SharedBuffer<geom_t *> points_SoA();
  SharedBuffer<geom_t> points_AoS();

  // Precision may vary based on compilation flags
  SharedBuffer<void *> jacobians_SoA(const block_idx_t block_id);
  SharedBuffer<void> jacobians_AoS(const block_idx_t block_id);
  SharedBuffer<void *> jacobian_adjugate_SoA(const block_idx_t block_id);
  SharedBuffer<void> jacobian_adjugate_AoS(const block_idx_t block_id);
  SharedBuffer<void> jacobian_determinant(const block_idx_t block_id);

  void set_num_blocks(const ptrdiff_t num_blocks);

  void set_elements_SoA(const block_idx_t block_id,
                        const SharedBuffer<idx_t *> &elements);
  void set_elements_AoS(const block_idx_t block_id,
                        const SharedBuffer<idx_t> &elements);
  void set_points_SoA(const SharedBuffer<geom_t *> &points);
  void set_points_AoS(const SharedBuffer<geom_t> &points);
  void set_jacobians_SoA(const block_idx_t block_id,
                         const SharedBuffer<void *> &jacobians);
  void set_jacobians_AoS(const block_idx_t block_id,
                         const SharedBuffer<void> &jacobians);
  void set_jacobian_adjugate_SoA(const block_idx_t block_id,
                                 const SharedBuffer<void *> &jacobian_adjugate);
  void set_jacobian_adjugate_AoS(const block_idx_t block_id,
                                 const SharedBuffer<void> &jacobian_adjugate);
  void set_jacobian_determinant(const block_idx_t block_id,
                                const SharedBuffer<void> &jacobian_determinant);

  int commit_to_device();

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

class Mesh final {
public:
  using NodeToNodeGraph = smesh::CRSGraph<count_t, idx_t>;
  using NodeToElementGraph = smesh::CRSGraph<count_t, element_idx_t>;

  class Block {
  public:
    Block();
    ~Block();
    Block(const std::string &name, enum ElemType element_type,
          SharedBuffer<idx_t *> elements);

    const std::string &name() const;
    enum ElemType element_type() const;
    int n_nodes_per_element() const;
    const SharedBuffer<idx_t *> &elements() const;

    void set_name(const std::string &name);
    void set_element_type(enum ElemType element_type);
    void set_elements(SharedBuffer<idx_t *> elements);
    ptrdiff_t n_elements() const;

  private:
    class Impl;
    std::unique_ptr<Impl> impl_;
  };

  Mesh();
  Mesh(const std::shared_ptr<Communicator> &comm);
  ~Mesh();

  Mesh(const std::shared_ptr<Communicator> &comm, enum ElemType element_type,
       SharedBuffer<idx_t *> elements, SharedBuffer<geom_t *> points);

  Mesh(const std::shared_ptr<Communicator> &comm,
       const std::vector<std::shared_ptr<Block>> &blocks,
       SharedBuffer<geom_t *> points);

  int init_compute_data(const int flags, const enum ExecutionSpace space);
  std::shared_ptr<ComputeData> compute_data() const;

  int read(const Path &path);
  int write(const Path &path) const;
  int initialize_node_to_node_graph();
  int convert_to_macro_element_mesh();
  const std::vector<std::shared_ptr<Block>> &blocks() const;
  std::vector<std::shared_ptr<Block>>
  blocks(const std::vector<std::string> &block_names) const;

  // Block-related methods
  size_t n_blocks() const;
  std::shared_ptr<const Block> block(size_t index) const;
  std::shared_ptr<Block> block(size_t index);
  std::shared_ptr<Block> find_block(const std::string &name) const;
  void add_block(const std::string &name, enum ElemType element_type,
                 SharedBuffer<idx_t *> elements);
  void add_block(const std::shared_ptr<Block> &block);
  void remove_block(size_t index);

  std::shared_ptr<Distributed> distributed() const;

  int spatial_dimension() const;
  ptrdiff_t n_nodes() const;
  ptrdiff_t n_elements() const;

  int n_nodes_per_element(block_idx_t block_id) const;
  ptrdiff_t n_elements(block_idx_t block_id) const;
  enum ElemType element_type(block_idx_t block_id) const;
  SharedBuffer<idx_t *> elements(block_idx_t block_id);
  SharedBuffer<idx_t *> elements(block_idx_t block_id) const;

  std::shared_ptr<NodeToNodeGraph> node_to_node_graph();
  std::shared_ptr<NodeToNodeGraph> node_to_node_graph_upper_triangular();
  std::shared_ptr<NodeToElementGraph> node_to_element_graph();
  SharedBuffer<element_idx_t> half_face_table();
  std::shared_ptr<NodeToNodeGraph>
  create_node_to_node_graph(const enum ElemType element_type);

  SharedBuffer<count_t> node_to_node_rowptr() const;
  SharedBuffer<idx_t> node_to_node_colidx() const;
  SharedBuffer<ptrdiff_t> node_offsets() const;
  SharedBuffer<idx_t> ghosts() const;
  SharedBuffer<int> node_owner() const;
  SharedBuffer<idx_t> node_mapping() const;

  SharedBuffer<geom_t *> points();
  SharedBuffer<geom_t *> points() const;

  void set_points(const SharedBuffer<geom_t *> &points);

  std::shared_ptr<Communicator> comm() const;

  inline static std::shared_ptr<Mesh>
  create_from_file(const std::shared_ptr<Communicator> &comm,
                   const Path &path) {
    auto ret = std::make_shared<Mesh>(comm);
    ret->read(path);
    return ret;
  }

  static std::shared_ptr<Mesh> create_hex8_reference_cube();

  static std::shared_ptr<Mesh>
  create_cube(const std::shared_ptr<Communicator> &comm,
              const enum ElemType element_type, const ptrdiff_t nx = 1,
              const ptrdiff_t ny = 1, const ptrdiff_t nz = 1,
              const geom_t xmin = 0, const geom_t ymin = 0,
              const geom_t zmin = 0, const geom_t xmax = 1,
              const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh> create_hex8_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 1,
      const ptrdiff_t ny = 1, const ptrdiff_t nz = 1, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh> create_semistructured_hex_cube(
      const std::shared_ptr<Communicator> &comm,
      const int micro_elements_per_dim = 2, const ptrdiff_t nx = 1,
      const ptrdiff_t ny = 1, const ptrdiff_t nz = 1, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh> create_tet4_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 1,
      const ptrdiff_t ny = 1, const ptrdiff_t nz = 1, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh>
  create_square(const std::shared_ptr<Communicator> &comm,
                const enum ElemType element_type, const ptrdiff_t nx = 1,
                const ptrdiff_t ny = 1, const geom_t xmin = 0,
                const geom_t ymin = 0, const geom_t xmax = 1,
                const geom_t ymax = 1);

  static std::shared_ptr<Mesh>
  create_tri3_square(const std::shared_ptr<Communicator> &comm,
                     const ptrdiff_t nx = 1, const ptrdiff_t ny = 1,
                     const geom_t xmin = 0, const geom_t ymin = 0,
                     const geom_t xmax = 1, const geom_t ymax = 1);

  static std::shared_ptr<Mesh>
  create_quad4_square(const std::shared_ptr<Communicator> &comm,
                      const ptrdiff_t nx = 1, const ptrdiff_t ny = 1,
                      const geom_t xmin = 0, const geom_t ymin = 0,
                      const geom_t xmax = 1, const geom_t ymax = 1);

  static std::shared_ptr<Mesh>
  create_quad4_ring(const std::shared_ptr<Communicator> &comm,
                    const geom_t inner_radius, const geom_t outer_radius,
                    const ptrdiff_t nlayers, const ptrdiff_t nelements);

  static std::shared_ptr<Mesh> create_hex8_checkerboard_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 2,
      const ptrdiff_t ny = 2, const ptrdiff_t nz = 2, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh> create_hex8_bidomain_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 2,
      const ptrdiff_t ny = 2, const ptrdiff_t nz = 2, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  std::vector<std::pair<block_idx_t, SharedBuffer<element_idx_t>>>
  select_elements(const std::function<bool(const geom_t, const geom_t,
                                           const geom_t)> &selector,
                  const std::vector<std::string> &block_names = {});

  int split_block(const SharedBuffer<element_idx_t> &elements,
                  const std::string &name);
  int split_boundary_layer();
  int renumber_nodes();
  int renumber_nodes(const SharedBuffer<idx_t> &node_mapping);
  void set_node_mapping(const SharedBuffer<idx_t> &node_mapping);
  void set_comm(const std::shared_ptr<Communicator> &comm);
  void set_element_type(const block_idx_t block_id,
                        const enum ElemType element_type);
  std::pair<SharedBuffer<geom_t>, SharedBuffer<geom_t>> compute_bounding_box();

  std::shared_ptr<Mesh> clone() const;

  void reorder_elements_from_tags(const block_idx_t block_id,
                                  const SharedBuffer<idx_t> &tags);

  void print(std::ostream &os = std::cout) const;

private:
  friend std::shared_ptr<Mesh>
  mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                    const std::shared_ptr<Sideset> &sideset);
  friend std::shared_ptr<Mesh>
  mesh_from_sideset_parallel(const std::shared_ptr<Mesh> &mesh,
                             const std::shared_ptr<Sideset> &sideset);

  class Impl;
  std::unique_ptr<Impl> impl_;
};

using SharedMesh = std::shared_ptr<Mesh>;
using SharedBlock = std::shared_ptr<Mesh::Block>;

std::shared_ptr<Mesh> convert_to(const enum ElemType element_type,
                                 const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh> promote_to(const enum ElemType element_type,
                                 const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                             const int levels = 1);
std::shared_ptr<Sideset> skin_sideset(const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh>
mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                  const std::shared_ptr<Sideset> &sideset);
std::shared_ptr<Mesh>
mesh_from_sideset_parallel(const std::shared_ptr<Mesh> &mesh,
                           const std::shared_ptr<Sideset> &sideset);
std::shared_ptr<Mesh> skin(const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                              const geom_t height, const ptrdiff_t nlayers);

} // namespace smesh

#endif // SMESH_MESH_HPP
