#ifndef SMESH_MESH_HPP
#define SMESH_MESH_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_crs_graph.hpp"
#include "smesh_forward_declarations.hpp"

// STL
#include <functional>
#include <vector>

namespace smesh {

class MeshTransformsDistributed;

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
  friend std::shared_ptr<Mesh>
  to_semistructured(const int level, const std::shared_ptr<Mesh> &mesh,
                    const bool hiearchical_ordering, const bool use_GLL);
  friend std::shared_ptr<Mesh> to_semistructured_distributed(
      const int level, const std::shared_ptr<Mesh> &mesh, const bool use_GLL,
      const bool hierarchical_ordering);
  friend std::shared_ptr<Mesh> to_semistructured_distributed_hex_tet(
      const int level, const std::shared_ptr<Mesh> &mesh,
      const bool hierarchical_ordering);
  friend std::shared_ptr<Mesh> to_semistructured_distributed_quad(
      const int level, const std::shared_ptr<Mesh> &mesh,
      const bool hierarchical_ordering);
  friend class MeshTransformsDistributed;
  friend std::shared_ptr<Mesh> sshex_to_hex8(const std::shared_ptr<Mesh> &sshex);
  friend std::shared_ptr<Mesh> derefine(const std::shared_ptr<Mesh> &mesh,
                                        const int to_level);
  friend std::shared_ptr<Mesh> convert_to(const enum ElemType element_type,
                                          const std::shared_ptr<Mesh> &mesh);
  friend std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                                      const int levels);
  friend std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                                       const geom_t height,
                                       const ptrdiff_t nlayers);

private:
  void set_nodes(ptrdiff_t n_global, ptrdiff_t n_owned, ptrdiff_t n_shared,
                 ptrdiff_t n_ghosts, ptrdiff_t n_aura,
                 SharedBuffer<large_idx_t> node_mapping,
                 SharedBuffer<int> node_owner,
                 SharedBuffer<ptrdiff_t> node_offsets,
                 SharedBuffer<idx_t> ghosts_and_aura);
  void set_elements(ptrdiff_t n_global, ptrdiff_t n_owned, ptrdiff_t n_shared,
                    ptrdiff_t n_ghosts,
                    SharedBuffer<large_idx_t> element_mapping,
                    SharedBuffer<large_idx_t> aura_element_mapping);

  class Impl;
  std::unique_ptr<Impl> impl_;
};

/// Per-block owned/shared/aura element layout. Nodes stay mesh-global on
/// `Distributed`. SoA order is `[owned-not-shared | shared | aura]`.
class DistributedBlock {
public:
  DistributedBlock();
  ~DistributedBlock();

  ptrdiff_t n_elements_local() const;
  ptrdiff_t n_elements_owned_not_shared() const;
  ptrdiff_t n_elements_owned() const;
  ptrdiff_t n_elements_shared() const;
  ptrdiff_t n_elements_ghosts() const;

  SharedBuffer<large_idx_t> element_mapping() const;
  SharedBuffer<large_idx_t> aura_element_mapping() const;

  void set_elements(ptrdiff_t n_owned, ptrdiff_t n_shared, ptrdiff_t n_ghosts,
                    SharedBuffer<large_idx_t> element_mapping,
                    SharedBuffer<large_idx_t> aura_element_mapping);

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

// enum KernelDataFlags {
//     GEO_NONE                  = 0,
//     GEO_ELEMENT_SOA           = 1 << 0,
//     GEO_ELEMENT_AOS           = 1 << 1,
//     GEO_POINT_AOS             = 1 << 2,
//     GEO_POINT_SOA             = 1 << 3,
//     GEO_JACOBIAN_SOA          = 1 << 0,
//     GEO_JACOBIAN_AOS          = 1 << 1,
//     GEO_JACOBIAN_ADJUGATE_SOA = 1 << 2,
//     GEO_JACOBIAN_ADJUGATE_AOS = 1 << 3,
//     GEO_JACOBIAN_DETERMINANT  = 1 << 4,
//     GEO_ALL_SOA = GEO_ELEMENT_SOA | GEO_POINT_SOA | GEO_JACOBIAN_SOA |
//     GEO_JACOBIAN_ADJUGATE_SOA | GEO_JACOBIAN_DETERMINANT, GEO_ALL_AOS =
//     GEO_ELEMENT_AOS | GEO_POINT_AOS | GEO_JACOBIAN_AOS |
//     GEO_JACOBIAN_ADJUGATE_AOS, GEO_ALL     = GEO_ALL_SOA | GEO_ALL_AOS
// };

// class KernelData {
// public:
//     KernelData();
//     ~KernelData();
//     SharedBuffer<idx_t *> elements_SoA(const block_idx_t block_id);
//     SharedBuffer<idx_t>   elements_AoS(const block_idx_t block_id);

//     SharedBuffer<geom_t *> points_SoA();
//     SharedBuffer<geom_t>   points_AoS();

//     // Precision may vary based on compilation flags
//     SharedBuffer<jacobian_t *> jacobians_SoA(const block_idx_t block_id);
//     SharedBuffer<jacobian_t>   jacobians_AoS(const block_idx_t block_id);
//     SharedBuffer<jacobian_t *> jacobian_adjugate_SoA(const block_idx_t
//     block_id); SharedBuffer<jacobian_t>   jacobian_adjugate_AoS(const
//     block_idx_t block_id); SharedBuffer<geom_t> jacobian_determinant(const
//     block_idx_t block_id);

//     void set_num_blocks(const ptrdiff_t num_blocks);

//     void set_elements_SoA(const block_idx_t block_id, const
//     SharedBuffer<idx_t *> &elements); void set_elements_AoS(const block_idx_t
//     block_id, const SharedBuffer<idx_t> &elements); void set_points_SoA(const
//     SharedBuffer<geom_t *> &points); void set_points_AoS(const
//     SharedBuffer<geom_t> &points); void set_jacobians_SoA(const block_idx_t
//     block_id, const SharedBuffer<jacobian_t *> &jacobians); void
//     set_jacobians_AoS(const block_idx_t block_id, const
//     SharedBuffer<jacobian_t> &jacobians); void
//     set_jacobian_adjugate_SoA(const block_idx_t block_id, const
//     SharedBuffer<jacobian_t *> &jacobian_adjugate); void
//     set_jacobian_adjugate_AoS(const block_idx_t block_id, const
//     SharedBuffer<jacobian_t> &jacobian_adjugate); void
//     set_jacobian_determinant(const block_idx_t block_id, const
//     SharedBuffer<geom_t> &jacobian_determinant);

//     int send_to_device();

// private:
//     class Impl;
//     std::unique_ptr<Impl> impl_;
// };

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

    std::shared_ptr<DistributedBlock> distributed() const;
    void set_distributed(const std::shared_ptr<DistributedBlock> &distributed);

    ptrdiff_t n_elements_owned() const;
    ptrdiff_t n_elements_shared() const;
    ptrdiff_t n_elements_ghosts() const;
    ptrdiff_t n_elements_owned_not_shared() const;
    SharedBuffer<large_idx_t> element_mapping() const;
    SharedBuffer<large_idx_t> aura_element_mapping() const;
    void set_distributed_elements(ptrdiff_t n_owned, ptrdiff_t n_shared,
                                  ptrdiff_t n_ghosts,
                                  SharedBuffer<large_idx_t> element_mapping,
                                  SharedBuffer<large_idx_t> aura_element_mapping);

    SharedBuffer<idx_t *> device_elements_SoA();
    SharedBuffer<idx_t> device_elements_AoS();
    void set_device_elements_SoA(const SharedBuffer<idx_t *> &elements);

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

  // std::shared_ptr<KernelData> create_kernel_data(const int flags, const enum
  // ExecutionSpace space);

  // int                         init_kernel_data(const int flags, const enum
  // ExecutionSpace space); std::shared_ptr<KernelData> kernel_data() const;

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
  bool is_distributed() const;

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
  std::shared_ptr<NodeToNodeGraph> edge_graph();
  std::shared_ptr<NodeToElementGraph> node_to_element_graph();
  SharedBuffer<block_idx_t> node_to_element_block_number() const;
  SharedBuffer<element_idx_t> half_face_table();
  SharedBuffer<element_idx_t> half_face_table(block_idx_t block_id);
  SharedBuffer<block_idx_t> half_face_neighbor_block(block_idx_t block_id);
  std::shared_ptr<NodeToNodeGraph>
  create_node_to_node_graph(const enum ElemType element_type);

  // SharedBuffer<count_t> node_to_node_rowptr() const;
  // SharedBuffer<idx_t> node_to_node_colidx() const;
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

#ifdef SMESH_ENABLE_RYAML
  static std::shared_ptr<Mesh>
  create_from_yaml(const std::shared_ptr<Communicator> &comm,
                   const ryml::NodeRef &node);
#endif

  static std::shared_ptr<Mesh> create_hex8_reference_cube();

  static std::shared_ptr<Mesh>
  create_cube(const std::shared_ptr<Communicator> &comm,
              const enum ElemType element_type, const ptrdiff_t nx = 1,
              const ptrdiff_t ny = 1, const ptrdiff_t nz = 1,
              const geom_t xmin = 0, const geom_t ymin = 0,
              const geom_t zmin = 0, const geom_t xmax = 1,
              const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh> create_wall_mounted_hump(
      const std::shared_ptr<Communicator> &comm,
      const enum ElemType element_type, const ptrdiff_t nx = 32,
      const ptrdiff_t ny = 12, const ptrdiff_t nz = 4,
      const geom_t length = 9.0, const geom_t height = 3.0,
      const geom_t width = 1.0, const geom_t hump_start = 0.65,
      const geom_t hump_length = 1.0, const geom_t hump_height = 0.128);

  static std::shared_ptr<Mesh>
  create_half_sphere(const std::shared_ptr<Communicator> &comm,
                     const enum ElemType element_type, const geom_t radius,
                     const ptrdiff_t nx = 8, const ptrdiff_t ny = 8,
                     const ptrdiff_t nz = 4);

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

  static std::shared_ptr<Mesh>
  create_tet4_half_sphere(const std::shared_ptr<Communicator> &comm,
                          const geom_t radius, const ptrdiff_t nx = 8,
                          const ptrdiff_t ny = 8, const ptrdiff_t nz = 4);

  static std::shared_ptr<Mesh>
  create_hex8_half_sphere(const std::shared_ptr<Communicator> &comm,
                          const geom_t radius, const ptrdiff_t nx = 8,
                          const ptrdiff_t ny = 8, const ptrdiff_t nz = 4);

  static std::shared_ptr<Mesh> create_hex8_checkerboard_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 2,
      const ptrdiff_t ny = 2, const ptrdiff_t nz = 2, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  /// HEX8 cube with the second half of hexes converted to 6 TET4 each.
  /// Blocks `"hex"` / `"tet"` share unique nodes.
  static std::shared_ptr<Mesh> create_hex8_tet4_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 2,
      const ptrdiff_t ny = 2, const ptrdiff_t nz = 2, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  static std::shared_ptr<Mesh> create_hex8_bidomain_cube(
      const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx = 2,
      const ptrdiff_t ny = 2, const ptrdiff_t nz = 2, const geom_t xmin = 0,
      const geom_t ymin = 0, const geom_t zmin = 0, const geom_t xmax = 1,
      const geom_t ymax = 1, const geom_t zmax = 1);

  /// Serial HEX-dominant unit: HEX8 + PYRAMID5 + WEDGE6 + TET4, shared nodes.
  static std::shared_ptr<Mesh>
  create_hex_dominant_serial(const std::shared_ptr<Communicator> &comm);

  /// HEX-dominant solid cylinder (axis +z, centerline at x=y=0).
  /// HEX8 Cartesian core (nodes on polar rays) + polar HEX8 annulus;
  /// WEDGE6 only for leftover mid-side polar nodes, not at square corners.
  /// Requires ntheta >= 8 and ntheta % 4 == 0.
  static std::shared_ptr<Mesh> create_hex_dominant_cylinder(
      const std::shared_ptr<Communicator> &comm, const geom_t radius,
      const geom_t height, const ptrdiff_t nr = 2, const ptrdiff_t ntheta = 16,
      const ptrdiff_t nz = 4, const geom_t zmin = 0);

  std::vector<std::pair<block_idx_t, SharedBuffer<element_idx_t>>>
  select_elements(const std::function<bool(const geom_t, const geom_t,
                                           const geom_t)> &selector,
                  const std::vector<std::string> &block_names = {});

  int split_block(const SharedBuffer<element_idx_t> &elements,
                  const std::string &name, block_idx_t block_id = 0);
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

  SharedBuffer<geom_t *> device_points_SoA();
  SharedBuffer<geom_t> device_points_AoS();

private:
  friend std::shared_ptr<Mesh>
  mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                    const std::shared_ptr<Sideset> &sideset);
  friend std::shared_ptr<Mesh>
  mesh_from_sideset_parallel(const std::shared_ptr<Mesh> &mesh,
                             const std::shared_ptr<Sideset> &sideset);
  friend std::shared_ptr<Mesh>
  to_semistructured(const int level, const std::shared_ptr<Mesh> &mesh,
                    const bool hiearchical_ordering, const bool use_GLL);
  friend std::shared_ptr<Mesh> to_semistructured_distributed(
      const int level, const std::shared_ptr<Mesh> &mesh, const bool use_GLL,
      const bool hierarchical_ordering);
  friend std::shared_ptr<Mesh> to_semistructured_distributed_hex_tet(
      const int level, const std::shared_ptr<Mesh> &mesh,
      const bool hierarchical_ordering);
  friend std::shared_ptr<Mesh> to_semistructured_distributed_quad(
      const int level, const std::shared_ptr<Mesh> &mesh,
      const bool hierarchical_ordering);
  friend class MeshTransformsDistributed;
  friend std::shared_ptr<Mesh> sshex_to_hex8(const std::shared_ptr<Mesh> &sshex);
  friend std::shared_ptr<Mesh> derefine(const std::shared_ptr<Mesh> &mesh,
                                        const int to_level);
  friend std::shared_ptr<Mesh> convert_to(const enum ElemType element_type,
                                          const std::shared_ptr<Mesh> &mesh);
  friend std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                                      const int levels);
  friend std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                                       const geom_t height,
                                       const ptrdiff_t nlayers);

  void set_distributed(const std::shared_ptr<Distributed> &distributed);

#ifdef SMESH_ENABLE_MPI
  static int adopt_parallel_arrays(Mesh *mesh, enum ElemType element_type,
                                   const char *block_name, int nnodesxelem,
                                   ptrdiff_t n_global_elements,
                                   ptrdiff_t n_owned_elements,
                                   ptrdiff_t n_shared_elements,
                                   ptrdiff_t n_ghost_elements,
                                   large_idx_t *element_mapping,
                                   large_idx_t *aura_element_mapping,
                                   idx_t **elements, int spatial_dim,
                                   ptrdiff_t n_global_nodes,
                                   ptrdiff_t n_owned_nodes,
                                   ptrdiff_t n_shared_nodes,
                                   ptrdiff_t n_ghost_nodes,
                                   ptrdiff_t n_aura_nodes,
                                   large_idx_t *node_mapping, geom_t **points,
                                   int *node_owner, ptrdiff_t *node_offsets,
                                   idx_t *ghosts);
  static std::shared_ptr<Mesh>
  wrap_create_parallel(const std::shared_ptr<Communicator> &comm,
                       enum ElemType element_type, int nnodesxelem,
                       ptrdiff_t n_local_elements, ptrdiff_t n_global_elements,
                       idx_t **elems, int spatial_dim, ptrdiff_t n_local_nodes,
                       ptrdiff_t n_global_nodes, geom_t **points);
  static std::shared_ptr<Mesh> with_nodal_distributed(
      const std::shared_ptr<Mesh> &src,
      const std::vector<std::shared_ptr<Block>> &blocks,
      ptrdiff_t n_elements_global, ptrdiff_t n_owned, ptrdiff_t n_shared,
      ptrdiff_t n_ghosts, SharedBuffer<large_idx_t> element_mapping,
      SharedBuffer<large_idx_t> aura_element_mapping);
  static std::shared_ptr<Mesh>
  split_hex8_checkerboard_distributed(const std::shared_ptr<Mesh> &hex_mesh,
                                      const ptrdiff_t nx, const ptrdiff_t ny);
  static std::shared_ptr<Mesh>
  split_hex8_tet4_distributed(const std::shared_ptr<Mesh> &hex_mesh,
                              const ptrdiff_t n_hex_all);
#endif

  class Impl;
  std::unique_ptr<Impl> impl_;
};

using SharedMesh = std::shared_ptr<Mesh>;
using SharedBlock = std::shared_ptr<Mesh::Block>;

#ifdef SMESH_ENABLE_MPI
class MeshTransformsDistributed {
public:
  static std::shared_ptr<Distributed> copy_distributed(const Distributed &src);
  static std::shared_ptr<Distributed>
  make_nodal_distributed(const std::shared_ptr<Communicator> &comm,
                         ptrdiff_t n_global, ptrdiff_t n_owned, ptrdiff_t n_shared,
                         ptrdiff_t n_ghosts, ptrdiff_t n_aura,
                         SharedBuffer<large_idx_t> node_mapping,
                         SharedBuffer<int> node_owner, ptrdiff_t n_elem_global,
                         ptrdiff_t n_elem_owned, ptrdiff_t n_elem_shared,
                         ptrdiff_t n_elem_ghosts,
                         SharedBuffer<large_idx_t> element_mapping,
                         SharedBuffer<large_idx_t> aura_element_mapping);
  static std::shared_ptr<Mesh>
  make_distributed_mesh(const std::shared_ptr<Communicator> &comm,
                        const std::vector<std::shared_ptr<Mesh::Block>> &blocks,
                        const SharedBuffer<geom_t *> &points,
                        const std::shared_ptr<Distributed> &dist);
  static void expand_element_maps(Distributed &dist, const int factor);
  static void expand_block_elements(Mesh::Block &block, const int factor);
  static int conversion_factor(const enum ElemType from, const enum ElemType to);
  static std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                                      const int levels);
  static std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                                       const geom_t height,
                                       const ptrdiff_t nlayers);
  static void clone_distributed(const Mesh &src, Mesh &dst);
  static int attach_convert_distributed(const Mesh &src, Mesh &dst);
  static std::shared_ptr<Mesh>
  attach_sshex_to_hex8(const std::shared_ptr<Mesh> &ss,
                       const std::shared_ptr<Mesh> &hex);
  static std::shared_ptr<Mesh>
  derefine(const std::shared_ptr<Mesh> &mesh,
           std::vector<std::shared_ptr<Mesh::Block>> &blocks);
};
#endif

std::shared_ptr<Mesh> convert_to(const enum ElemType element_type,
                                 const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh> promote_to(const enum ElemType element_type,
                                 const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh> refine(const std::shared_ptr<Mesh> &mesh,
                             const int levels = 1);
std::shared_ptr<Sideset> skin_sideset(const std::shared_ptr<Mesh> &mesh);
std::vector<std::shared_ptr<Sideset>>
skin_sidesets(const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh>
mesh_from_sideset(const std::shared_ptr<Mesh> &mesh,
                  const std::shared_ptr<Sideset> &sideset);
std::shared_ptr<Mesh>
mesh_from_sideset_parallel(const std::shared_ptr<Mesh> &mesh,
                           const std::shared_ptr<Sideset> &sideset);
std::shared_ptr<Mesh> skin(const std::shared_ptr<Mesh> &mesh);
std::shared_ptr<Mesh> extrude(const std::shared_ptr<Mesh> &mesh,
                              const geom_t height, const ptrdiff_t nlayers);

std::shared_ptr<Mesh> concatenate(const std::shared_ptr<Mesh> &mesh1,
                                  const std::shared_ptr<Mesh> &mesh2);

bool surface_is_closed(const std::shared_ptr<Mesh> &mesh);

} // namespace smesh

#endif // SMESH_MESH_HPP
