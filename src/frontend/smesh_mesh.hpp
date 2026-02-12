#ifndef SMESH_MESH_HPP
#define SMESH_MESH_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_crs_graph.hpp"

// STL
#include <functional>

namespace smesh {
    class Mesh final {
    public:
        using CRSGraph = smesh::CRSGraph<count_t, idx_t>;

        class Block {
        public:
            Block();
            ~Block();

            const std::string           &name() const;
            enum ElemType                element_type() const;
            int                          n_nodes_per_element() const;
            const SharedBuffer<idx_t *> &elements() const;

            // Setters for internal use
            void      set_name(const std::string &name);
            void      set_element_type(enum ElemType element_type);
            void      set_elements(SharedBuffer<idx_t *> elements);
            ptrdiff_t n_elements() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

        Mesh();
        Mesh(const std::shared_ptr<Communicator> &comm);
        ~Mesh();

        Mesh(const std::shared_ptr<Communicator> &comm,
             enum ElemType                        element_type,
             SharedBuffer<idx_t *>                elements,
             SharedBuffer<geom_t *>               points);

        Mesh(const std::shared_ptr<Communicator> &comm,
             const std::vector<std::shared_ptr<Block>> &blocks,
             SharedBuffer<geom_t *> points);

        friend class FunctionSpace;
        friend class Op;
        // friend class NeumannConditions;

        int                                        read(const Path &path);
        int                                        write(const Path &path) const;
        int                                        initialize_node_to_node_graph();
        int                                        convert_to_macro_element_mesh();
        const std::vector<std::shared_ptr<Block>> &blocks() const;
        std::vector<std::shared_ptr<Block>>        blocks(const std::vector<std::string> &block_names) const;

        // Block-related methods
        size_t                       n_blocks() const;
        std::shared_ptr<const Block> block(size_t index) const;
        std::shared_ptr<Block>       block(size_t index);
        std::shared_ptr<Block>       find_block(const std::string &name) const;
        void add_block(const std::string &name, enum ElemType element_type, SharedBuffer<idx_t *> elements);
        void remove_block(size_t index);

        int           spatial_dimension() const;
        int           n_nodes_per_element() const;
        ptrdiff_t     n_nodes() const;
        ptrdiff_t     n_elements() const;
        enum ElemType element_type() const;
        ptrdiff_t     n_owned_nodes() const;
        ptrdiff_t     n_owned_nodes_with_ghosts() const;
        ptrdiff_t     n_owned_elements() const;
        ptrdiff_t     n_owned_elements_with_ghosts() const;
        ptrdiff_t     n_shared_elements() const;

        std::shared_ptr<CRSGraph>   node_to_node_graph();
        std::shared_ptr<CRSGraph>   node_to_node_graph_upper_triangular();
        SharedBuffer<element_idx_t> half_face_table();
        std::shared_ptr<CRSGraph>   create_node_to_node_graph(const enum ElemType element_type);

        SharedBuffer<count_t> node_to_node_rowptr() const;
        SharedBuffer<idx_t>   node_to_node_colidx() const;
        SharedBuffer<idx_t>   node_offsets() const;
        SharedBuffer<idx_t>   ghosts() const;
        SharedBuffer<int>     node_owner() const;
        SharedBuffer<idx_t>   node_mapping() const;
        SharedBuffer<idx_t>   element_mapping() const;

        const geom_t * points(const int coord) const;
        const idx_t *  idx(const int node_num) const;

        SharedBuffer<geom_t *> points();
        SharedBuffer<idx_t *>  elements();
        SharedBuffer<idx_t *>  default_elements();  // For backward compatibility

        std::shared_ptr<Communicator> comm() const;

        inline static std::shared_ptr<Mesh> create_from_file(const std::shared_ptr<Communicator> &comm, const Path &path) {
            auto ret = std::make_shared<Mesh>(comm);
            ret->read(path);
            return ret;
        }

        static std::shared_ptr<Mesh> create_hex8_reference_cube();

        static std::shared_ptr<Mesh> create_cube(const std::shared_ptr<Communicator> &comm,
                                                 const enum ElemType                  element_type,
                                                 const int                            nx   = 1,
                                                 const int                            ny   = 1,
                                                 const int                            nz   = 1,
                                                 const geom_t                         xmin = 0,
                                                 const geom_t                         ymin = 0,
                                                 const geom_t                         zmin = 0,
                                                 const geom_t                         xmax = 1,
                                                 const geom_t                         ymax = 1,
                                                 const geom_t                         zmax = 1);

        static std::shared_ptr<Mesh> create_hex8_cube(const std::shared_ptr<Communicator> &comm,
                                                      const int                            nx   = 1,
                                                      const int                            ny   = 1,
                                                      const int                            nz   = 1,
                                                      const geom_t                         xmin = 0,
                                                      const geom_t                         ymin = 0,
                                                      const geom_t                         zmin = 0,
                                                      const geom_t                         xmax = 1,
                                                      const geom_t                         ymax = 1,
                                                      const geom_t                         zmax = 1);

        static std::shared_ptr<Mesh> create_tet4_cube(const std::shared_ptr<Communicator> &comm,
                                                      const int                            nx   = 1,
                                                      const int                            ny   = 1,
                                                      const int                            nz   = 1,
                                                      const geom_t                         xmin = 0,
                                                      const geom_t                         ymin = 0,
                                                      const geom_t                         zmin = 0,
                                                      const geom_t                         xmax = 1,
                                                      const geom_t                         ymax = 1,
                                                      const geom_t                         zmax = 1);

        static std::shared_ptr<Mesh> create_square(const std::shared_ptr<Communicator> &comm,
                                                   const enum ElemType                  element_type,
                                                   const int                            nx   = 1,
                                                   const int                            ny   = 1,
                                                   const geom_t                         xmin = 0,
                                                   const geom_t                         ymin = 0,
                                                   const geom_t                         xmax = 1,
                                                   const geom_t                         ymax = 1);

        static std::shared_ptr<Mesh> create_tri3_square(const std::shared_ptr<Communicator> &comm,
                                                        const int                            nx   = 1,
                                                        const int                            ny   = 1,
                                                        const geom_t                         xmin = 0,
                                                        const geom_t                         ymin = 0,
                                                        const geom_t                         xmax = 1,
                                                        const geom_t                         ymax = 1);

        static std::shared_ptr<Mesh> create_quad4_square(const std::shared_ptr<Communicator> &comm,
                                                         const int                            nx   = 1,
                                                         const int                            ny   = 1,
                                                         const geom_t                         xmin = 0,
                                                         const geom_t                         ymin = 0,
                                                         const geom_t                         xmax = 1,
                                                         const geom_t                         ymax = 1);

        static std::shared_ptr<Mesh> create_hex8_checkerboard_cube(const std::shared_ptr<Communicator> &comm,
                                                                   const int                            nx   = 2,
                                                                   const int                            ny   = 2,
                                                                   const int                            nz   = 2,
                                                                   const geom_t                         xmin = 0,
                                                                   const geom_t                         ymin = 0,
                                                                   const geom_t                         zmin = 0,
                                                                   const geom_t                         xmax = 1,
                                                                   const geom_t                         ymax = 1,
                                                                   const geom_t                         zmax = 1);

        static std::shared_ptr<Mesh> create_hex8_bidomain_cube(const std::shared_ptr<Communicator> &comm,
                                                               const int                            nx   = 2,
                                                               const int                            ny   = 2,
                                                               const int                            nz   = 2,
                                                               const geom_t                         xmin = 0,
                                                               const geom_t                         ymin = 0,
                                                               const geom_t                         zmin = 0,
                                                               const geom_t                         xmax = 1,
                                                               const geom_t                         ymax = 1,
                                                               const geom_t                         zmax = 1);

        std::vector<std::pair<block_idx_t, SharedBuffer<element_idx_t>>> select_elements(
                const std::function<bool(const geom_t, const geom_t, const geom_t)> &selector,
                const std::vector<std::string>                                      &block_names = {});

        int  split_block(const SharedBuffer<element_idx_t> &elements, const std::string &name);
        int  split_boundary_layer();
        int  renumber_nodes();
        int  renumber_nodes(const SharedBuffer<idx_t> &node_mapping);
        void set_node_mapping(const SharedBuffer<idx_t> &node_mapping);
        void set_comm(const std::shared_ptr<Communicator> &comm);
        void set_element_type(const enum ElemType element_type);
        std::pair<SharedBuffer<geom_t>, SharedBuffer<geom_t>> compute_bounding_box();

        std::shared_ptr<Mesh> clone() const;

        void reorder_elements_from_tags(const SharedBuffer<idx_t> &tags);

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    using SharedMesh  = std::shared_ptr<Mesh>;
    using SharedBlock = std::shared_ptr<Mesh::Block>;

    std::shared_ptr<Mesh> convert_to(const enum ElemType element_type, const std::shared_ptr<Mesh> &mesh);
    std::shared_ptr<Mesh> promote_to(const enum ElemType element_type, const std::shared_ptr<Mesh> &mesh);
}  // namespace smesh

#endif  // SMESH_MESH_HPP
