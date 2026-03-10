#ifndef SMESH_SEMI_STRUCTURED_MESH_HPP
#define SMESH_SEMI_STRUCTURED_MESH_HPP

// C Includes
#include "smesh_base.hpp"

// C++ Includes
#include "smesh_buffer.hpp"
#include "smesh_crs_graph.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_mesh.hpp"

#include <memory>
#include <vector>

namespace smesh {
    void sshex_block_to_hex8_block(const Mesh::Block &block, Mesh::Block &new_block);
    void ssquad_block_to_quad4_block(const Mesh::Block &block, Mesh::Block &new_block);

    int semistructured_hierarchical_renumbering(const enum ElemType          element_type,
                                                const int                    level,
                                                const ptrdiff_t              n_nodes,
                                                const SharedBuffer<idx_t *> &elements);

    int semistructured_hierarchical_renumbering(const enum ElemType           element_type,
                                                const int                     level,
                                                const ptrdiff_t               n_nodes,
                                                const SharedBuffer<idx_t *>  &elements,
                                                const SharedBuffer<geom_t *> &points);

    std::shared_ptr<Mesh> to_semistructured(const int                    level,
                                            const std::shared_ptr<Mesh> &mesh,
                                            const bool                   hiearchical_ordering = false,
                                            const bool                   use_GLL              = false);
    std::shared_ptr<Mesh> sshex_to_hex8(const std::shared_ptr<Mesh> &sshex);
    std::shared_ptr<Mesh> derefine(const std::shared_ptr<Mesh> &mesh, const int to_level);

    inline int semistructured_level(const Mesh &mesh) { return proteus_hex_micro_elements_per_dim(mesh.element_type(0)); }

}  // namespace smesh

#endif  // SMESH_SEMI_STRUCTURED_MESH_HPP
