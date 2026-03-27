#ifndef SMESH_SEMI_STRUCTURED_MESH_HPP
#define SMESH_SEMI_STRUCTURED_MESH_HPP

// C Includes
#include "smesh_base.hpp"

// C++ Includes
#include "smesh_buffer.hpp"
#include "smesh_crs_graph.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sshex8_graph.hpp"

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

    inline int semistructured_level(const Mesh &mesh) { return semistructured_level(mesh.element_type(0)); }

    inline ptrdiff_t semistructured_interior_start(const Mesh &mesh) {
        const int       level                  = semistructured_level(mesh);
        const ptrdiff_t n_interior_per_element = level > 1 ? static_cast<ptrdiff_t>(level - 1) * (level - 1) * (level - 1) : 0;
        return mesh.n_nodes() - mesh.n_elements() * n_interior_per_element;
    }

    inline std::vector<int> derefinement_levels(const Mesh &mesh) {
        const int        level   = semistructured_level(mesh);
        const int        nlevels = smesh::sshex8_hierarchical_n_levels(level);
        std::vector<int> levels(nlevels);
        smesh::sshex8_hierarchical_mesh_levels(level, nlevels, levels.data());
        return levels;
    }

    inline int semistructured_export_as_standard(const std::shared_ptr<Mesh> &mesh, const char *path) {
        auto standard_mesh = smesh::sshex_to_hex8(mesh);
        if (!standard_mesh) {
            return SMESH_FAILURE;
        }

        return standard_mesh->write(smesh::Path(path));
    }

}  // namespace smesh

#endif  // SMESH_SEMI_STRUCTURED_MESH_HPP
