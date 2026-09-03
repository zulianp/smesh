#ifndef SMESH_SEMI_STRUCTURED_MESH_HPP
#define SMESH_SEMI_STRUCTURED_MESH_HPP

// C Includes
#include "smesh_base.hpp"

// C++ Includes
#include "smesh_buffer.hpp"
#include "smesh_crs_graph.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_sspyramid.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sswedge.hpp"

#include <memory>
#include <vector>

namespace smesh {
    void sshex_block_to_hex8_block(const Mesh::Block &block, Mesh::Block &new_block);
    void ssquad_block_to_quad4_block(const Mesh::Block &block, Mesh::Block &new_block);
    void sstet_block_to_tet4_block(const Mesh::Block &block, Mesh::Block &new_block);
    void sswedge_block_to_wedge6_block(const Mesh::Block &block, Mesh::Block &new_block);

    int semistructured_hierarchical_renumbering(const enum ElemType          element_type,
                                                const int                    level,
                                                const ptrdiff_t              n_nodes,
                                                const SharedBuffer<idx_t *> &elements,
                                                const bool                   preserve_corner_ordering);

    int semistructured_hierarchical_renumbering(const enum ElemType           element_type,
                                                const int                     level,
                                                const ptrdiff_t               n_nodes,
                                                const SharedBuffer<idx_t *>  &elements,
                                                const SharedBuffer<geom_t *> &points,
                                                const bool                    preserve_corner_ordering);

    std::shared_ptr<Mesh> to_semistructured(const int                    level,
                                            const std::shared_ptr<Mesh> &mesh,
                                            const bool                   hiearchical_ordering = false,
                                            const bool                   use_GLL              = false);

    std::shared_ptr<Mesh> to_semistructured_distributed(const int                    level,
                                                        const std::shared_ptr<Mesh> &mesh,
                                                        const bool                   use_GLL,
                                                        const bool                   hierarchical_ordering = false);
    std::shared_ptr<Mesh> sshex_to_hex8(const std::shared_ptr<Mesh> &sshex);
    std::shared_ptr<Mesh> sstet_to_tet4(const std::shared_ptr<Mesh> &sstet);
    std::shared_ptr<Mesh> ssquad_to_quad4(const std::shared_ptr<Mesh> &ssquad);
    std::shared_ptr<Mesh> sswedge_to_wedge6(const std::shared_ptr<Mesh> &sswedge);
    /// Per-block explode of mixed HEX/QUAD/WEDGE SS to linear types, shared points.
    /// PYRAMID and TET SS are rejected (no pyramid explode; TET would be Kuhn, not Bey).
    std::shared_ptr<Mesh> ss_to_linear(const std::shared_ptr<Mesh> &ss);

    /// Linear type after exploding a QUAD SS block (QUADSHELL* stays QUADSHELL4).
    inline enum ElemType ssquad_linear_type(const enum ElemType ss_type) {
        if (!is_quad_ss_family(ss_type)) {
            return INVALID;
        }
        return (shell_type(ss_type) == ss_type) ? QUADSHELL4 : QUAD4;
    }
    std::shared_ptr<Mesh> derefine(const std::shared_ptr<Mesh> &mesh, const int to_level);

    /// Device SoA view into a finer SSHEX8 connectivity: uploads only a pointer table
    /// (nxe device row pointers), keeping row data owned by \p fine_device_soa.
    SharedBuffer<idx_t *> sshex8_device_elements_view(const SharedBuffer<idx_t *> &fine_device_soa,
                                                      const int                    from_level,
                                                      const int                    to_level);

    /// HEX8 coarsest view of SSHEX8 corners (level-1 conversion): 8 device row pointers
    /// into \p fine_device_soa, matching \c sshex8_to_standard_hex8_mesh ordering.
    SharedBuffer<idx_t *> sshex8_to_hex8_device_elements_view(const SharedBuffer<idx_t *> &fine_device_soa,
                                                              const int                    from_level);

    /// HEX8 / TET4 / QUAD4 / WEDGE6 / PYRAMID5 local node → SS SoA row (VTK/family order).
    inline bool ss_source_family_corners(const enum ElemType family, const int L, int *const corners, int *const n_corners) {
        if (!corners || !n_corners || L < 1) {
            return false;
        }
        if (family == HEX8) {
            corners[0] = sshex8_lidx(L, 0, 0, 0);
            corners[1] = sshex8_lidx(L, L, 0, 0);
            corners[2] = sshex8_lidx(L, L, L, 0);
            corners[3] = sshex8_lidx(L, 0, L, 0);
            corners[4] = sshex8_lidx(L, 0, 0, L);
            corners[5] = sshex8_lidx(L, L, 0, L);
            corners[6] = sshex8_lidx(L, L, L, L);
            corners[7] = sshex8_lidx(L, 0, L, L);
            *n_corners = 8;
            return true;
        }
        if (family == TET4) {
            corners[0] = sstet4_lidx(L, 0, 0, 0);
            corners[1] = sstet4_lidx(L, L, 0, 0);
            corners[2] = sstet4_lidx(L, 0, L, 0);
            corners[3] = sstet4_lidx(L, 0, 0, L);
            *n_corners = 4;
            return true;
        }
        if (family == QUAD4) {
            corners[0] = ssquad4_lidx(L, 0, 0);
            corners[1] = ssquad4_lidx(L, L, 0);
            corners[2] = ssquad4_lidx(L, L, L);
            corners[3] = ssquad4_lidx(L, 0, L);
            *n_corners = 4;
            return true;
        }
        if (family == WEDGE6) {
            corners[0] = sswedge_lidx(L, 0, 0, 0);
            corners[1] = sswedge_lidx(L, L, 0, 0);
            corners[2] = sswedge_lidx(L, 0, L, 0);
            corners[3] = sswedge_lidx(L, 0, 0, L);
            corners[4] = sswedge_lidx(L, L, 0, L);
            corners[5] = sswedge_lidx(L, 0, L, L);
            *n_corners = 6;
            return true;
        }
        if (family == PYRAMID5) {
            corners[0] = sspyramid_lidx(L, 0, 0, 0);
            corners[1] = sspyramid_lidx(L, L, 0, 0);
            corners[2] = sspyramid_lidx(L, L, L, 0);
            corners[3] = sspyramid_lidx(L, 0, L, 0);
            corners[4] = sspyramid_lidx(L, 0, 0, L);
            *n_corners = 5;
            return true;
        }
        return false;
    }

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

    inline int semistructured_export_as_standard(const std::shared_ptr<Mesh> &mesh, const smesh::Path &path) {
        if (!mesh) {
            return SMESH_FAILURE;
        }

        bool all_hex = true;
        bool all_tet = true;
        bool all_quad = true;
        bool all_wedge = true;
        for (size_t b = 0; b < mesh->n_blocks(); ++b) {
            const auto fam = ss_source_family(mesh->element_type(static_cast<block_idx_t>(b)));
            all_hex &= fam == HEX8;
            all_tet &= fam == TET4;
            all_quad &= fam == QUAD4;
            all_wedge &= fam == WEDGE6;
        }

        std::shared_ptr<Mesh> standard;
        if (all_hex) {
            standard = sshex_to_hex8(mesh);
        } else if (all_tet) {
            standard = sstet_to_tet4(mesh);
        } else if (all_quad) {
            standard = ssquad_to_quad4(mesh);
        } else if (all_wedge) {
            standard = sswedge_to_wedge6(mesh);
        } else {
            fprintf(stderr, "semistructured_export_as_standard: mixed SS families are not supported\n");
            return SMESH_FAILURE;
        }
        if (!standard) {
            return SMESH_FAILURE;
        }

        return standard->write(path);
    }

}  // namespace smesh

#endif  // SMESH_SEMI_STRUCTURED_MESH_HPP
