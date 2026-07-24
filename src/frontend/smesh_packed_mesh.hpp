#ifndef SMESH_PACKED_HPP
#define SMESH_PACKED_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_mesh.hpp"

#include <memory>

namespace smesh {
    template <typename pack_idx_t>
    class PackedMesh {
    public:
        PackedMesh();
        ~PackedMesh();

        std::shared_ptr<Mesh> mesh() const;
        static std::shared_ptr<PackedMesh> create(const std::shared_ptr<Mesh>    &mesh,
                                              const std::vector<std::string> &block_names = {},
                                              const bool                      modify_mesh = false,
                                             const int pack_size = 0);

        void map_to_packed(const real_t *const SMESH_RESTRICT values,
                           real_t *const SMESH_RESTRICT       out_values,
                           const int                         block_size = 1) const;
        //

        void                   map_to_unpacked(const real_t *const SMESH_RESTRICT values,
                                               real_t *const SMESH_RESTRICT       out_values,
                                               const int                         block_size = 1) const;
        SharedBuffer<geom_t *> points();

        ptrdiff_t                  n_blocks() const;
        std::string                block_name(const int block_idx) const;
        SharedBuffer<pack_idx_t *> elements(const int block_idx) const;
        SharedBuffer<ptrdiff_t>    owned_nodes_ptr(const int block_idx) const;
        SharedBuffer<ptrdiff_t>    n_shared(const int block_idx) const;
        SharedBuffer<ptrdiff_t>    ghost_ptr(const int block_idx) const;
        SharedBuffer<idx_t>        ghost_idx(const int block_idx) const;
        /// CSR gather graph for atomics-free ghost reduction (two-pass apply).
        /// Row r reduces into global dof ghost_reduce_dest()[r] by summing
        /// ghost buffer slots ghost_reduce_idx()[ghost_reduce_ptr()[r] .. ghost_reduce_ptr()[r+1]).
        SharedBuffer<ptrdiff_t>    ghost_reduce_ptr(const int block_idx) const;
        SharedBuffer<ptrdiff_t>    ghost_reduce_idx(const int block_idx) const;
        SharedBuffer<idx_t>        ghost_reduce_dest(const int block_idx) const;
        ptrdiff_t                  n_ghost_entries(const int block_idx) const;
        ptrdiff_t                  n_ghost_reduce_rows(const int block_idx) const;
        ptrdiff_t                  n_packs(const int block_idx) const;
        ptrdiff_t                  n_elements_per_pack(const int block_idx) const;

        ptrdiff_t max_nodes_per_pack() const;

        int write(const Path &path);

        void print(std::ostream &os = std::cout, const int verbosity = 0) const;

    private:
        class Block;
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace smesh

#endif  // SMESH_PACKED_HPP
