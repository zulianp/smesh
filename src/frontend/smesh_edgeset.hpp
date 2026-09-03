#ifndef SMESH_EDGESET_HPP
#define SMESH_EDGESET_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"

#include <iostream>
#include <memory>
#include <vector>

namespace smesh {

    /// Lightweight `(block_id, parent[], lei[])` edge set.
    /// `lei` is the 0-based local edge of the parent element. On TRI/QUAD/shell,
    /// local edges coincide with local sides. Semistructured level changes keep
    /// the same `(parent, lei)` valid. Unstructured refine/promote does not remap
    /// edgesets.
    class Edgeset final {
    public:
        int read(const std::shared_ptr<Communicator> &comm, const Path &path, block_idx_t block_id = 0);

        int read_and_redistibute(const std::shared_ptr<Mesh> &mesh, const Path &path, block_idx_t block_id = 0);
        int redistribute(const std::shared_ptr<Mesh> &mesh);

        SharedBuffer<element_idx_t>     parent() const;
        SharedBuffer<i16>               lei() const;
        block_idx_t                     block_id() const;
        SharedBuffer<large_idx_t>       element_mapping() const;
        static std::shared_ptr<Edgeset> create_from_file(const std::shared_ptr<Communicator> &comm,
                                                         const Path                          &path,
                                                         block_idx_t                          block_id = 0);

        ptrdiff_t                     size() const;
        std::shared_ptr<Communicator> comm() const;
        int                           write(const Path &path) const;

        Edgeset(const std::shared_ptr<Communicator> &comm,
                const SharedBuffer<element_idx_t>   &parent,
                const SharedBuffer<i16>             &lei,
                block_idx_t                          block_id        = 0,
                const SharedBuffer<large_idx_t>     &element_mapping = nullptr);
        Edgeset();
        ~Edgeset();

        static std::shared_ptr<Edgeset> create(const std::shared_ptr<Communicator> &comm,
                                               const SharedBuffer<element_idx_t>   &parent,
                                               const SharedBuffer<i16>             &lei,
                                               block_idx_t                          block_id        = 0,
                                               const SharedBuffer<large_idx_t>     &element_mapping = nullptr);

        /// `parent[i] <- old_to_new[parent[i]]`. `n` is the local element count.
        int remap_parents(const element_idx_t *old_to_new, ptrdiff_t n);

        /// `new_to_old[new] = old` (gather permutation from element reorder).
        int apply_element_permutation(const element_idx_t *new_to_old, ptrdiff_t n);

        void print(std::ostream &os = std::cout) const;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    int remap_edgesets(const std::vector<std::shared_ptr<Edgeset>> &edgesets,
                       block_idx_t                                  block_id,
                       const element_idx_t                         *old_to_new,
                       ptrdiff_t                                    n);

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_edges_from_edgeset(
            const std::shared_ptr<Mesh>    &mesh,
            const std::shared_ptr<Edgeset> &edgeset);

}  // namespace smesh

#endif  // SMESH_EDGESET_HPP
