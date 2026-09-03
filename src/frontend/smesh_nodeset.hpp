#ifndef SMESH_NODESET_HPP
#define SMESH_NODESET_HPP

#include "smesh_edgeset.hpp"

#include <iostream>
#include <memory>
#include <vector>

namespace smesh {

    /// Mesh-global node ids. Optional `node_mapping` is used for MPI write /
    /// redistribute (serial GIDs on disk). Unstructured refine/promote does not
    /// remap nodesets; `Mesh::renumber_nodes` and SFC reorder do.
    class Nodeset final {
    public:
        int read(const std::shared_ptr<Communicator> &comm, const Path &path);

        int read_and_redistibute(const std::shared_ptr<Mesh> &mesh, const Path &path);
        int redistribute(const std::shared_ptr<Mesh> &mesh);

        SharedBuffer<idx_t>             nodes() const;
        SharedBuffer<large_idx_t>       node_mapping() const;
        static std::shared_ptr<Nodeset> create_from_file(const std::shared_ptr<Communicator> &comm,
                                                         const Path                          &path);

        ptrdiff_t                     size() const;
        std::shared_ptr<Communicator> comm() const;
        int                           write(const Path &path) const;

        Nodeset(const std::shared_ptr<Communicator> &comm,
                const SharedBuffer<idx_t>           &nodes,
                const SharedBuffer<large_idx_t>     &node_mapping = nullptr);
        Nodeset();
        ~Nodeset();

        static std::shared_ptr<Nodeset> create(const std::shared_ptr<Communicator> &comm,
                                               const SharedBuffer<idx_t>           &nodes,
                                               const SharedBuffer<large_idx_t>     &node_mapping = nullptr);

        /// `nodes[i] <- old_to_new[nodes[i]]`. `n` is the local node count.
        int remap_nodes(const idx_t *old_to_new, ptrdiff_t n);

        void print(std::ostream &os = std::cout) const;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    std::shared_ptr<Nodeset> create_nodeset_from_edgeset(const std::shared_ptr<Mesh>    &mesh,
                                                         const std::shared_ptr<Edgeset> &edgeset);

}  // namespace smesh

#endif  // SMESH_NODESET_HPP
