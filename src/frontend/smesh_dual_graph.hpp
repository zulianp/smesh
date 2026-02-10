#ifndef SMESH_DUAL_GRAPH_HPP
#define SMESH_DUAL_GRAPH_HPP

#include "smesh_base.hpp"

#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"

namespace smesh {

    class DualGraph final{
    public:
        DualGraph();
        ~DualGraph();

        SharedBuffer<count_t> adj_ptr();
        SharedBuffer<element_idx_t> adj_idx();
        static std::shared_ptr<DualGraph> create(const std::shared_ptr<Mesh> &mesh);
    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };
    
} // namespace smesh

#endif // SMESH_DUAL_GRAPH_HPP