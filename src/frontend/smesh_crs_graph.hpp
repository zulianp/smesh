#ifndef SMESH_CRS_GRAPH_HPP
#define SMESH_CRS_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_device_buffer.hpp"

#include <iostream>
#include <memory>

namespace smesh {
    template <typename count_t, typename idx_t>
    class CRSGraph final {
    public:
        CRSGraph();
        ~CRSGraph();
        CRSGraph(const SharedBuffer<count_t> &rowptr, const SharedBuffer<idx_t> &colidx);
        ptrdiff_t                 n_nodes() const;
        ptrdiff_t                 nnz() const;
        SharedBuffer<count_t>     rowptr() const;
        SharedBuffer<idx_t>       colidx() const;
        std::shared_ptr<CRSGraph> block_to_scalar(const int block_size);
        void                      print(std::ostream &os = std::cout) const;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    template <typename count_t, typename idx_t>
    inline std::shared_ptr<CRSGraph<count_t, idx_t>> to_device(const std::shared_ptr<CRSGraph<count_t, idx_t>> &in) {
        if(!in) {
            SMESH_ERROR("to_device(CRSGraph): grap is null!");
            return nullptr;
        }
        
        if (in->rowptr()->mem_space() == MEMORY_SPACE_DEVICE) {
            return in;
        }

        return std::make_shared<CRSGraph<count_t, idx_t>>(to_device(in->rowptr()), to_device(in->colidx()));
    }

}  // namespace smesh

#endif  // SMESH_CRS_GRAPH_HPP
