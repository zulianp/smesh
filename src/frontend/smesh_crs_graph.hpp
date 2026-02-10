#ifndef SMESH_CRS_GRAPH_HPP
#define SMESH_CRS_GRAPH_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"

#include <memory>
#include <iostream>

namespace smesh {
    template<typename count_t, typename idx_t>
	class CRSGraph final {
	public:
	    CRSGraph();
	    ~CRSGraph();
	    CRSGraph(const SharedBuffer<count_t>& rowptr, const SharedBuffer<idx_t>& colidx);
	    ptrdiff_t                        n_nodes() const;
	    ptrdiff_t                        nnz() const;
	    SharedBuffer<count_t>            rowptr() const;
	    SharedBuffer<idx_t>              colidx() const;
	    std::shared_ptr<CRSGraph>        block_to_scalar(const int block_size);
	    void print(std::ostream &os = std::cout) const;
	private:
	    class Impl;
	    std::unique_ptr<Impl> impl_;
	};

} // namespace smesh

#endif // SMESH_CRS_GRAPH_HPP
