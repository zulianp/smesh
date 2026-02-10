#include "smesh_crs_graph.impl.hpp"

namespace smesh {
    template class CRSGraph<i64, i64>;
    template class CRSGraph<i64, i32>;
    template class CRSGraph<i32, i32>;
} // namespace smesh
