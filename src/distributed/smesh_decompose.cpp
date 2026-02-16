#include "smesh_decompose.impl.hpp"

namespace smesh {
template int create_n2e<int, int, int>(MPI_Comm comm, const ptrdiff_t n_local_elements,
               const ptrdiff_t n_global_elements, const ptrdiff_t n_local_nodes,
               const ptrdiff_t n_global_nodes, const int nnodesxelem,
               const int *const SMESH_RESTRICT *const SMESH_RESTRICT elems,
               int **out_n2eptr, int **out_elindex);
}
