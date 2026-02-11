
#include "smesh_multiblock_graph.impl.hpp"

// TODO: explicit instantiation of functions in smesh_multiblock_graph.impl.hpp
#define SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(COUNT_T, IDX_T)                \
  template int create_multiblock_n2e<COUNT_T, IDX_T>(                                    \
      const block_idx_t n_blocks, const enum ElemType element_types[],                 \
      const ptrdiff_t n_elements[], idx_t **const elements[],                           \
      const ptrdiff_t n_nodes, block_idx_t **out_block_number, count_t **out_n2eptr,     \
      element_idx_t **out_elindex);

SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(i32, i32);
SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(i64, i32);
SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH(i64, i64);

#undef SMESH_EXPLICIT_INSTANTIATE_MULTIBLOCK_GRAPH