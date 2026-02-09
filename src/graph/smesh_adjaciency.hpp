#ifndef SMESH_ADJACENCY_HPP
#define SMESH_ADJACENCY_HPP

#include "smesh_base.hpp"
#include "smesh_elem_type.hpp"

namespace smesh {

struct LocalSideTable {
  static constexpr int MAX_NUM_SIDES = 6;
  static constexpr int MAX_NUM_NODES_PER_SIDE = 6;
  int nnxs{-1};
  int table[MAX_NUM_SIDES * MAX_NUM_NODES_PER_SIDE];

  inline void fill(enum ElemType element_type) {
    enum ElemType st = side_type(element_type);
    this->nnxs = elem_num_nodes(st);
    if (element_type == TET10 || element_type == TET4) {
      (*this)(0, 0) = 1 - 1;
      (*this)(0, 1) = 2 - 1;
      (*this)(0, 2) = 4 - 1;

      (*this)(1, 0) = 2 - 1;
      (*this)(1, 1) = 3 - 1;
      (*this)(1, 2) = 4 - 1;

      (*this)(2, 0) = 1 - 1;
      (*this)(2, 1) = 4 - 1;
      (*this)(2, 2) = 3 - 1;

      (*this)(3, 0) = 1 - 1;
      (*this)(3, 1) = 3 - 1;
      (*this)(3, 2) = 2 - 1;

      if (element_type == TET10) {
        (*this)(0, 3) = 5 - 1;
        (*this)(0, 4) = 9 - 1;
        (*this)(0, 5) = 8 - 1;

        (*this)(1, 3) = 6 - 1;
        (*this)(1, 4) = 10 - 1;
        (*this)(1, 5) = 9 - 1;

        (*this)(2, 3) = 8 - 1;
        (*this)(2, 4) = 10 - 1;
        (*this)(2, 5) = 7 - 1;

        (*this)(3, 3) = 7 - 1;
        (*this)(3, 4) = 6 - 1;
        (*this)(3, 5) = 5 - 1;
      }

    } else if (element_type == TRI3) {
      (*this)(0, 0) = 1 - 1;
      (*this)(0, 1) = 2 - 1;

      (*this)(1, 0) = 2 - 1;
      (*this)(1, 1) = 3 - 1;

      (*this)(2, 0) = 3 - 1;
      (*this)(2, 1) = 1 - 1;
    } else if (element_type == TRI6) {
      (*this)(0, 0) = 1 - 1;
      (*this)(0, 1) = 2 - 1;
      (*this)(0, 2) = 4 - 1;

      (*this)(1, 0) = 2 - 1;
      (*this)(1, 1) = 3 - 1;
      (*this)(1, 2) = 5 - 1;

      (*this)(2, 0) = 3 - 1;
      (*this)(2, 1) = 1 - 1;
      (*this)(2, 2) = 6 - 1;
    } else if (element_type == QUAD4) {
      (*this)(0, 0) = 1 - 1;
      (*this)(0, 1) = 2 - 1;

      (*this)(1, 0) = 2 - 1;
      (*this)(1, 1) = 3 - 1;

      (*this)(2, 0) = 3 - 1;
      (*this)(2, 1) = 4 - 1;

      (*this)(3, 0) = 4 - 1;
      (*this)(3, 1) = 1 - 1;
    } else if (element_type == HEX8) {
      (*this)(0, 0) = 1 - 1;
      (*this)(0, 1) = 2 - 1;
      (*this)(0, 2) = 6 - 1;
      (*this)(0, 3) = 5 - 1;

      (*this)(1, 0) = 2 - 1;
      (*this)(1, 1) = 3 - 1;
      (*this)(1, 2) = 7 - 1;
      (*this)(1, 3) = 6 - 1;

      (*this)(2, 0) = 3 - 1;
      (*this)(2, 1) = 4 - 1;
      (*this)(2, 2) = 8 - 1;
      (*this)(2, 3) = 7 - 1;

      (*this)(3, 0) = 4 - 1;
      (*this)(3, 1) = 1 - 1;
      (*this)(3, 2) = 5 - 1;
      (*this)(3, 3) = 8 - 1;

      (*this)(4, 0) = 4 - 1;
      (*this)(4, 1) = 3 - 1;
      (*this)(4, 2) = 2 - 1;
      (*this)(4, 3) = 1 - 1;

      (*this)(5, 0) = 5 - 1;
      (*this)(5, 1) = 6 - 1;
      (*this)(5, 2) = 7 - 1;
      (*this)(5, 3) = 8 - 1;
    } else {
      SMESH_ERROR("fill_local_side_table");
    }
  }

  inline int &operator()(const int side, const int node) {
    return table[side * nnxs + node];
  }

  inline int operator()(const int side, const int node) const {
    return table[side * nnxs + node];
  }
};

} // namespace smesh

#endif // SMESH_ADJACENCY_HPP