#ifndef SMESH_SORT_HPP
#define SMESH_SORT_HPP

#include <algorithm>
#include <iterator>

#ifdef SFEM_ENABLE_AVX512_SORT
#include "avx512-16bit-qsort.hpp"
#include "avx512-32bit-qsort.hpp"
#include "avx512-64bit-qsort.hpp"
#else
#ifdef SFEM_ENABLE_AVX2_SORT
#include "avx2sort.h"
#endif
#endif

namespace smesh {

template <typename idx_t> size_t sort_and_unique(idx_t *arr, const size_t size) {
#ifdef SFEM_ENABLE_AVX512_SORT
  avx512_qsort<idx_t>(arr, size);
#else
#ifdef SFEM_ENABLE_AVX2_SORT
  avx2::quicksort(arr, size);
#else
  std::sort(arr, arr + size);
#endif
#endif
  auto it = std::unique(arr, arr + size);
  return std::distance(arr, it);
}

} // namespace smesh

#endif // SMESH_SORT_HPP
