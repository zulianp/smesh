#ifndef SMESH_SORT_HPP
#define SMESH_SORT_HPP

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>

#ifdef SMESH_ENABLE_AVX512_SORT
#include "avx512-16bit-qsort.hpp"
#include "avx512-32bit-qsort.hpp"
#include "avx512-64bit-qsort.hpp"
#else
#ifdef SMESH_ENABLE_AVX2_SORT
#include "avx2sort.h"
#endif
#endif

namespace smesh {

/// In-place sort of a primitive array. Uses AVX512/AVX2 qsort when enabled.
template <typename T>
void sort(T *arr, const size_t n) {
    if (n <= 1) {
        return;
    }
#ifdef SMESH_ENABLE_AVX512_SORT
    avx512_qsort<T>(arr, n);
#else
#ifdef SMESH_ENABLE_AVX2_SORT
    avx2::quicksort(arr, n);
#else
    std::sort(arr, arr + n);
#endif
#endif
}

/// In-place sort with a comparator. `std::sort` (AVX kernels are key-only).
template <typename T, typename Cmp>
void sort(T *arr, const size_t n, Cmp cmp) {
    if (n <= 1) {
        return;
    }
    std::sort(arr, arr + n, cmp);
}

template <typename T>
size_t unique(T *arr, const size_t n) {
    if (n <= 1) {
        return n;
    }
    return static_cast<size_t>(std::distance(arr, std::unique(arr, arr + n)));
}

template <typename T, typename Eq>
size_t unique(T *arr, const size_t n, Eq eq) {
    if (n <= 1) {
        return n;
    }
    return static_cast<size_t>(std::distance(arr, std::unique(arr, arr + n, eq)));
}

template <typename idx_t>
size_t sort_and_unique(idx_t *arr, const size_t size) {
    sort(arr, size);
    return unique(arr, size);
}

template <typename K, typename I>
void argsort(const ptrdiff_t n, const K *key, I *idx) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        idx[i] = static_cast<I>(i);
    }
    if (n <= 1) {
        return;
    }
    std::sort(idx, idx + n, [key](const I l, const I r) { return key[l] < key[r]; });
}

template <typename I, typename Cmp, typename = std::enable_if_t<!std::is_pointer<std::decay_t<Cmp>>::value>>
void argsort(const ptrdiff_t n, I *idx, Cmp cmp) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        idx[i] = static_cast<I>(i);
    }
    if (n <= 1) {
        return;
    }
    std::sort(idx, idx + n, cmp);
}

/// Insertion index of `query` in `key[idx[0..n)]` (idx is an argsort of key).
template <typename K, typename I>
ptrdiff_t argsort_lower_bound(const ptrdiff_t n, const K *key, const I *idx, const K query) {
    ptrdiff_t lo = 0;
    ptrdiff_t hi = n;
    while (lo < hi) {
        const ptrdiff_t mid = lo + (hi - lo) / 2;
        if (key[idx[mid]] < query) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

template <typename T>
ptrdiff_t lower_bound(const T *arr, const ptrdiff_t n, const T key) {
    return static_cast<ptrdiff_t>(std::lower_bound(arr, arr + n, key) - arr);
}

}  // namespace smesh

#endif  // SMESH_SORT_HPP
