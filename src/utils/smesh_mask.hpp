#ifndef SMESH_MASK_HPP
#define SMESH_MASK_HPP

#include "smesh_base.hpp"
#include "smesh_types.hpp"

namespace smesh {

mask_t *mask_create(ptrdiff_t n);
void mask_destroy(mask_t *ptr);

inline ptrdiff_t mask_count(ptrdiff_t n) {
  return (n + sizeof(mask_t) - 1) / sizeof(mask_t);
}

inline int mask_get(ptrdiff_t i, const mask_t *const mask) {
  const static ptrdiff_t nbits = static_cast<ptrdiff_t>(sizeof(mask_t) * 8);
  const ptrdiff_t idx = i / nbits;
  const ptrdiff_t shift = i - idx * nbits;

  SMESH_ASSERT(shift < nbits);
  SMESH_ASSERT(shift >= 0);

  mask_t m;
#pragma omp atomic read
  m = mask[idx];
  mask_t q = static_cast<mask_t>(1 << shift);
  return !!(m & q);
}

inline void mask_set(ptrdiff_t i, mask_t *const mask) {
  const static ptrdiff_t nbits = static_cast<ptrdiff_t>(sizeof(mask_t) * 8);
  const ptrdiff_t idx = i / nbits;
  const ptrdiff_t shift = i - idx * nbits;
  mask_t q = static_cast<mask_t>(1 << shift);

#pragma omp atomic update
  mask[idx] |= q;
}

inline void mask_unset(ptrdiff_t i, mask_t *const mask) {
  const static ptrdiff_t nbits = static_cast<ptrdiff_t>(sizeof(mask_t) * 8);
  const ptrdiff_t idx = i / nbits;
  const ptrdiff_t shift = i - idx * nbits;
  mask_t q = static_cast<mask_t>(1 << shift);

  {
    mask_t m;
#pragma omp atomic read
    m = mask[idx];

    if (!(q & m)) {
      return;
    }
  }

#pragma omp atomic update
  mask[idx] = q ^ mask[idx];
}

void mask_print(const ptrdiff_t n, mask_t *const mask);

} // namespace smesh

#endif // SMESH_MASK_HPP
