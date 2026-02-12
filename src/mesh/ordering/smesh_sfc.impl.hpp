#ifndef SMESH_SFC_IMPL_HPP
#define SMESH_SFC_IMPL_HPP

#include "smesh_common.hpp"
#include "smesh_sfc.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace smesh {

template <typename geom_t>
static inline u32 quantize_to_u32(geom_t v, geom_t vmin, geom_t vmax, u32 Q) {
  const geom_t range = vmax - vmin;
  if (!(range > geom_t(0)))
    return 0u;

  geom_t t = (v - vmin) / range;
  t = clamp(t, geom_t(0), geom_t(1));

  // round-to-nearest; deterministic alternative: floor(t*Q + 0.5)
  const long double qf = (long double)t * (long double)Q;
  long long qi = llround(qf);

  if (qi < 0)
    qi = 0;
  if (qi > (long long)Q)
    qi = (long long)Q;
  return (u32)qi;
}

// Dilate 21-bit integer so bits occupy positions 0,3,6,... in 64-bit word.
static inline u64 part1by2_21(u32 x) {
  u64 v = (u64)(x & 0x1FFFFFu);
  v = (v | (v << 32)) & 0x1F00000000FFFFULL;
  v = (v | (v << 16)) & 0x1F0000FF0000FFULL;
  v = (v | (v << 8)) & 0x100F00F00F00F00FULL;
  v = (v | (v << 4)) & 0x10C30C30C30C30C3ULL;
  v = (v | (v << 2)) & 0x1249249249249249ULL;
  return v;
}

template <typename geom_t>
int encode_morton3(const ptrdiff_t n_points,
                   const geom_t *const SMESH_RESTRICT x,
                   const geom_t *const SMESH_RESTRICT y,
                   const geom_t *const SMESH_RESTRICT z,
                   u32 *const SMESH_RESTRICT encoding) {
  if (n_points <= 0)
    return 0;
  if (!x || !y || !z || !encoding)
    return SMESH_FAILURE;

  constexpr int DIG = std::numeric_limits<u32>::digits;
  constexpr int BITS_PER_AXIS = (DIG / 3 > 21) ? 21 : (DIG / 3);
  constexpr u32 Q = (u32(1) << BITS_PER_AXIS) - 1u;

  geom_t x_min, x_max, y_min, y_max, z_min, z_max;
  minmax(n_points, x, &x_min, &x_max);
  minmax(n_points, y, &y_min, &y_max);
  minmax(n_points, z, &z_min, &z_max);

  const geom_t x_range = x_max - x_min;
  const geom_t y_range = y_max - y_min;
  const geom_t z_range = z_max - z_min;

  const bool x_deg = !(x_range > geom_t(0));
  const bool y_deg = !(y_range > geom_t(0));
  const bool z_deg = !(z_range > geom_t(0));

  for (ptrdiff_t i = 0; i < n_points; ++i) {
    const u32 ix = x_deg ? 0u : quantize_to_u32<geom_t>(x[i], x_min, x_max, Q);
    const u32 iy = y_deg ? 0u : quantize_to_u32<geom_t>(y[i], y_min, y_max, Q);
    const u32 iz = z_deg ? 0u : quantize_to_u32<geom_t>(z[i], z_min, z_max, Q);

    const u32 code =
        part1by2_21(ix) | (part1by2_21(iy) << 1) | (part1by2_21(iz) << 2);

    encoding[i] = (u32)code;
  }

  return SMESH_SUCCESS;
}

// Convert 3D axes -> Hilbert transpose (in-place).
static inline void hilbert3_axes_to_transpose(const unsigned nBits, u32 X[3]) {
  // This is the standard general-nDims algorithm specialized to nDims=3.
  // Preconditions: nBits >= 1, X[d] in [0, 2^nBits-1]
  u32 M = 1u << (nBits - 1);
  u32 P, Q, t;

  // "Inverse undo"
  for (Q = M; Q > 1; Q >>= 1) {
    P = Q - 1;
    // i = 0..2
    for (unsigned i = 0; i < 3; ++i) {
      if (X[i] & Q) {
        X[0] ^= P;
      } else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }

  // Gray encode
  X[1] ^= X[0];
  X[2] ^= X[1];

  // Prefix transformation
  t = 0;
  for (Q = M; Q > 1; Q >>= 1) {
    if (X[2] & Q)
      t ^= (Q - 1);
  }
  X[0] ^= t;
  X[1] ^= t;
  X[2] ^= t;
}

// Pack transpose words into scalar Hilbert index by bit-interleaving
// (MSB-first).
static inline u32 hilbert3_transpose_to_index_u32(const unsigned nBits,
                                                  const u32 X[3]) {
  // nBits must satisfy 3*nBits <= 32 for u32 output.
  u32 idx = 0;
  for (int b = (int)nBits - 1; b >= 0; --b) {
    idx = (idx << 1) | ((X[0] >> b) & 1u);
    idx = (idx << 1) | ((X[1] >> b) & 1u);
    idx = (idx << 1) | ((X[2] >> b) & 1u);
  }
  return idx;
}

static inline u32 hilbert3D_u32(const unsigned nBits, u32 x, u32 y, u32 z) {
  u32 X[3] = {x, y, z};
  hilbert3_axes_to_transpose(nBits, X);
  return hilbert3_transpose_to_index_u32(nBits, X);
}

template <typename geom_t>
int encode_hilbert3(const ptrdiff_t n_points,
                    const geom_t *const SMESH_RESTRICT x,
                    const geom_t *const SMESH_RESTRICT y,
                    const geom_t *const SMESH_RESTRICT z,
                    u32 *const SMESH_RESTRICT encoding) {
  if (n_points <= 0)
    return 0;
  if (!x || !y || !z || !encoding)
    return SMESH_FAILURE;

  // Choose bits-per-axis that *fit* in u32: digits(u32)=32 => floor(32/3)=10
  constexpr int DIG = std::numeric_limits<u32>::digits;
  constexpr unsigned BITS_PER_AXIS =
      (unsigned)((DIG / 3 > 21) ? 21 : (DIG / 3));
  static_assert(BITS_PER_AXIS >= 1, "u32 too small for 3D Hilbert encoding");

  constexpr u32 Q = (u32(1) << BITS_PER_AXIS) - 1u;

  geom_t x_min, x_max, y_min, y_max, z_min, z_max;
  minmax(n_points, x, &x_min, &x_max);
  minmax(n_points, y, &y_min, &y_max);
  minmax(n_points, z, &z_min, &z_max);

  const geom_t x_range = x_max - x_min;
  const geom_t y_range = y_max - y_min;
  const geom_t z_range = z_max - z_min;

  const bool x_deg = !(x_range > geom_t(0));
  const bool y_deg = !(y_range > geom_t(0));
  const bool z_deg = !(z_range > geom_t(0));

  for (ptrdiff_t i = 0; i < n_points; ++i) {
    const u32 ix = x_deg ? 0u : quantize_to_u32<geom_t>(x[i], x_min, x_max, Q);
    const u32 iy = y_deg ? 0u : quantize_to_u32<geom_t>(y[i], y_min, y_max, Q);
    const u32 iz = z_deg ? 0u : quantize_to_u32<geom_t>(z[i], z_min, z_max, Q);

    const u32 code = hilbert3D_u32(BITS_PER_AXIS, ix, iy, iz);
    encoding[i] = code;
  }

  return SMESH_SUCCESS;
}

static inline bool valid_perm3(int a, int b, int c) {
  if ((unsigned)a > 2u || (unsigned)b > 2u || (unsigned)c > 2u)
    return false;
  if (a == b || a == c || b == c)
    return false;
  return true;
}

template <typename geom_t>
int encode_cartesian3(const ptrdiff_t n_points,
                      const geom_t *const SMESH_RESTRICT x,
                      const geom_t *const SMESH_RESTRICT y,
                      const geom_t *const SMESH_RESTRICT z, int fast, int mid,
                      int slow, u32 *const SMESH_RESTRICT encoding) {
  if (n_points <= 0)
    return 0;
  if (!x || !y || !z || !encoding)
    return SMESH_FAILURE;

  auto valid_perm3 = [](int a, int b, int c) {
    if ((unsigned)a > 2u || (unsigned)b > 2u || (unsigned)c > 2u)
      return false;
    if (a == b || a == c || b == c)
      return false;
    return true;
  };

  if (!valid_perm3(fast, mid, slow))
    return SMESH_FAILURE;

  // pick bits-per-axis that fit in u32 code: 3*BITS_PER_AXIS <= 32
  constexpr int DIG = std::numeric_limits<u32>::digits;          // 32
  constexpr int BITS_PER_AXIS = (DIG / 3 > 21) ? 21 : (DIG / 3); // 10 for u32
  static_assert(BITS_PER_AXIS >= 1, "u32 too small for 3D cartesian encoding");

  constexpr u32 Q = (u32(1) << BITS_PER_AXIS) - 1u; // per-axis max

  // bounding box
  geom_t x_min, x_max, y_min, y_max, z_min, z_max;
  minmax(n_points, x, &x_min, &x_max);
  minmax(n_points, y, &y_min, &y_max);
  minmax(n_points, z, &z_min, &z_max);

  const geom_t x_range = x_max - x_min;
  const geom_t y_range = y_max - y_min;
  const geom_t z_range = z_max - z_min;

  const bool x_deg = !(x_range > geom_t(0));
  const bool y_deg = !(y_range > geom_t(0));
  const bool z_deg = !(z_range > geom_t(0));

  // order[0]=fastest, [1]=mid, [2]=slowest
  const int order[3] = {fast, mid, slow};

  for (ptrdiff_t i = 0; i < n_points; ++i) {
    const u32 ix = x_deg ? 0u : quantize_to_u32<geom_t>(x[i], x_min, x_max, Q);
    const u32 iy = y_deg ? 0u : quantize_to_u32<geom_t>(y[i], y_min, y_max, Q);
    const u32 iz = z_deg ? 0u : quantize_to_u32<geom_t>(z[i], z_min, z_max, Q);

    // gather into idx[dim]
    const u32 idx[3] = {ix, iy, iz};

    // encode: order[0] fastest => code = a + (Q+1)*(b + (Q+1)*c)
    const u32 base = Q + 1u;

    const u32 a = idx[order[0]];
    const u32 b = idx[order[1]];
    const u32 c = idx[order[2]];

    // 3*BITS_PER_AXIS <= 32 ensures this fits in u32
    const u32 code = a + base * (b + base * c);

    encoding[i] = code;
  }

  return SMESH_SUCCESS;
}

} // namespace smesh

#endif // SMESH_SFC_IMPL_HPP
