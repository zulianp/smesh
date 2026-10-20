#ifndef SMESH_EXCHANGE_KERNELS_HPP
#define SMESH_EXCHANGE_KERNELS_HPP

#include "smesh_types.hpp"

#include <cstddef>

namespace smesh {
namespace device {

/// out[i * block_size + c] = in[idx[i] * block_size + c] for i < n, all on the device: a gather's
/// pack. Elements are copied as `element_size` bytes (1, 2, 4 or 8).
int gather_by_index(const ptrdiff_t n, const ptrdiff_t block_size,
                    const size_t element_size, const idx_t *const idx,
                    const void *const in, void *const out);

/// out[idx[i] * block_size + c] = in[i * block_size + c]: a gather's unpack.
int scatter_by_index(const ptrdiff_t n, const ptrdiff_t block_size,
                     const size_t element_size, const idx_t *const idx,
                     const void *const in, void *const out);

/// Waits for the work queued on the device, as before handing its memory to MPI.
int synchronize();

} // namespace device
} // namespace smesh

#endif // SMESH_EXCHANGE_KERNELS_HPP
