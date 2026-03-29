#ifndef SMESH_DISTRIBUTED_ARRAY_HPP
#define SMESH_DISTRIBUTED_ARRAY_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_communicator.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_spaces.hpp"

#include "matrixio_array.h"

#include <string_view>

namespace smesh {

template <typename T>
std::shared_ptr<Buffer<T>>
create_buffer_from_file(const std::shared_ptr<Communicator> &comm,
                        const std::string_view &path) {
  T *data{nullptr};
  ptrdiff_t local_size{0}, global_size{0};
  array_create_from_file_convert_from_extension(
      comm->get(), path.data(), (void **)&data, &local_size, &global_size);
  return manage_host_buffer<T>(local_size, data);
}

} // namespace smesh

#endif // SMESH_DISTRIBUTED_ARRAY_HPP