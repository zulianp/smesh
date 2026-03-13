#ifndef SMESH_DEVICE_SIDESET_HPP
#define SMESH_DEVICE_SIDESET_HPP

#include "smesh_sideset.hpp"
#include <memory>
#include <vector>

namespace smesh {

std::shared_ptr<Sideset> to_device(const std::shared_ptr<Sideset> &ss);

std::vector<std::shared_ptr<Sideset>>
to_device(const std::vector<std::shared_ptr<Sideset>> &ss);

} // namespace smesh

#endif // SMESH_DEVICE_SIDESET_HPP
