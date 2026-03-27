#include "smesh_device_sideset.hpp"
#include "smesh_device_buffer.hpp"

namespace smesh {

    std::shared_ptr<Sideset> to_device(const std::shared_ptr<Sideset> &sideset) {
        if (!sideset) {
            SMESH_ERROR("Input sideset is null");
        }

        return std::make_shared<Sideset>(sideset->comm(), to_device(sideset->parent()), to_device(sideset->lfi()));
    }

    std::vector<std::shared_ptr<Sideset>> to_device(const std::vector<std::shared_ptr<Sideset>> &ss) {
        std::vector<std::shared_ptr<Sideset>> ret;
        for (const auto &s : ss) {
            if (s) ret.push_back(to_device(s));
        }
        return ret;
    }
}  // namespace smesh
