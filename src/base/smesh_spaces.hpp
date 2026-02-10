#ifndef SMESH_SPACES_HPP
#define SMESH_SPACES_HPP

#include "smesh_base.hpp"
#include <string_view>

namespace smesh {
    enum ExecutionSpace { EXECUTION_SPACE_HOST = 0, EXECUTION_SPACE_DEVICE = 1, EXECUTION_SPACE_INVALID = -1 };

    static inline ExecutionSpace execution_space_from_string(const std::string_view &str) {
        if (str == "host") {
            return EXECUTION_SPACE_HOST;
        }

        if (str == "device") {
            return EXECUTION_SPACE_DEVICE;
        }

        SMESH_ERROR("Invalid ExecutionSpace: %s\n", str.data());
        return EXECUTION_SPACE_INVALID;
    }

    enum MemorySpace {
        MEMORY_SPACE_HOST    = EXECUTION_SPACE_HOST,
        MEMORY_SPACE_DEVICE  = EXECUTION_SPACE_DEVICE,
        MEMORY_SPACE_INVALID = EXECUTION_SPACE_INVALID
    };

} // namespace smesh

#endif // SMESH_SPACES_HPP