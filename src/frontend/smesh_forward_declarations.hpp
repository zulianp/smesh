#ifndef SMESH_FORWARD_DECLARATIONS_HPP
#define SMESH_FORWARD_DECLARATIONS_HPP

#include "smesh_base.hpp"

namespace smesh {
    class Communicator;
    template <typename idx_t, typename geom_t> class Mesh;
    template <typename element_idx_t> class Sideset;
    template <typename T> class Buffer;
} // namespace smesh

#endif // SMESH_FORWARD_DECLARATIONS_HPP
