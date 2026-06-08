#ifndef SMESH_FORWARD_DECLARATIONS_HPP
#define SMESH_FORWARD_DECLARATIONS_HPP

#include "smesh_base.hpp"

namespace smesh {
    class Communicator;
    class Mesh;
    class Sideset;
    template <typename T>
    class Buffer;

    template <typename pack_idx_t>
    class PackedMesh;
}  // namespace smesh

#ifdef SMESH_ENABLE_RYAML
namespace c4 {
    namespace yml {
        class NodeRef;
        class ConstNodeRef;
    }  // namespace yml
}  // namespace c4

namespace ryml {
    using c4::yml::ConstNodeRef;
    using c4::yml::NodeRef;
}  // namespace ryml
#endif  // SMESH_ENABLE_RYAML

#endif  // SMESH_FORWARD_DECLARATIONS_HPP
