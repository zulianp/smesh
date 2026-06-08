#ifndef SMESH_IO_YAML_HPP
#define SMESH_IO_YAML_HPP

#include "smesh_mesh.hpp"

#ifdef SMESH_ENABLE_RYAML
namespace smesh {
    std::shared_ptr<Mesh> mesh_from_yaml(const std::shared_ptr<Communicator> &comm, const ryml::NodeRef &node);
}
#endif

#endif  // SMESH_IO_YAML_HPP
