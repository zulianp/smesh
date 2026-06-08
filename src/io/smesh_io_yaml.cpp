#include "smesh_io_yaml.hpp"

#ifdef SMESH_ENABLE_RYAML

#include "smesh_alloc.hpp"
#include "smesh_buffer.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_path.hpp"
#include "smesh_read.hpp"

#include <ryml.hpp>

#include <memory>
#include <string>

// TODO: FIXME! Do not assume naming conventions, read the actual nodes and paths from the yaml files.

namespace smesh {
    static std::string yaml_scalar_to_string(const ryml::ConstNodeRef &node) {
        auto v = node.val();
        return std::string(v.str, v.len);
    }

    template <typename T>
    static bool yaml_read(const ryml::ConstNodeRef &node, const char *key, T &value) {
        if (!node.has_child(key)) {
            return false;
        }

        node[key] >> value;
        return true;
    }

    static Path yaml_resolve_path(const ryml::ConstNodeRef &node, const std::string &path) {
        if (path.empty() || path[0] == '/' || (path.size() > 1 && path[1] == ':')) {
            return Path(path);
        }

        if (node.has_child("path")) {
            return Path(yaml_scalar_to_string(node["path"])) / path;
        }

        if (node.has_child("folder")) {
            return Path(yaml_scalar_to_string(node["folder"])) / path;
        }

        return Path(path);
    }

    static bool yaml_sequence_entry_path(const ryml::ConstNodeRef &sequence,
                                         const size_t              index,
                                         const char               *key,
                                         std::string              &path) {
        if (index >= sequence.num_children()) {
            return false;
        }

        auto entry = sequence[index];
        if (!entry.has_child(key)) {
            return false;
        }

        path = yaml_scalar_to_string(entry[key]);
        return true;
    }

    static SharedBuffer<geom_t *> yaml_read_points(const ryml::ConstNodeRef &node) {
        int       spatial_dim = 0;
        ptrdiff_t n_nodes     = 0;
        if (!yaml_read(node, "spatial_dimension", spatial_dim) || !yaml_read(node, "n_nodes", n_nodes) ||
            !node.has_child("points")) {
            SMESH_ERROR("Mesh::create_from_yaml: Missing point metadata\n");
            return nullptr;
        }

        static constexpr char xyz[3][2] = {"x", "y", "z"};

        geom_t **points = (geom_t **)SMESH_CALLOC(spatial_dim, sizeof(geom_t *));
        auto     pnode  = node["points"];
        for (int d = 0; d < spatial_dim; ++d) {
            std::string path;
            if (!yaml_sequence_entry_path(pnode, d, xyz[d], path)) {
                SMESH_ERROR("Mesh::create_from_yaml: Missing point path for %s\n", xyz[d]);
                return nullptr;
            }

            ptrdiff_t n_read = 0;
            if (array_read_convert_from_extension(yaml_resolve_path(node, path), &points[d], &n_read) != SMESH_SUCCESS) {
                for (int i = 0; i < spatial_dim; ++i) {
                    SMESH_FREE(points[i]);
                }
                SMESH_FREE(points);
                return nullptr;
            }

            if (n_read != n_nodes) {
                SMESH_ERROR("Mesh::create_from_yaml: Inconsistent point count %ld != %ld\n", (long)n_read, (long)n_nodes);
                return nullptr;
            }
        }

        return manage_host_buffer<geom_t>(spatial_dim, n_nodes, points);
    }

    static SharedBuffer<idx_t *> yaml_read_elements(const ryml::ConstNodeRef &root,
                                                    const ryml::ConstNodeRef &block,
                                                    const int                 n_nodes_per_element,
                                                    const ptrdiff_t           n_elements) {
        if (!block.has_child("elements")) {
            SMESH_ERROR("Mesh::create_from_yaml: Missing element metadata\n");
            return nullptr;
        }

        idx_t **elements = (idx_t **)SMESH_CALLOC(n_nodes_per_element, sizeof(idx_t *));
        auto    enode    = block["elements"];
        for (int d = 0; d < n_nodes_per_element; ++d) {
            const std::string key = "i" + std::to_string(d);
            std::string       path;
            if (!yaml_sequence_entry_path(enode, d, key.c_str(), path)) {
                SMESH_ERROR("Mesh::create_from_yaml: Missing element path for %s\n", key.c_str());
                return nullptr;
            }

            ptrdiff_t n_read = 0;
            if (array_read_convert_from_extension(yaml_resolve_path(root, path), &elements[d], &n_read) != SMESH_SUCCESS) {
                for (int i = 0; i < n_nodes_per_element; ++i) {
                    SMESH_FREE(elements[i]);
                }
                SMESH_FREE(elements);
                return nullptr;
            }

            if (n_read != n_elements) {
                SMESH_ERROR("Mesh::create_from_yaml: Inconsistent element count %ld != %ld\n", (long)n_read, (long)n_elements);
                return nullptr;
            }
        }

        return manage_host_buffer<idx_t>(n_nodes_per_element, n_elements, elements);
    }

    std::shared_ptr<Mesh> mesh_from_yaml(const std::shared_ptr<Communicator> &comm, const ryml::NodeRef &node) {
        auto ret = std::make_shared<Mesh>(comm);

        ret->set_points(yaml_read_points(node));

        if (node.has_child("blocks")) {
            auto         blocks   = node["blocks"];
            const size_t n_blocks = blocks.num_children();
            for (size_t b = 0; b < n_blocks; ++b) {
                auto blk = blocks[b];

                int       n_nodes_per_element = 0;
                ptrdiff_t n_elements          = 0;
                if (!yaml_read(blk, "elem_num_nodes", n_nodes_per_element) || !yaml_read(blk, "n_elements", n_elements) ||
                    !blk.has_child("name") || !blk.has_child("element_type")) {
                    SMESH_ERROR("Mesh::create_from_yaml: Missing block metadata\n");
                    return nullptr;
                }

                const std::string name         = yaml_scalar_to_string(blk["name"]);
                const std::string element_type = yaml_scalar_to_string(blk["element_type"]);
                ret->add_block(name,
                               type_from_string(element_type.c_str()),
                               yaml_read_elements(node, blk, n_nodes_per_element, n_elements));
            }
        } else {
            int       n_nodes_per_element = 0;
            ptrdiff_t n_elements          = 0;
            if (!yaml_read(node, "elem_num_nodes", n_nodes_per_element) || !yaml_read(node, "n_elements", n_elements) ||
                !node.has_child("element_type")) {
                SMESH_ERROR("Mesh::create_from_yaml: Missing element metadata\n");
                return nullptr;
            }

            const std::string element_type = yaml_scalar_to_string(node["element_type"]);
            ret->add_block("default",
                           type_from_string(element_type.c_str()),
                           yaml_read_elements(node, node, n_nodes_per_element, n_elements));
        }

        return ret;
    }

}  // namespace smesh

#endif  // SMESH_ENABLE_RYAML
