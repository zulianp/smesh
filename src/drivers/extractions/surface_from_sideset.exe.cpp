#include "smesh_buffer.hpp"
#include "smesh_common.hpp"
#include "smesh_context.hpp"
#include "smesh_glob.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_sideset.hpp"
#include "smesh_tracer.hpp"
#include "smesh_semistructured.hpp"

#include "smesh_sidesets.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_ssquad4_mesh.hpp"

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc != 4) {
        fprintf(stderr, "usage: %s <mesh> <sideset> <output_folder>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int SMESH_ELEMENT_REFINE_LEVEL = Env::read("SMESH_ELEMENT_REFINE_LEVEL", 1);
    bool SMESH_EXTRACT_NODESET = Env::read("SMESH_EXTRACT_NODESET", false);
    bool SMESH_CONVERT_TO_STD_MESH = Env::read("SMESH_CONVERT_TO_STD_MESH", true);

    auto        path_mesh    = Path(argv[1]);
    auto        m            = Mesh::create_from_file(ctx->communicator(), path_mesh);
    auto path_sideset = Path(argv[2]);
    auto        s            = Sideset::create_from_file(ctx->communicator(), path_sideset);
    const auto  elements     = m->elements(0)->data();

    // Make sure the folder exists
    create_directory(Path(argv[3]));

    ElemType    element_type       = m->element_type(0);
    std::string path_output_format = Path(argv[3]) / ("i%d." + str(TypeToString<idx_t>::value()));

    if (SMESH_ELEMENT_REFINE_LEVEL <= 1) {
        int  nnxs       = elem_num_nodes(side_type(element_type));
        auto surf_elems = create_host_buffer<idx_t>(nnxs, s->parent()->size());

        {
            SMESH_TRACE_SCOPE("extract_surface_from_sideset");
            if (extract_surface_from_sideset(
                        element_type, elements, s->parent()->size(), s->parent()->data(), s->lfi()->data(), surf_elems->data()) !=
                SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }
        }

        if (surf_elems->to_files(Path(path_output_format)) != SMESH_SUCCESS) {
            SMESH_ERROR("Unable to write files!\n");
        }

        if (SMESH_EXTRACT_NODESET) {
            ptrdiff_t n_nodes{0};
            idx_t    *nodes{nullptr};

            {
                SMESH_TRACE_SCOPE("extract_nodeset_from_sideset");
                if (extract_nodeset_from_sideset(element_type,
                                                 elements,
                                                 s->parent()->size(),
                                                 s->parent()->data(),
                                                 s->lfi()->data(),
                                                 &n_nodes,
                                                 &nodes) != SMESH_SUCCESS) {
                    SMESH_ERROR("Unable to extract nodeset from sideset!\n");
                }
            }

            std::string path_nodes = std::string(argv[3]) + "/nodeset." + std::string(TypeToString<idx_t>::value());
            auto        nodeset    = manage_host_buffer(n_nodes, nodes);
            nodeset->to_file(Path(path_nodes));
        }

    } else {
        if (element_type != HEX8) {
            SMESH_ERROR("Element %s not supported for semi-structured discretization\n", type_to_string(element_type));
        }

        auto      ss    = to_semistructured(SMESH_ELEMENT_REFINE_LEVEL, m, true, false);
        const int level = semistructured_level(*ss);

        std::shared_ptr<Buffer<idx_t *>> surf_elems;

        if (SMESH_CONVERT_TO_STD_MESH) {
            const int nnxs     = 4;
            const int nexs     = level * level;
            surf_elems         = create_host_buffer<idx_t>(nnxs, s->parent()->size() * nexs);
            auto ss_surf_elems = create_host_buffer<idx_t>((level + 1) * (level + 1), s->parent()->size());

            SMESH_TRACE_SCOPE("sshex8_extract_surface_from_sideset");
            if (sshex8_extract_surface_from_sideset(level,
                                                    ss->elements(0)->data(),
                                                    s->parent()->size(),
                                                    s->parent()->data(),
                                                    s->lfi()->data(),
                                                    ss_surf_elems->data()) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }

            ssquad4_to_standard_quad4_mesh(level, s->parent()->size(), ss_surf_elems->data(), surf_elems->data());

        } else {
            int nnxs   = (level + 1) * (level + 1);
            surf_elems = create_host_buffer<idx_t>(nnxs, s->parent()->size());

            SMESH_TRACE_SCOPE("sshex8_extract_surface_from_sideset");
            if (sshex8_extract_surface_from_sideset(level,
                                                    ss->elements(0)->data(),
                                                    s->parent()->size(),
                                                    s->parent()->data(),
                                                    s->lfi()->data(),
                                                    surf_elems->data()) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }
        }

        if (surf_elems->to_files(Path(path_output_format)) != SMESH_SUCCESS) {
            SMESH_ERROR("Unable to write files!\n");
        }

        if (SMESH_EXTRACT_NODESET) {
            ptrdiff_t n_nodes{0};
            idx_t    *nodes{nullptr};

            {
                SMESH_TRACE_SCOPE("sshex8_extract_nodeset_from_sideset");
                if (sshex8_extract_nodeset_from_sideset(level,
                                                        ss->elements(0)->data(),
                                                        s->parent()->size(),
                                                        s->parent()->data(),
                                                        s->lfi()->data(),
                                                        &n_nodes,
                                                        &nodes) != SMESH_SUCCESS) {
                    SMESH_ERROR("Unable to extract nodeset from sideset!\n");
                }
            }

            std::string path_nodes = std::string(argv[3]) + "/nodeset." + std::string(TypeToString<idx_t>::value());
            auto        nodeset    = manage_host_buffer(n_nodes, nodes);
            nodeset->to_file(Path(path_nodes));
        }
    }

    return SMESH_SUCCESS;
}
