#include "smesh_context.hpp"
#include "smesh_output.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("io_test.exe");

  auto ctx = smesh::initialize(argc, argv);

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input_folder> <output_folder>\n", argv[0]);
    return SMESH_FAILURE;
  }

  {
    auto comm = ctx->communicator();

    comm->barrier();
    SMESH_TRACE_SCOPE("body");

    auto mesh = Mesh::create_from_file(comm, Path(argv[1]));
    auto output = Output::create(mesh, Path(argv[2]));

    if (comm->size() == 1) {
      auto geom_type = TypeToEnum<geom_t>::value();
      output->write_nodal("nodal_data", geom_type, mesh->points()->data()[0],
                          1);

      auto idx_type = TypeToEnum<idx_t>::value();
      output->write_elemental("elemental_data", idx_type,
                              mesh->elements()->data()[0], 1);
    } else {
      auto dist = mesh->distributed();
      auto node_mapping = mesh->distributed()->node_mapping();
      auto element_mapping = mesh->distributed()->element_mapping();
      comm->print_callback([&](std::ostream &os) {
        os << "n_nodes_global: " << dist->n_nodes_global() << "\n";
        os << "n_nodes_local: " << dist->n_nodes_local() << "\n";
        os << "n_nodes_owned: " << dist->n_nodes_owned() << "\n";
        os << "n_nodes_shared: " << dist->n_nodes_shared() << "\n";
        os << "n_nodes_ghosts: " << dist->n_nodes_ghosts() << "\n";

        os << "n_elements_global: " << dist->n_elements_global() << "\n";
        os << "n_elements_local: " << dist->n_elements_local() << "\n";
        os << "n_elements_owned: " << dist->n_elements_owned() << "\n";
        os << "n_elements_ghosts: " << dist->n_elements_ghosts() << "\n";
        os << "n_elements_shared: " << dist->n_elements_shared() << "\n";

        os << "element_type: " << type_to_string(mesh->element_type()) << "\n";
        os << "elem_num_nodes: " << elem_num_nodes(mesh->element_type())
           << "\n";
      });

      comm->barrier();

      output->write_nodal("nodal_data", node_mapping);
      //   output->write_elemental("elemental_data", element_idx_type,
      //   mesh->element_mapping()->data(), 1);

      comm->barrier();
    }
  }

  return SMESH_SUCCESS;
}