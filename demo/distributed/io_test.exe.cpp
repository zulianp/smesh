#include "smesh_context.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_output.hpp"
#include "smesh_packed_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <cstring>
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
      output->write_nodal("nodal_data", TypeToEnum<geom_t>::value(),
                          mesh->points()->data()[0], 1);

      output->write_elemental("elemental_data", TypeToEnum<idx_t>::value(),
                              mesh->elements(0)->data()[0], 1);
    } else {
      auto dist = mesh->distributed();
      auto node_mapping = mesh->distributed()->node_mapping();
      auto element_mapping = mesh->distributed()->element_mapping();

      if (0) {
        comm->print_callback([&](std::ostream &os) {
          os << "n_nodes_global: " << dist->n_nodes_global() << "\n";
          os << "n_nodes_local: " << dist->n_nodes_local() << "\n";
          os << "n_nodes_owned: " << dist->n_nodes_owned() << "\n";
          os << "n_nodes_shared: " << dist->n_nodes_shared() << "\n";
          os << "n_nodes_ghosts: " << dist->n_nodes_ghosts() << "\n";
          os << "n_nodes_aura: " << dist->n_nodes_aura() << "\n";

          os << "n_elements_global: " << dist->n_elements_global() << "\n";
          os << "n_elements_local: " << dist->n_elements_local() << "\n";
          os << "n_elements_owned: " << dist->n_elements_owned() << "\n";
          os << "n_elements_ghosts: " << dist->n_elements_ghosts() << "\n";
          os << "n_elements_shared: " << dist->n_elements_shared() << "\n";

          os << "element_type: " << type_to_string(mesh->element_type(0))
             << "\n";
          os << "elem_num_nodes: " << elem_num_nodes(mesh->element_type(0))
             << "\n";
        });
      }

      comm->barrier();

      output->write_nodal("nodal_data", node_mapping);
      output->write_elemental("elemental_data",
                              TypeToEnum<large_idx_t>::value(),
                              dist->element_mapping()->data(), 1);

      ptrdiff_t element_offset = 0;
      const ptrdiff_t n_owned_elements = dist->n_elements_owned();
      SMESH_MPI_CATCH(MPI_Exscan(&n_owned_elements, &element_offset, 1,
                                 mpi_type<ptrdiff_t>(), MPI_SUM, comm->get()));
      if (comm->rank() == 0) {
        element_offset = 0;
      }

      auto element_global_numbering =
          create_host_buffer<large_idx_t>(n_owned_elements);
      auto element_global_numbering_data = element_global_numbering->data();
      for (ptrdiff_t i = 0; i < n_owned_elements; ++i) {
        element_global_numbering_data[i] =
            static_cast<large_idx_t>(element_offset + i);
      }

      output->write_elemental("sorted_ids", TypeToEnum<large_idx_t>::value(),
                              element_global_numbering_data, 1);

      comm->barrier();

      if (0) {
        auto local_path = Path("test_" + std::to_string(comm->rank()));
        auto local_mesh = std::make_shared<Mesh>(
            Communicator::self(), mesh->blocks(), mesh->points());
        local_mesh->write(local_path);

        auto aura_element_mapping = dist->aura_element_mapping();
        auto local_element_mapping =
            create_host_buffer<large_idx_t>(local_mesh->n_elements());
        std::memcpy(local_element_mapping->data(), element_mapping->data(),
                    element_mapping->nbytes());
        std::memcpy(local_element_mapping->data() + element_mapping->size(),
                    aura_element_mapping->data(),
                    aura_element_mapping->nbytes());

        auto local_output = Output::create(local_mesh, local_path);
        local_output->write_elemental("elemental_data", local_element_mapping);
      }
    }

    comm->barrier();
  }

  return SMESH_SUCCESS;
}
