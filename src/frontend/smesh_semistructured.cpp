#include "smesh_semistructured.hpp"

// STL
#include <fstream>
#include <sstream>

#include "smesh_glob.hpp"
#include "smesh_line_quadrature.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_sshex8_mesh.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_tracer.hpp"

namespace smesh {

int semistructured_hierarchical_renumbering(
    const enum ElemType element_type, const int level, const ptrdiff_t n_nodes,
    const SharedBuffer<idx_t *> &elements) {
  SMESH_UNUSED(element_type); // FIXME: harcoded for sshex8
  const int nlevels = sshex8_hierarchical_n_levels(level);
  std::vector<int> levels(nlevels);

  // FiXME harcoded for sshex8
  sshex8_hierarchical_mesh_levels(level, nlevels, levels.data());

  return sshex8_hierarchical_renumbering(level, nlevels, levels.data(),
                                         elements->extent(1), n_nodes,
                                         elements->data());
}

// FIXME hardcoded for sshex8
std::shared_ptr<Mesh> to_semistructured(const int level,
                                        const std::shared_ptr<Mesh> &mesh,
                                        const bool hiearchical_ordering,
                                        const bool use_GLL) {
  SMESH_TRACE_SCOPE("to_semistructured");

  if (mesh->n_blocks() == 1) {
    auto block = mesh->block(0);

    auto default_block = std::make_shared<Mesh::Block>();
    default_block->set_name(block->name());

    enum ElemType element_type =
        semistructured_type(block->element_type(), level);
    default_block->set_element_type(element_type);

    const int nxe = sshex8_nxe(level);
    auto elements = create_host_buffer<idx_t>(nxe, mesh->n_elements());

    ptrdiff_t n_unique_nodes{-1};
    ptrdiff_t interior_start{-1};
    sshex8_generate_elements(level, mesh->n_elements(), mesh->n_nodes(),
                             mesh->elements(0)->data(), elements->data(),
                             &n_unique_nodes, &interior_start);

    default_block->set_elements(elements);
    std::vector<std::shared_ptr<Mesh::Block>> blocks;
    blocks.push_back(default_block);

    if (hiearchical_ordering) {
      semistructured_hierarchical_renumbering(element_type, level,
                                              n_unique_nodes, elements);
    }

    auto p = smesh::create_host_buffer<geom_t>(mesh->spatial_dimension(),
                                               n_unique_nodes);
    auto macro_p = mesh->points()->data();

    if (use_GLL) {
      const double *qx{nullptr};
      switch (level) {
      case 1: {
        qx = line_GL_q2_x;
        break;
      }
      case 2: {
        qx = line_GL_q3_x;
        break;
      }
      case 4: {
        qx = line_GL_q5_x;
        break;
      }
      case 8: {
        qx = line_GL_q9_x;
        break;
      }
      case 16: {
        qx = line_GL_q17_x;
        break;
      }
      default: {
        SMESH_ERROR("Unsupported order %d!", level);
      }
      }

      sshex8_fill_points_1D_map(level, mesh->n_elements(), elements->data(),
                                macro_p, qx, p->data());
    } else {
      sshex8_fill_points(level, mesh->n_elements(), elements->data(), macro_p,
                         p->data());
    }

    return std::make_shared<Mesh>(mesh->comm(), blocks, p);

  } else {
    SMESH_ERROR("hex8_to_sshex does not support multi-block meshes");
    return nullptr;
  }
}

void sshex_block_to_hex8_block(const Mesh::Block &block,
                               Mesh::Block &new_block) {
  auto elements = block.elements();
  const int level = proteus_hex_micro_elements_per_dim(block.element_type());
  ptrdiff_t n_micro_elements = block.n_elements() * level * level * level;
  auto hex8_elements = create_host_buffer<idx_t>(8, n_micro_elements);
  sshex8_to_standard_hex8_mesh(level, block.n_elements(), elements->data(),
                               hex8_elements->data());

  new_block.set_name(block.name());
  new_block.set_elements(hex8_elements);
  new_block.set_element_type(HEX8);
}

std::shared_ptr<Mesh> sshex_to_hex8(const std::shared_ptr<Mesh> &sshex) {
  std::vector<std::shared_ptr<Mesh::Block>> blocks;
  for (auto &block : sshex->blocks()) {
    auto new_block = std::make_shared<Mesh::Block>();
    sshex_block_to_hex8_block(*block, *new_block);
    blocks.push_back(new_block);
  }

  return std::make_shared<Mesh>(sshex->comm(), blocks, sshex->points());
}

std::shared_ptr<Mesh> derefine(const std::shared_ptr<Mesh> &mesh,
                               const int to_level) {
  if (mesh->n_blocks() > 1) {
    SMESH_ERROR("derefine is not supported for multiblock meshes");
    return nullptr;
  }

  std::vector<std::shared_ptr<Mesh::Block>> blocks;
  ptrdiff_t n_unique_nodes{-1};
  for (auto &block : mesh->blocks()) {
    
    const int from_level =
        proteus_hex_micro_elements_per_dim(block->element_type());
    const int step_factor = from_level / to_level;
    const int nxe = (to_level + 1) * (to_level + 1) * (to_level + 1);

    auto elements = block->elements();
    auto view = std::make_shared<Buffer<idx_t *>>(
        nxe, block->n_elements(), (idx_t **)malloc(nxe * sizeof(idx_t *)),
        [keep_alive = elements](int, void **v) {
          (void)keep_alive;
          free(v);
        },
        elements->mem_space());

    for (int zi = 0; zi <= to_level; zi++) {
      for (int yi = 0; yi <= to_level; yi++) {
        for (int xi = 0; xi <= to_level; xi++) {
          const int from_lidx = sshex8_lidx(from_level, xi * step_factor,
                                            yi * step_factor, zi * step_factor);
          const int to_lidx = sshex8_lidx(to_level, xi, yi, zi);
          view->data()[to_lidx] = elements->data()[from_lidx];
        }
      }
    }

    enum ElemType element_type =
        semistructured_type(macro_base_elem(block->element_type()), to_level);

    auto derefined_block = std::make_shared<Mesh::Block>();
    derefined_block->set_name(block->name());
    derefined_block->set_element_type(element_type);
    derefined_block->set_elements(view);
    blocks.push_back(derefined_block);

    // Find max node id
    {
      auto vv = view->data();
      const ptrdiff_t nelements = block->n_elements();
      for (size_t v = 0; v < view->extent(0); v++) {
        for (ptrdiff_t e = 0; e < nelements; e++) {
          n_unique_nodes =
              std::max(static_cast<ptrdiff_t>(vv[v][e]), n_unique_nodes);
        }
      }
    }
  }

  n_unique_nodes += 1;

  int sdim = mesh->spatial_dimension();
  auto points = smesh::view(mesh->points(), 0, sdim, 0, n_unique_nodes);

  return std::make_shared<Mesh>(mesh->comm(), blocks, points);
}

} // namespace smesh
