#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_glob.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include <stdio.h>

using namespace smesh;

int main(int argc, char **argv) {
  SMESH_TRACE_SCOPE("label_to_mesh.exe");

  auto ctx = smesh::initialize_serial(argc, argv);
  SMESH_UNUSED(ctx);
  {
    if (argc != 6) {
      fprintf(stderr, "Usage: %s <xc> <yc> <zc> <label.uint8> <output_mesh>\n",
              argv[0]);
      return SMESH_FAILURE;
    }

    const ptrdiff_t nx = std::atoll(argv[1]);
    const ptrdiff_t ny = std::atoll(argv[2]);
    const ptrdiff_t nz = std::atoll(argv[3]);
    const Path label_path = Path(argv[4]);
    const Path output_mesh = Path(argv[5]);

    const u8 label = static_cast<u8>(Env::read<int>("SMESH_LABEL", 1));
    const i64 z_chunk_size = Env::read<i64>("SMESH_Z_CHUNK_SIZE", 128);
    const ptrdiff_t chunk_depth = std::max<i64>(1, z_chunk_size);

    if (nx <= 0 || ny <= 0 || nz <= 0) {
      fprintf(stderr, "Image dimensions must be positive.\n");
      return SMESH_FAILURE;
    }

    const ptrdiff_t nxy = nx * ny;
    const ptrdiff_t nxy_nodes = (nx + 1) * (ny + 1);

    FILE *label_file = fopen(label_path.c_str(), "rb");
    if (!label_file) {
      fprintf(stderr, "Failed to open file %s\n", label_path.c_str());
      return SMESH_FAILURE;
    }

    fseek(label_file, 0, SEEK_END);
    const ptrdiff_t n_labels = ftell(label_file) / (ptrdiff_t)sizeof(u8);
    fseek(label_file, 0, SEEK_SET);
    if (n_labels != nxy * nz) {
      fprintf(stderr,
              "Label array size mismatch: got %ld cells, expected %ld.\n",
              (long)n_labels, (long)(nxy * nz));
      fclose(label_file);
      return SMESH_FAILURE;
    }

    if (create_directory(output_mesh) != SMESH_SUCCESS) {
      fclose(label_file);
      return SMESH_FAILURE;
    }

    std::vector<FILE *> elem_files(8, nullptr);
    for (int v = 0; v < 8; ++v) {
      const Path path =
          output_mesh / ("i" + std::to_string(v) + "." +
                         std::string(TypeToString<idx_t>::value()));
      elem_files[v] = fopen(path.c_str(), "wb");
      if (!elem_files[v]) {
        fprintf(stderr, "Failed to open file %s\n", path.c_str());
        for (int i = 0; i < v; ++i) {
          fclose(elem_files[i]);
        }
        fclose(label_file);
        return SMESH_FAILURE;
      }
    }

    const Path x_out_path =
        output_mesh / ("x." + std::string(TypeToString<geom_t>::value()));
    const Path y_out_path =
        output_mesh / ("y." + std::string(TypeToString<geom_t>::value()));
    const Path z_out_path =
        output_mesh / ("z." + std::string(TypeToString<geom_t>::value()));

    FILE *x_out = fopen(x_out_path.c_str(), "wb");
    FILE *y_out = fopen(y_out_path.c_str(), "wb");
    FILE *z_out = fopen(z_out_path.c_str(), "wb");

    if (!x_out || !y_out || !z_out) {
      if (x_out)
        fclose(x_out);
      if (y_out)
        fclose(y_out);
      if (z_out)
        fclose(z_out);
      for (FILE *f : elem_files) {
        fclose(f);
      }
      fclose(label_file);
      return SMESH_FAILURE;
    }

    const ptrdiff_t ix[8] = {0, 1, 1, 0, 0, 1, 1, 0};
    const ptrdiff_t iy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    const int iz[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    const idx_t invalid = invalid_idx<idx_t>();

    std::vector<idx_t> current_layer(nxy_nodes, invalid);
    std::vector<idx_t> next_layer(nxy_nodes, invalid);
    std::vector<u8> labels(chunk_depth * nxy);
    std::vector<idx_t> elem_buffers[8];
    std::vector<geom_t> x_buffer;
    std::vector<geom_t> y_buffer;
    std::vector<geom_t> z_buffer;

    const ptrdiff_t max_chunk_cells = chunk_depth * nxy;
    for (int v = 0; v < 8; ++v) {
      elem_buffers[v].reserve(max_chunk_cells);
    }
    x_buffer.reserve(chunk_depth * nxy_nodes);
    y_buffer.reserve(chunk_depth * nxy_nodes);
    z_buffer.reserve(chunk_depth * nxy_nodes);

    i64 next_node_id = 0;
    ptrdiff_t n_elements = 0;

    for (ptrdiff_t z_begin = 0; z_begin < nz; z_begin += chunk_depth) {
      const ptrdiff_t chunk_nz = std::min<ptrdiff_t>(chunk_depth, nz - z_begin);
      const size_t chunk_size = (size_t)(chunk_nz * nxy);

      if (fread(labels.data(), sizeof(u8), chunk_size, label_file) !=
          chunk_size) {
        fprintf(stderr, "Failed to read file %s\n", label_path.c_str());
        for (FILE *f : elem_files) {
          fclose(f);
        }
        fclose(x_out);
        fclose(y_out);
        fclose(z_out);
        fclose(label_file);
        return SMESH_FAILURE;
      }

      for (int v = 0; v < 8; ++v) {
        elem_buffers[v].clear();
      }
      x_buffer.clear();
      y_buffer.clear();
      z_buffer.clear();

      for (ptrdiff_t z_local = 0; z_local < chunk_nz; ++z_local) {
        const ptrdiff_t z_idx = z_begin + z_local;
        std::fill(next_layer.begin(), next_layer.end(), invalid);
        const u8 *slice = labels.data() + z_local * nxy;

        for (ptrdiff_t yi = 0; yi < ny; ++yi) {
          const ptrdiff_t row_offset = yi * nx;
          const ptrdiff_t node_row = yi * (nx + 1);
          const ptrdiff_t node_row_top = (yi + 1) * (nx + 1);

          for (ptrdiff_t xi = 0; xi < nx; ++xi) {
            if (slice[row_offset + xi] != label) {
              continue;
            }

            idx_t local_ids[8];
            for (int v = 0; v < 8; ++v) {
              const ptrdiff_t x_node = xi + ix[v];
              const ptrdiff_t y_node = yi + iy[v];
              const ptrdiff_t layer_idx =
                  (y_node == yi ? node_row : node_row_top) + x_node;
              std::vector<idx_t> &layer =
                  iz[v] == 0 ? current_layer : next_layer;
              idx_t &id = layer[layer_idx];

              if (id == invalid) {
                if (next_node_id > std::numeric_limits<idx_t>::max()) {
                  fprintf(stderr, "Node index overflow.\n");
                  for (FILE *f : elem_files) {
                    fclose(f);
                  }
                  fclose(x_out);
                  fclose(y_out);
                  fclose(z_out);
                  fclose(label_file);
                  return SMESH_FAILURE;
                }

                id = static_cast<idx_t>(next_node_id++);
                x_buffer.push_back(static_cast<geom_t>(x_node));
                y_buffer.push_back(static_cast<geom_t>(y_node));
                z_buffer.push_back(static_cast<geom_t>(z_idx + iz[v]));
              }

              local_ids[v] = id;
            }

            for (int v = 0; v < 8; ++v) {
              elem_buffers[v].push_back(local_ids[v]);
            }

            ++n_elements;
          }
        }

        current_layer.swap(next_layer);
      }

      for (int v = 0; v < 8; ++v) {
        if (!elem_buffers[v].empty() &&
            fwrite(elem_buffers[v].data(), sizeof(idx_t),
                   elem_buffers[v].size(),
                   elem_files[v]) != elem_buffers[v].size()) {
          fprintf(stderr, "Failed to write connectivity chunk.\n");
          for (FILE *f : elem_files) {
            fclose(f);
          }
          fclose(x_out);
          fclose(y_out);
          fclose(z_out);
          fclose(label_file);
          return SMESH_FAILURE;
        }
      }

      if ((!x_buffer.empty() &&
           fwrite(x_buffer.data(), sizeof(geom_t), x_buffer.size(), x_out) !=
               x_buffer.size()) ||
          (!y_buffer.empty() &&
           fwrite(y_buffer.data(), sizeof(geom_t), y_buffer.size(), y_out) !=
               y_buffer.size()) ||
          (!z_buffer.empty() &&
           fwrite(z_buffer.data(), sizeof(geom_t), z_buffer.size(), z_out) !=
               z_buffer.size())) {
        fprintf(stderr, "Failed to write point chunk.\n");
        for (FILE *f : elem_files) {
          fclose(f);
        }
        fclose(x_out);
        fclose(y_out);
        fclose(z_out);
        fclose(label_file);
        return SMESH_FAILURE;
      }
    }

    for (FILE *f : elem_files) {
      fclose(f);
    }
    fclose(x_out);
    fclose(y_out);
    fclose(z_out);
    fclose(label_file);

    if (mesh_write_yaml_basic(output_mesh, HEX8, n_elements, 3, next_node_id) !=
        SMESH_SUCCESS) {
      return SMESH_FAILURE;
    }
  }

  return SMESH_SUCCESS;
}
